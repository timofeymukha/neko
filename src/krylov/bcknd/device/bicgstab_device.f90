! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Provides a device implementation of the BiCGStab method.
module bicgstab_device
  use num_types, only : rp
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : NEKO_EPS
  use device, only : device_map, device_unmap, device_get_ptr, &
       device_event_create, device_event_destroy, device_event_sync
  use device_math, only : device_rzero, device_copy, device_glsc3, &
       device_add2s1, device_add2s2
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  !> Device implementation of the right-preconditioned BiCGStab method.
  !!
  !! The solver starts from a zero initial guess and uses the preconditioner
  !! provided by [ksp_t](#krylov::ksp_t). The coupled interface solves the
  !! three components independently and does not apply a coupled operator.
  type, public, extends(ksp_t) :: bicgstab_device_t
     !> Search direction \f$p\f$.
     real(kind=rp), allocatable :: p(:)
     !> Preconditioned search direction \f$\hat{p} = M^{-1}p\f$.
     real(kind=rp), allocatable :: p_hat(:)
     !> Residual \f$r\f$.
     real(kind=rp), allocatable :: r(:)
     !> Intermediate residual \f$s = r - \alpha v\f$.
     real(kind=rp), allocatable :: s(:)
     !> Preconditioned intermediate residual \f$\hat{s} = M^{-1}s\f$.
     real(kind=rp), allocatable :: s_hat(:)
     !> Operator action \f$t = A\hat{s}\f$.
     real(kind=rp), allocatable :: t(:)
     !> Operator action \f$v = A\hat{p}\f$.
     real(kind=rp), allocatable :: v(:)
     type(c_ptr) :: p_d = C_NULL_PTR
     type(c_ptr) :: p_hat_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: s_d = C_NULL_PTR
     type(c_ptr) :: s_hat_d = C_NULL_PTR
     type(c_ptr) :: t_d = C_NULL_PTR
     type(c_ptr) :: v_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     !> Initialise a device BiCGStab solver.
     procedure, pass(this) :: init => bicgstab_device_init
     !> Free a device BiCGStab solver.
     procedure, pass(this) :: free => bicgstab_device_free
     !> Solve a linear system with the device BiCGStab method.
     procedure, pass(this) :: solve => bicgstab_device_solve
     !> Solve three independent systems with the device BiCGStab method.
     procedure, pass(this) :: solve_coupled => bicgstab_device_solve_coupled
  end type bicgstab_device_t

contains

  !> Initialise a device BiCGStab solver.
  subroutine bicgstab_device_init(this, n, max_iter, M, rel_tol, abs_tol, &
       monitor)
    class(bicgstab_device_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%p(n), this%p_hat(n), this%r(n), this%s(n), &
         this%s_hat(n), this%t(n), this%v(n))
    call device_map(this%p, this%p_d, n)
    call device_map(this%p_hat, this%p_hat_d, n)
    call device_map(this%r, this%r_d, n)
    call device_map(this%s, this%s_d, n)
    call device_map(this%s_hat, this%s_hat_d, n)
    call device_map(this%t, this%t_d, n)
    call device_map(this%v, this%v_d, n)

    if (present(M)) then
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(monitor) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, monitor = monitor)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else if (present(monitor)) then
       call this%ksp_init(max_iter, monitor = monitor)
    else
       call this%ksp_init(max_iter)
    end if

    call device_event_create(this%gs_event, 2)

  end subroutine bicgstab_device_init

  !> Free a device BiCGStab solver.
  subroutine bicgstab_device_free(this)
    class(bicgstab_device_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%p)) then
       if (c_associated(this%p_d)) call device_unmap(this%p, this%p_d)
       deallocate(this%p)
    end if
    if (allocated(this%p_hat)) then
       if (c_associated(this%p_hat_d)) then
          call device_unmap(this%p_hat, this%p_hat_d)
       end if
       deallocate(this%p_hat)
    end if
    if (allocated(this%r)) then
       if (c_associated(this%r_d)) call device_unmap(this%r, this%r_d)
       deallocate(this%r)
    end if
    if (allocated(this%s)) then
       if (c_associated(this%s_d)) call device_unmap(this%s, this%s_d)
       deallocate(this%s)
    end if
    if (allocated(this%s_hat)) then
       if (c_associated(this%s_hat_d)) then
          call device_unmap(this%s_hat, this%s_hat_d)
       end if
       deallocate(this%s_hat)
    end if
    if (allocated(this%t)) then
       if (c_associated(this%t_d)) call device_unmap(this%t, this%t_d)
       deallocate(this%t)
    end if
    if (allocated(this%v)) then
       if (c_associated(this%v_d)) call device_unmap(this%v, this%v_d)
       deallocate(this%v)
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine bicgstab_device_free

  !> Solve a linear system with the device BiCGStab method.
  function bicgstab_device_solve(this, Ax, x, f, n, coef, blst, gs_h, &
       niter) result(ksp_results)
    class(bicgstab_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t) :: ksp_results
    integer :: iter, max_iter
    real(kind=rp) :: rnorm, rtr, norm_fac, gamma
    real(kind=rp) :: r_norm, s_norm, shadow_norm, t_norm, v_norm
    real(kind=rp) :: sts, ftv, vtv, stt, ttt
    real(kind=rp) :: beta, alpha, omega, rho_1, rho_2
    type(c_ptr) :: f_d

    f_d = device_get_ptr(f)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(p => this%p, p_hat => this%p_hat, r => this%r, s => this%s, &
         s_hat => this%s_hat, t => this%t, v => this%v, &
         p_d => this%p_d, p_hat_d => this%p_hat_d, r_d => this%r_d, &
         s_d => this%s_d, s_hat_d => this%s_hat_d, t_d => this%t_d, &
         v_d => this%v_d)

      call device_rzero(x%x_d, n)
      call device_copy(r_d, f_d, n)

      rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
      r_norm = bicgstab_device_sqrt(rtr, 'initial residual')
      shadow_norm = r_norm
      rnorm = r_norm * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      if (r_norm .le. 0.0_rp .or. rnorm .lt. this%abs_tol .or. &
           rnorm .lt. gamma) then
         ksp_results%converged = .true.
         return
      end if

      call this%monitor_start('BiCGStab')
      do iter = 1, max_iter
         rho_1 = device_glsc3(f_d, coef%mult_d, r_d, n)
         call bicgstab_device_check_inner_product(rho_1, shadow_norm, &
              r_norm, 'rho inner product')

         if (iter .eq. 1) then
            call device_copy(p_d, r_d, n)
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            if (.not. ieee_is_finite(beta)) then
               call neko_error('BiCGStab failure: non-finite beta')
            end if
            call device_add2s2(p_d, v_d, -omega, n)
            call device_add2s1(p_d, r_d, beta, n)
         end if

         call this%M%solve(p_hat, p, n)
         call Ax%compute(v, p_hat, coef, x%msh, x%Xh)
         call gs_h%op(v, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call blst%apply(v, n)

         ftv = device_glsc3(f_d, coef%mult_d, v_d, n)
         vtv = device_glsc3(v_d, coef%mult_d, v_d, n)
         v_norm = bicgstab_device_sqrt(vtv, 'operator result v')
         call bicgstab_device_check_inner_product(ftv, shadow_norm, v_norm, &
              'alpha denominator')
         alpha = rho_1 / ftv
         if (.not. ieee_is_finite(alpha)) then
            call neko_error('BiCGStab failure: non-finite alpha')
         end if

         call device_copy(s_d, r_d, n)
         call device_add2s2(s_d, v_d, -alpha, n)
         sts = device_glsc3(s_d, coef%mult_d, s_d, n)
         s_norm = bicgstab_device_sqrt(sts, 'intermediate residual')
         rnorm = s_norm * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call device_add2s2(x%x_d, p_hat_d, alpha, n)
            call this%monitor_iter(iter, rnorm)
            exit
         end if

         call this%M%solve(s_hat, s, n)
         call Ax%compute(t, s_hat, coef, x%msh, x%Xh)
         call gs_h%op(t, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call blst%apply(t, n)

         stt = device_glsc3(s_d, coef%mult_d, t_d, n)
         ttt = device_glsc3(t_d, coef%mult_d, t_d, n)
         t_norm = bicgstab_device_sqrt(ttt, 'operator result t')
         if (t_norm .le. 0.0_rp) then
            call neko_error('BiCGStab breakdown: zero omega denominator')
         end if
         if (.not. ieee_is_finite(stt)) then
            call neko_error('BiCGStab failure: non-finite omega numerator')
         end if
         omega = stt / ttt
         if (.not. ieee_is_finite(omega)) then
            call neko_error('BiCGStab failure: non-finite omega')
         end if

         call device_add2s2(x%x_d, p_hat_d, alpha, n)
         call device_add2s2(x%x_d, s_hat_d, omega, n)
         call device_copy(r_d, s_d, n)
         call device_add2s2(r_d, t_d, -omega, n)

         rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
         r_norm = bicgstab_device_sqrt(rtr, 'recursive residual')
         rnorm = r_norm * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         call bicgstab_device_check_inner_product(stt, t_norm, s_norm, &
              'omega numerator')
         rho_2 = rho_1
      end do
    end associate

    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)

  end function bicgstab_device_solve

  !> Check an inner product for a BiCGStab breakdown.
  subroutine bicgstab_device_check_inner_product(inner_product, norm_a, &
       norm_b, quantity)
    real(kind=rp), intent(in) :: inner_product
    real(kind=rp), intent(in) :: norm_a
    real(kind=rp), intent(in) :: norm_b
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: large_norm, small_norm

    if (.not. ieee_is_finite(inner_product)) then
       call neko_error('BiCGStab failure: non-finite ' // trim(quantity))
    end if

    large_norm = max(norm_a, norm_b)
    small_norm = min(norm_a, norm_b)
    if (large_norm .le. 0.0_rp .or. &
         abs(inner_product) / large_norm .le. NEKO_EPS * small_norm) then
       call neko_error('BiCGStab breakdown: near-zero ' // trim(quantity))
    end if

  end subroutine bicgstab_device_check_inner_product

  !> Return the square root of a valid squared norm.
  function bicgstab_device_sqrt(value, quantity) result(root)
    real(kind=rp), intent(in) :: value
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: root

    if (.not. ieee_is_finite(value) .or. value .lt. 0.0_rp) then
       call neko_error('BiCGStab failure: invalid ' // trim(quantity) // &
            ' norm')
    end if
    root = sqrt(value)

  end function bicgstab_device_sqrt

  !> Solve three independent systems with the device BiCGStab method.
  function bicgstab_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(bicgstab_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t), dimension(3) :: ksp_results

    ksp_results(1) = this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function bicgstab_device_solve_coupled

end module bicgstab_device
