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
!> Provides a component-fused coupled device implementation of BiCGStab.
module bicgstab_cpld_device
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
  use device_math, only : device_rzero, device_copy, device_vdot3, &
       device_glsc2, device_add2s1
  use device_mathops, only : device_opadd2cm
  use operators, only : rotate_cyc
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter :: gdim = 3

  !> Component-fused coupled right-preconditioned device BiCGStab method.
  !!
  !! All three components share one Krylov recurrence and all inner products
  !! are combined over the complete vector system. Three-component device
  !! operations are used for products and updates where the portable device
  !! interface provides them.
  type, public, extends(ksp_t) :: bicgstab_cpld_device_t
     !> Three-component search direction \f$p\f$.
     real(kind=rp), allocatable :: p(:, :)
     !> Preconditioned search direction \f$\hat{p} = M^{-1}p\f$.
     real(kind=rp), allocatable :: p_hat(:, :)
     !> Three-component residual \f$r\f$.
     real(kind=rp), allocatable :: r(:, :)
     !> Intermediate residual \f$s = r - \alpha v\f$.
     real(kind=rp), allocatable :: s(:, :)
     !> Preconditioned intermediate residual \f$\hat{s} = M^{-1}s\f$.
     real(kind=rp), allocatable :: s_hat(:, :)
     !> Coupled operator action \f$t = A\hat{s}\f$.
     real(kind=rp), allocatable :: t(:, :)
     !> Coupled operator action \f$v = A\hat{p}\f$.
     real(kind=rp), allocatable :: v(:, :)
     !> Pointwise contributions for combined three-component products.
     real(kind=rp), allocatable :: tmp(:)
     type(c_ptr) :: p_d(gdim) = C_NULL_PTR
     type(c_ptr) :: p_hat_d(gdim) = C_NULL_PTR
     type(c_ptr) :: r_d(gdim) = C_NULL_PTR
     type(c_ptr) :: s_d(gdim) = C_NULL_PTR
     type(c_ptr) :: s_hat_d(gdim) = C_NULL_PTR
     type(c_ptr) :: t_d(gdim) = C_NULL_PTR
     type(c_ptr) :: v_d(gdim) = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     !> Initialise a coupled device BiCGStab solver.
     procedure, pass(this) :: init => bicgstab_cpld_device_init
     !> Free a coupled device BiCGStab solver.
     procedure, pass(this) :: free => bicgstab_cpld_device_free
     !> Reject scalar solves.
     procedure, pass(this) :: solve => bicgstab_cpld_device_solve_scalar
     !> Solve a three-component coupled system.
     procedure, pass(this) :: solve_coupled => bicgstab_cpld_device_solve
  end type bicgstab_cpld_device_t

contains

  !> Initialise a coupled device BiCGStab solver.
  subroutine bicgstab_cpld_device_init(this, n, max_iter, M, rel_tol, &
       abs_tol, monitor)
    class(bicgstab_cpld_device_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor
    integer :: i

    call this%free()

    allocate(this%p(n, gdim), this%p_hat(n, gdim), this%r(n, gdim), &
         this%s(n, gdim), this%s_hat(n, gdim), this%t(n, gdim), &
         this%v(n, gdim), this%tmp(n))

    do i = 1, gdim
       call device_map(this%p(:, i), this%p_d(i), n)
       call device_map(this%p_hat(:, i), this%p_hat_d(i), n)
       call device_map(this%r(:, i), this%r_d(i), n)
       call device_map(this%s(:, i), this%s_d(i), n)
       call device_map(this%s_hat(:, i), this%s_hat_d(i), n)
       call device_map(this%t(:, i), this%t_d(i), n)
       call device_map(this%v(:, i), this%v_d(i), n)
    end do
    call device_map(this%tmp, this%tmp_d, n)

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

  end subroutine bicgstab_cpld_device_init

  !> Free a coupled device BiCGStab solver.
  subroutine bicgstab_cpld_device_free(this)
    class(bicgstab_cpld_device_t), intent(inout) :: this
    integer :: i

    call this%ksp_free()

    do i = 1, gdim
       if (allocated(this%p)) then
          if (c_associated(this%p_d(i))) then
             call device_unmap(this%p(:, i), this%p_d(i))
          end if
       end if
       if (allocated(this%p_hat)) then
          if (c_associated(this%p_hat_d(i))) then
             call device_unmap(this%p_hat(:, i), this%p_hat_d(i))
          end if
       end if
       if (allocated(this%r)) then
          if (c_associated(this%r_d(i))) then
             call device_unmap(this%r(:, i), this%r_d(i))
          end if
       end if
       if (allocated(this%s)) then
          if (c_associated(this%s_d(i))) then
             call device_unmap(this%s(:, i), this%s_d(i))
          end if
       end if
       if (allocated(this%s_hat)) then
          if (c_associated(this%s_hat_d(i))) then
             call device_unmap(this%s_hat(:, i), this%s_hat_d(i))
          end if
       end if
       if (allocated(this%t)) then
          if (c_associated(this%t_d(i))) then
             call device_unmap(this%t(:, i), this%t_d(i))
          end if
       end if
       if (allocated(this%v)) then
          if (c_associated(this%v_d(i))) then
             call device_unmap(this%v(:, i), this%v_d(i))
          end if
       end if
    end do

    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%p_hat)) deallocate(this%p_hat)
    if (allocated(this%r)) deallocate(this%r)
    if (allocated(this%s)) deallocate(this%s)
    if (allocated(this%s_hat)) deallocate(this%s_hat)
    if (allocated(this%t)) deallocate(this%t)
    if (allocated(this%v)) deallocate(this%v)
    if (allocated(this%tmp)) then
       if (c_associated(this%tmp_d)) then
          call device_unmap(this%tmp, this%tmp_d)
       end if
       deallocate(this%tmp)
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine bicgstab_cpld_device_free

  !> Reject a scalar solve with a coupled device BiCGStab solver.
  function bicgstab_cpld_device_solve_scalar(this, Ax, x, f, n, coef, &
       blst, gs_h, niter) result(ksp_results)
    class(bicgstab_cpld_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t) :: ksp_results

    call neko_error('Coupled BiCGStab is only defined for coupled solves')

    ksp_results%res_start = 0.0_rp
    ksp_results%res_final = 0.0_rp
    ksp_results%iter = 0
    ksp_results%converged = .false.

  end function bicgstab_cpld_device_solve_scalar

  !> Solve a three-component coupled system with device BiCGStab.
  function bicgstab_cpld_device_solve(this, Ax, x, y, z, fx, fy, fz, n, &
       coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(bicgstab_cpld_device_t), intent(inout) :: this
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
    type(ksp_monitor_t), dimension(gdim) :: ksp_results
    integer :: iter, max_iter
    real(kind=rp) :: alpha, beta, omega, rho_1, rho_2
    real(kind=rp) :: rnorm, norm_fac, gamma
    real(kind=rp) :: r_norm, s_norm, shadow_norm, t_norm, v_norm
    real(kind=rp) :: rtr, sts, ftv, vtv, stt, ttt
    type(c_ptr) :: f_d(gdim)

    f_d(1) = device_get_ptr(fx)
    f_d(2) = device_get_ptr(fy)
    f_d(3) = device_get_ptr(fz)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(p => this%p, p_hat => this%p_hat, r => this%r, &
         s => this%s, s_hat => this%s_hat, t => this%t, v => this%v, &
         p_d => this%p_d, p_hat_d => this%p_hat_d, r_d => this%r_d, &
         s_d => this%s_d, s_hat_d => this%s_hat_d, t_d => this%t_d, &
         v_d => this%v_d)

      call device_rzero(x%x_d, n)
      call device_rzero(y%x_d, n)
      call device_rzero(z%x_d, n)
      call bicgstab_cpld_device_copy(r_d, f_d, n)

      rtr = bicgstab_cpld_device_product(this%tmp_d, r_d, r_d, &
           coef%mult_d, n)
      r_norm = bicgstab_cpld_device_sqrt(rtr, 'initial residual')
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

      call this%monitor_start('Coupled device BiCGStab')
      do iter = 1, max_iter
         rho_1 = bicgstab_cpld_device_product(this%tmp_d, f_d, r_d, &
              coef%mult_d, n)
         call bicgstab_cpld_device_check_inner_product(rho_1, shadow_norm, &
              r_norm, 'rho inner product')

         if (iter .eq. 1) then
            call bicgstab_cpld_device_copy(p_d, r_d, n)
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            if (.not. ieee_is_finite(beta)) then
               call neko_error(&
                    'Coupled BiCGStab failure: non-finite beta')
            end if
            call device_opadd2cm(p_d(1), p_d(2), p_d(3), v_d(1), &
                 v_d(2), v_d(3), -omega, n, gdim)
            call device_add2s1(p_d(1), r_d(1), beta, n)
            call device_add2s1(p_d(2), r_d(2), beta, n)
            call device_add2s1(p_d(3), r_d(3), beta, n)
         end if

         call this%M%solve(p_hat(:, 1), p(:, 1), n)
         call this%M%solve(p_hat(:, 2), p(:, 2), n)
         call this%M%solve(p_hat(:, 3), p(:, 3), n)

         call Ax%compute_vector(v(:, 1), v(:, 2), v(:, 3), p_hat(:, 1), &
              p_hat(:, 2), p_hat(:, 3), coef, x%msh, x%Xh)
         call bicgstab_cpld_device_assemble(v, v_d, n, coef, blstx, &
              blsty, blstz, gs_h, this%gs_event)

         ftv = bicgstab_cpld_device_product(this%tmp_d, f_d, v_d, &
              coef%mult_d, n)
         vtv = bicgstab_cpld_device_product(this%tmp_d, v_d, v_d, &
              coef%mult_d, n)
         v_norm = bicgstab_cpld_device_sqrt(vtv, 'operator result v')
         call bicgstab_cpld_device_check_inner_product(ftv, shadow_norm, &
              v_norm, 'alpha denominator')
         alpha = rho_1 / ftv
         if (.not. ieee_is_finite(alpha)) then
            call neko_error('Coupled BiCGStab failure: non-finite alpha')
         end if

         call bicgstab_cpld_device_copy(s_d, r_d, n)
         call device_opadd2cm(s_d(1), s_d(2), s_d(3), v_d(1), v_d(2), &
              v_d(3), -alpha, n, gdim)
         sts = bicgstab_cpld_device_product(this%tmp_d, s_d, s_d, &
              coef%mult_d, n)
         s_norm = bicgstab_cpld_device_sqrt(sts, 'intermediate residual')
         rnorm = s_norm * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call device_opadd2cm(x%x_d, y%x_d, z%x_d, p_hat_d(1), &
                 p_hat_d(2), p_hat_d(3), alpha, n, gdim)
            call this%monitor_iter(iter, rnorm)
            exit
         end if

         call this%M%solve(s_hat(:, 1), s(:, 1), n)
         call this%M%solve(s_hat(:, 2), s(:, 2), n)
         call this%M%solve(s_hat(:, 3), s(:, 3), n)

         call Ax%compute_vector(t(:, 1), t(:, 2), t(:, 3), s_hat(:, 1), &
              s_hat(:, 2), s_hat(:, 3), coef, x%msh, x%Xh)
         call bicgstab_cpld_device_assemble(t, t_d, n, coef, blstx, &
              blsty, blstz, gs_h, this%gs_event)

         stt = bicgstab_cpld_device_product(this%tmp_d, s_d, t_d, &
              coef%mult_d, n)
         ttt = bicgstab_cpld_device_product(this%tmp_d, t_d, t_d, &
              coef%mult_d, n)
         t_norm = bicgstab_cpld_device_sqrt(ttt, 'operator result t')
         if (t_norm .le. 0.0_rp) then
            call neko_error(&
                 'Coupled BiCGStab breakdown: zero omega denominator')
         end if
         if (.not. ieee_is_finite(stt)) then
            call neko_error(&
                 'Coupled BiCGStab failure: non-finite omega numerator')
         end if
         omega = stt / ttt
         if (.not. ieee_is_finite(omega)) then
            call neko_error('Coupled BiCGStab failure: non-finite omega')
         end if

         call device_opadd2cm(x%x_d, y%x_d, z%x_d, p_hat_d(1), &
              p_hat_d(2), p_hat_d(3), alpha, n, gdim)
         call device_opadd2cm(x%x_d, y%x_d, z%x_d, s_hat_d(1), &
              s_hat_d(2), s_hat_d(3), omega, n, gdim)
         call bicgstab_cpld_device_copy(r_d, s_d, n)
         call device_opadd2cm(r_d(1), r_d(2), r_d(3), t_d(1), t_d(2), &
              t_d(3), -omega, n, gdim)

         rtr = bicgstab_cpld_device_product(this%tmp_d, r_d, r_d, &
              coef%mult_d, n)
         r_norm = bicgstab_cpld_device_sqrt(rtr, 'recursive residual')
         rnorm = r_norm * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         call bicgstab_cpld_device_check_inner_product(stt, s_norm, t_norm, &
              'omega numerator')
         rho_2 = rho_1
      end do
    end associate

    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)

  end function bicgstab_cpld_device_solve

  !> Copy all three components of a device vector.
  subroutine bicgstab_cpld_device_copy(a_d, b_d, n)
    type(c_ptr), dimension(gdim), intent(inout) :: a_d
    type(c_ptr), dimension(gdim), intent(in) :: b_d
    integer, intent(in) :: n

    call device_copy(a_d(1), b_d(1), n)
    call device_copy(a_d(2), b_d(2), n)
    call device_copy(a_d(3), b_d(3), n)

  end subroutine bicgstab_cpld_device_copy

  !> Compute a combined weighted inner product over three components.
  function bicgstab_cpld_device_product(tmp_d, a_d, b_d, mult_d, n) &
       result(product)
    type(c_ptr), intent(in), value :: tmp_d
    type(c_ptr), dimension(gdim), intent(in) :: a_d
    type(c_ptr), dimension(gdim), intent(in) :: b_d
    type(c_ptr), intent(in), value :: mult_d
    integer, intent(in) :: n
    real(kind=rp) :: product

    call device_vdot3(tmp_d, a_d(1), a_d(2), a_d(3), b_d(1), b_d(2), &
         b_d(3), n)
    product = device_glsc2(tmp_d, mult_d, n)

  end function bicgstab_cpld_device_product

  !> Assemble a coupled operator result and apply component boundary data.
  subroutine bicgstab_cpld_device_assemble(a, a_d, n, coef, blstx, blsty, &
       blstz, gs_h, gs_event)
    integer, intent(in) :: n
    real(kind=rp), dimension(n, gdim), intent(inout) :: a
    type(c_ptr), dimension(gdim), intent(in) :: a_d
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(c_ptr), intent(inout) :: gs_event
    type(c_ptr) :: a1_d, a2_d, a3_d

    a1_d = a_d(1)
    a2_d = a_d(2)
    a3_d = a_d(3)

    call rotate_cyc(a1_d, a2_d, a3_d, 1, coef)
    call gs_h%op(a(:, 1), a(:, 2), a(:, 3), n, GS_OP_ADD, gs_event)
    call device_event_sync(gs_event)
    call rotate_cyc(a1_d, a2_d, a3_d, 0, coef)

    call blstx%apply(a(:, 1), n)
    call blsty%apply(a(:, 2), n)
    call blstz%apply(a(:, 3), n)

  end subroutine bicgstab_cpld_device_assemble

  !> Check a combined inner product for a BiCGStab breakdown.
  subroutine bicgstab_cpld_device_check_inner_product(inner_product, &
       norm_a, norm_b, quantity)
    real(kind=rp), intent(in) :: inner_product
    real(kind=rp), intent(in) :: norm_a
    real(kind=rp), intent(in) :: norm_b
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: large_norm, small_norm

    if (.not. ieee_is_finite(inner_product)) then
       call neko_error('Coupled BiCGStab failure: non-finite ' // &
            trim(quantity))
    end if

    large_norm = max(norm_a, norm_b)
    small_norm = min(norm_a, norm_b)
    if (large_norm .le. 0.0_rp .or. &
         abs(inner_product) / large_norm .le. NEKO_EPS * small_norm) then
       call neko_error('Coupled BiCGStab breakdown: near-zero ' // &
            trim(quantity))
    end if

  end subroutine bicgstab_cpld_device_check_inner_product

  !> Return the square root of a valid combined squared norm.
  function bicgstab_cpld_device_sqrt(value, quantity) result(root)
    real(kind=rp), intent(in) :: value
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: root

    if (.not. ieee_is_finite(value) .or. value .lt. 0.0_rp) then
       call neko_error('Coupled BiCGStab failure: invalid ' // &
            trim(quantity) // ' norm')
    end if
    root = sqrt(value)

  end function bicgstab_cpld_device_sqrt

end module bicgstab_cpld_device
