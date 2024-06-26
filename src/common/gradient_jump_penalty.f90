! Copyright (c) 2024, The Neko Authors
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
!
!> Implements `gradient_jump_penalty_t`.
module gradient_jump_penalty
  use num_types, only : rp
  use speclib, only : pnleg, pndleg
  use utils, only : neko_error
  use math
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use scratch_registry, only : neko_scratch_registry
  use element, only : element_t
  use hex
  use quad
  use operators, only : dudxyz
  use gs_ops, only : GS_OP_ADD

  implicit none
  private

  !> Implements the gradient jump penalty.
  !! @note Reference DOI: 10.1016/j.cma.2021.114200
  type, public :: gradient_jump_penalty_t
     !> Coefficient of the penalty.
     real :: tau
     !> Polynomial order
     integer :: p
     !> Penalty terms
     real(kind=rp), allocatable :: penalty(:, :, :, :)
     !> Work array to store penalty terms on each element
     real(kind=rp), allocatable :: penalty_el(:, :, :)
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Gradient jump on the elementary interface (zero inside each element)
     type(field_t), pointer :: G
     !> Flux of the quantity
     type(field_t), pointer :: flux1, flux2, flux3
     !> Expanded array of facet normal
     type(field_t), pointer :: n1, n2, n3
     !> Number of facet in elements and its maximum
     integer, allocatable :: n_facet(:)
     integer :: n_facet_max
     !> Length scale for element regarding a facet
     real(kind=rp), allocatable :: h(:)
     !> Polynomial Quadrature
     real(kind=rp), allocatable :: w(:), w2(:,:)
     !> Polynomial evaluated at collocation points
     real(kind=rp), allocatable :: phi(:,:)
     !> The first derivative of polynomial at two ends of the interval
     real(kind=rp), allocatable :: dphidxi(:,:)

  contains
     !> Constructor.
     procedure, pass(this) :: init => gradient_jump_penalty_init
     !> Destructor
     procedure, pass(this) :: free => gradient_jump_penalty_free
     !> Compute gradient jump penalty term.
     procedure, pass(this) :: compute => gradient_jump_penalty_compute

  end type gradient_jump_penalty_t

contains
  !> Constructor.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  subroutine gradient_jump_penalty_init(this, dofmap, coef)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    
    class(element_t), pointer :: ep
    integer :: temp_indices(4), i, j
    real(kind=rp), allocatable :: zg(:) ! Quadrature points
    
    call this%free()

    this%p = dofmap%xh%lx - 1
    if (this%p .gt. 1) then
       this%tau = 0.8_rp * (this%p + 1) ** (-4)
    else
       this%tau = 0.02_rp
    end if
    this%coef => coef

    allocate(this%penalty(this%p + 1, this%p + 1 , this%p + 1 , this%coef%msh%nelv))
    allocate(this%penalty_el(this%p + 1, this%p + 1 , this%p + 1))
    
    allocate(this%n_facet(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          this%n_facet(i) = 6
       type is (quad_t)
          call neko_error("Only Hexahedral element is supported now for gradient jump penalty")
       end select
    end do
    this%n_facet_max = maxval(this%n_facet)

    allocate(this%h(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call eval_h_hex(this%h(i), ep)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size evaluation is not supported for quad_t")
       end select
    end do
    
    allocate(zg(dofmap%xh%lx))
    allocate(this%w(dofmap%xh%lx))
    allocate(this%w2(dofmap%xh%lx, dofmap%xh%lx))
    allocate(this%phi(this%p + 1, this%p + 1))
    allocate(this%dphidxi(this%p + 1, this%p + 1))

    zg = dofmap%xh%zg(:,1)
    this%w = dofmap%xh%wx
    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%w2(i, j) = this%w(i) * this%w(j)
       end do
    end do

    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%phi(j,i) = pnleg(zg(j),i-1)
          this%dphidxi(j,i) = pndleg(zg(j),i-1)
       end do
    end do

    call neko_scratch_registry%request_field(this%G, temp_indices(1))
    call neko_scratch_registry%request_field(this%flux1, temp_indices(2))
    call neko_scratch_registry%request_field(this%flux2, temp_indices(3))
    call neko_scratch_registry%request_field(this%flux3, temp_indices(4))
    call neko_scratch_registry%request_field(this%n1, temp_indices(2))
    call neko_scratch_registry%request_field(this%n2, temp_indices(3))
    call neko_scratch_registry%request_field(this%n3, temp_indices(4))
    
    ! formulate n1, n2 and n3
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call expand_normal_hex_el(this%n1%x(:, :, :, i), this%coef%nx(:, :, :, i), this%p + 1)
          call expand_normal_hex_el(this%n2%x(:, :, :, i), this%coef%ny(:, :, :, i), this%p + 1)
          call expand_normal_hex_el(this%n3%x(:, :, :, i), this%coef%nz(:, :, :, i), this%p + 1)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size evaluation is not supported for quad_t")
       end select
    end do

  end subroutine gradient_jump_penalty_init
  
  !> Evaluate h for each element
  subroutine eval_h_hex(h_el, ep)
    real(kind=rp), intent(inout) :: h_el
    type(hex_t), pointer, intent(in) :: ep
    
    !! todo: estimation of the length scale of the mesh could be more elegant
    h_el = ep%diameter()

  end subroutine eval_h_hex

  !> Expanding facet normal array into a field of (lx, lx, lx) for each element
  subroutine expand_normal_hex_el(n1_el, coefnx_el, lx)
    integer :: lx
    real(kind=rp), intent(inout) :: n1_el(lx, lx, lx)
    real(kind=rp), intent(in) :: coefnx_el(lx, lx, 6)

    n1_el(1,:,:) = coefnx_el(:,:,1)
    n1_el(lx,:,:) = coefnx_el(:,:,2)
    n1_el(:,1,:) = coefnx_el(:,:,3)
    n1_el(:,lx,:) = coefnx_el(:,:,4)
    n1_el(:,:,1) = coefnx_el(:,:,5)
    n1_el(:,:,lx) = coefnx_el(:,:,6)

  end subroutine expand_normal_hex_el

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    
    if (allocated(this%penalty)) then
       deallocate(this%penalty)
    end if
    if (allocated(this%penalty_el)) then
       deallocate(this%penalty_el)
    end if
    if (allocated(this%h)) then
       deallocate(this%h)
    end if
    if (allocated(this%n_facet)) then
       deallocate(this%n_facet)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%w2)) then
       deallocate(this%w2)
    end if
    if (allocated(this%phi)) then
       deallocate(this%phi)
    end if
    if (allocated(this%dphidxi)) then
       deallocate(this%dphidxi)
    end if
    nullify(this%coef)
    nullify(this%G)
    nullify(this%flux1)
    nullify(this%flux2)
    nullify(this%flux3)
    nullify(this%n1)
    nullify(this%n2)
    nullify(this%n3)
    
  end subroutine gradient_jump_penalty_free

  !> Compute the gradient jump penalty term.
  !! @param u x-velocity
  !! @param v y-velocity
  !! @param w z-velocity
  !! @param s quantity of interest
  subroutine gradient_jump_penalty_compute(this, u, v, w, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w, s

    class(element_t), pointer :: ep
    integer :: i

    call G_compute(this, s)
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call gradient_jump_penalty_compute_hex_el(this, u, v, w, s)
          this%penalty(:, :, :, i) = this%penalty_el
       type is (quad_t)
          call neko_error("Only Hexahedral element is supported now for gradient jump penalty")
       end select
    end do

  end subroutine gradient_jump_penalty_compute

  !> Compute the gradient jump penalty term for a single hexatedral element.
  !! @param u x-velocity
  !! @param v y-velocity
  !! @param w z-velocity
  !! @param s The quantity of interest
  subroutine gradient_jump_penalty_compute_hex_el(this, u, v, w, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w, s

    integer :: i

    write(*,*) "call gradient_jump_penalty_compute_hex_el"

  end subroutine gradient_jump_penalty_compute_hex_el

  !> Compute the average of the flux over facets
  !! @param s The quantity of interest
  subroutine G_compute(this, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: s

    integer :: i

    call dudxyz(this%flux1%x, s%x, this%coef%drdx, this%coef%dsdx, this%coef%dtdx, this%coef)
    call dudxyz(this%flux2%x, s%x, this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
    call dudxyz(this%flux3%x, s%x, this%coef%drdz, this%coef%dsdz, this%coef%dtdz, this%coef)

    call col2(this%flux1%x, this%n1%x, this%coef%dof%size())
    call col2(this%flux2%x, this%n2%x, this%coef%dof%size())
    call col2(this%flux3%x, this%n3%x, this%coef%dof%size())

    call this%coef%gs_h%op(this%flux1%x, this%coef%dof%size(), GS_OP_ADD)
    call this%coef%gs_h%op(this%flux2%x, this%coef%dof%size(), GS_OP_ADD)
    call this%coef%gs_h%op(this%flux3%x, this%coef%dof%size(), GS_OP_ADD)

    call add3(this%G%x, this%flux1%x, this%flux2%x, this%coef%dof%size())
    call add2(this%G%x, this%flux3%x, this%coef%dof%size())

  end subroutine G_compute

  ! !> Integrate over a facet
  ! !! @param f Integrant
  ! subroutine integrate_over_facet(f)

  ! end subroutine integrate_over_facet

end module gradient_jump_penalty