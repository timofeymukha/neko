! Copyright (c) 2025, The Neko Authors
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
!> Implements the PnPn-2 pressure operator on the lower-order pressure space.
module pnpn2_prs_ax
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use field, only : field_t
  use interpolation, only : interpolator_t
  use math, only : add2, col2
  use mesh, only : mesh_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use operators, only : cdtp, opgrad
  use scratch_registry, only : neko_scratch_registry
  use space, only : space_t
  use utils, only : neko_error
  implicit none
  private

  type, public, extends(ax_t) :: pnpn2_prs_ax_t
   contains
     procedure, nopass :: compute => pnpn2_prs_ax_compute
     procedure, pass(this) :: compute_vector => pnpn2_prs_ax_compute_vector
  end type pnpn2_prs_ax_t

  type(space_t), pointer, save :: velocity_space => null()
  type(space_t), pointer, save :: pressure_space => null()
  type(coef_t), pointer, save :: velocity_coef => null()
  type(interpolator_t), pointer, save :: pressure_interpolator => null()

  public :: pnpn2_prs_ax_init, pnpn2_prs_ax_clear

contains

  !> Bind the PnPn-2 pressure operator to the active unequal-order hierarchy.
  subroutine pnpn2_prs_ax_init(Xh, Yh, c_Xh, prs_interp)
    type(space_t), target, intent(inout) :: Xh
    type(space_t), target, intent(inout) :: Yh
    type(coef_t), target, intent(inout) :: c_Xh
    type(interpolator_t), target, intent(inout) :: prs_interp

    velocity_space => Xh
    pressure_space => Yh
    velocity_coef => c_Xh
    pressure_interpolator => prs_interp
  end subroutine pnpn2_prs_ax_init

  !> Clear the module state used by the PnPn-2 pressure operator.
  subroutine pnpn2_prs_ax_clear()
    nullify(velocity_space)
    nullify(pressure_space)
    nullify(velocity_coef)
    nullify(pressure_interpolator)
  end subroutine pnpn2_prs_ax_clear

  !> Apply the PnPn-2 projection operator \f$ D B^{-1} G \f$ on \f$ Y_h \f$.
  subroutine pnpn2_prs_ax_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(field_t), pointer :: p_Xh, gx, gy, gz, div_Xh
    integer :: scratch_ids(5)
    integer :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2_prs_ax_t is currently CPU-only.')
    end if

    if ((.not. associated(velocity_space)) .or. &
         (.not. associated(pressure_space)) .or. &
         (.not. associated(velocity_coef)) .or. &
         (.not. associated(pressure_interpolator))) then
      call neko_error('PnPn-2 pressure operator used before initialization.')
    end if

    n = velocity_coef%dof%size()

    call neko_scratch_registry%request_field(p_Xh, scratch_ids(1), .false.)
    call neko_scratch_registry%request_field(gx, scratch_ids(2), .false.)
    call neko_scratch_registry%request_field(gy, scratch_ids(3), .false.)
    call neko_scratch_registry%request_field(gz, scratch_ids(4), .false.)
    call neko_scratch_registry%request_field(div_Xh, scratch_ids(5), .false.)

    call pressure_interpolator%map(p_Xh%x, u, msh%nelv, velocity_space)

    call opgrad(gx%x, gy%x, gz%x, p_Xh%x, velocity_coef)

    call col2(gx%x, velocity_coef%Binv, n)
    call col2(gy%x, velocity_coef%Binv, n)
    call col2(gz%x, velocity_coef%Binv, n)

    call cdtp(div_Xh%x, gx%x, velocity_coef%drdx, velocity_coef%dsdx, &
         velocity_coef%dtdx, velocity_coef)
    call cdtp(gy%x, gy%x, velocity_coef%drdy, velocity_coef%dsdy, &
         velocity_coef%dtdy, velocity_coef)
    call cdtp(gz%x, gz%x, velocity_coef%drdz, velocity_coef%dsdz, &
         velocity_coef%dtdz, velocity_coef)

    call add2(div_Xh%x, gy%x, n)
    call add2(div_Xh%x, gz%x, n)

    call pressure_interpolator%map(w, div_Xh%x, msh%nelv, pressure_space)

    call neko_scratch_registry%relinquish_field(scratch_ids)
  end subroutine pnpn2_prs_ax_compute

  !> Apply the scalar PnPn-2 pressure operator componentwise.
  subroutine pnpn2_prs_ax_compute_vector(this, au, av, aw, u, v, w, coef, msh, &
       Xh)
    class(pnpn2_prs_ax_t), intent(in) :: this
    type(space_t), intent(in) :: Xh
    type(mesh_t), intent(in) :: msh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call pnpn2_prs_ax_compute(au, u, coef, msh, Xh)
    call pnpn2_prs_ax_compute(av, v, coef, msh, Xh)
    call pnpn2_prs_ax_compute(aw, w, coef, msh, Xh)
  end subroutine pnpn2_prs_ax_compute_vector

end module pnpn2_prs_ax
