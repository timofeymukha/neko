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
  use bc_list, only : bc_list_t
  use coefs, only : coef_t
  use field, only : field_t
  use gs_ops, only : GS_OP_ADD
  use math, only : cmult, col2
  use mesh, only : mesh_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use pnpn2_mixed_ops, only : pnpn2_mixed_ops_t
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

  type(pnpn2_mixed_ops_t), pointer, save :: mixed_ops => null()
  type(bc_list_t), pointer, save :: bclst_x => null()
  type(bc_list_t), pointer, save :: bclst_y => null()
  type(bc_list_t), pointer, save :: bclst_z => null()

  public :: pnpn2_prs_ax_init, pnpn2_prs_ax_clear

contains

  !> Bind the pressure operator to the active hierarchy and velocity masks.
  subroutine pnpn2_prs_ax_init(mixed_space_ops, bclst_du, bclst_dv, bclst_dw)
    type(pnpn2_mixed_ops_t), target, intent(inout) :: mixed_space_ops
    type(bc_list_t), target, intent(inout) :: bclst_du
    type(bc_list_t), target, intent(inout) :: bclst_dv
    type(bc_list_t), target, intent(inout) :: bclst_dw

    mixed_ops => mixed_space_ops
    bclst_x => bclst_du
    bclst_y => bclst_dv
    bclst_z => bclst_dw
  end subroutine pnpn2_prs_ax_init

  !> Clear the module state used by the PnPn-2 pressure operator.
  subroutine pnpn2_prs_ax_clear()
    nullify(mixed_ops)
    nullify(bclst_x)
    nullify(bclst_y)
    nullify(bclst_z)
  end subroutine pnpn2_prs_ax_clear

  !> Apply \f$ D (\rho B)^{-1} G \f$ on \f$ Y_h \f$.
  subroutine pnpn2_prs_ax_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(field_t), pointer :: gx, gy, gz
    integer :: scratch_ids(3)
    integer :: n
    real(kind=rp) :: rho_inv

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2_prs_ax_t is currently CPU-only.')
    end if

    if (.not. associated(mixed_ops) .or. .not. associated(bclst_x) .or. &
         .not. associated(bclst_y) .or. .not. associated(bclst_z)) then
      call neko_error('PnPn-2 pressure operator used before initialization.')
    end if

    n = mixed_ops%coef_Xh%dof%size()
    rho_inv = coef%h1(1,1,1,1)

    call neko_scratch_registry%request_field(gx, scratch_ids(1), .false.)
    call neko_scratch_registry%request_field(gy, scratch_ids(2), .false.)
    call neko_scratch_registry%request_field(gz, scratch_ids(3), .false.)

    call mixed_ops%opgradt(gx%x, gy%x, gz%x, u)

    ! Nek's opbinv first restricts the transpose gradient to the homogeneous
    ! velocity-increment space, then assembles it.
    call bclst_x%apply_scalar(gx%x, n)
    call bclst_y%apply_scalar(gy%x, n)
    call bclst_z%apply_scalar(gz%x, n)

    call mixed_ops%coef_Xh%gs_h%op(gx, GS_OP_ADD)
    call mixed_ops%coef_Xh%gs_h%op(gy, GS_OP_ADD)
    call mixed_ops%coef_Xh%gs_h%op(gz, GS_OP_ADD)

    call col2(gx%x, mixed_ops%coef_Xh%Binv, n)
    call col2(gy%x, mixed_ops%coef_Xh%Binv, n)
    call col2(gz%x, mixed_ops%coef_Xh%Binv, n)
    call cmult(gx%x, rho_inv, n)
    call cmult(gy%x, rho_inv, n)
    call cmult(gz%x, rho_inv, n)

    call mixed_ops%opdiv(w, gx%x, gy%x, gz%x)
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
