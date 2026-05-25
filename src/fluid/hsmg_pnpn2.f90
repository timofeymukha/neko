! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   * Neither the name of the authors nor the names of its contributors may be
!     used to endorse or promote products derived from this software without
!     specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!> Pressure-space HSMG preconditioner for the Pn/Pn-2 pressure equation.
module hsmg_pnpn2
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_math, only : field_rzero
  use gather_scatter, only : gs_t, GS_OP_ADD
  use interpolation, only : interpolator_t
  use json_module, only : json_file
  use krylov, only : KSP_MAX_ITER, ksp_monitor_t, ksp_t, krylov_solver_factory
  use ax_product, only : ax_t, ax_helm_factory
  use jacobi, only : jacobi_t
  use math, only : add2, col2, copy
  use num_types, only : rp
  use precon, only : pc_t
  use schwarz, only : schwarz_t
  use space, only : GLL, space_t
  use utils, only : neko_error
  use zero_dirichlet, only : zero_dirichlet_t
  implicit none
  private

  !> Nek-style pressure-space HSMG preconditioner for Pn/Pn-2.
  !!
  !! This type owns the Pn/Pn-2 pressure-specific path instead of delegating to
  !! Neko's generic continuous-space HSMG. The pressure residual remains a local
  !! GL field; only the embedded mesh-1 Schwarz/FDM work uses the continuous Xh
  !! gather-scatter, matching the top local-solve part of Nek's `hsmg_solve`.
  type, public, extends(pc_t) :: hsmg_pnpn2_t
     type(coef_t), pointer :: c_Xh => null()
     type(space_t), pointer :: Yh => null()
     type(dofmap_t), pointer :: dm_Yh => null()
     type(dofmap_t), pointer :: dm_Xh => null()
     type(gs_t), pointer :: gs_Xh => null()
     type(bc_list_t) :: bclst_Xh
     type(schwarz_t) :: schwarz_Xh
     type(space_t) :: Xh_mg, Xh_crs
     type(dofmap_t) :: dm_mg, dm_crs
     type(gs_t) :: gs_mg, gs_crs
     type(coef_t) :: c_mg, c_crs
     type(zero_dirichlet_t) :: bc_mg, bc_crs
     type(bc_list_t) :: bclst_mg, bclst_crs
     type(schwarz_t) :: schwarz_mg
     type(interpolator_t) :: interp_top_mg, interp_mg_crs
     type(field_t) :: r_Yh, w_Yh
     type(field_t) :: r_Xh
     type(field_t) :: z_Xh
     type(field_t) :: r_mg, e_mg, w_mg
     type(field_t) :: r_crs, e_crs
     type(jacobi_t) :: pc_crs
     class(ksp_t), allocatable :: crs_solver
     class(ax_t), allocatable :: ax_crs
   contains
     procedure, pass(this) :: init => hsmg_pnpn2_init
     procedure, pass(this) :: free => hsmg_pnpn2_free
     procedure, pass(this) :: solve => hsmg_pnpn2_solve
     procedure, pass(this) :: update => hsmg_pnpn2_update
  end type hsmg_pnpn2_t

contains

  !> Initialize the pressure-space HSMG cycle on the Pn/Pn-2 pressure grid.
  subroutine hsmg_pnpn2_init(this, c_Xh, Yh, dm_Yh, dm_Xh, gs_Xh, bcs_prs_Xh, &
       params)
    class(hsmg_pnpn2_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c_Xh
    type(space_t), target, intent(inout) :: Yh
    type(dofmap_t), target, intent(inout) :: dm_Yh
    type(dofmap_t), target, intent(inout) :: dm_Xh
    type(gs_t), target, intent(inout) :: gs_Xh
    type(bc_list_t), target, intent(inout) :: bcs_prs_Xh
    type(json_file), intent(inout) :: params
    class(bc_t), pointer :: bc_i
    integer :: i, lx_crs, lx_mg

    call this%free()

    this%c_Xh => c_Xh
    this%Yh => Yh
    this%dm_Yh => dm_Yh
    this%dm_Xh => dm_Xh
    this%gs_Xh => gs_Xh

    lx_crs = 2
    lx_mg = hsmg_pnpn2_mid_lx(c_Xh%Xh%lx)

    call this%bclst_Xh%init()
    do i = 1, bcs_prs_Xh%size()
      if (bcs_prs_Xh%strong(i)) then
        bc_i => bcs_prs_Xh%get(i)
        call this%bclst_Xh%append(bc_i)
      end if
    end do

    call this%r_Xh%init(dm_Xh, 'hsmg_pnpn2_r_Xh')
    call this%z_Xh%init(dm_Xh, 'hsmg_pnpn2_z_Xh')
    call this%schwarz_Xh%init(c_Xh%Xh, dm_Xh, gs_Xh, this%bclst_Xh, c_Xh%msh)

    call this%r_Yh%init(dm_Yh, 'hsmg_pnpn2_r_Yh')
    call this%w_Yh%init(dm_Yh, 'hsmg_pnpn2_w_Yh')

    call this%Xh_mg%init(GLL, lx_mg, lx_mg, lx_mg)
    call this%dm_mg%init(c_Xh%msh, this%Xh_mg)
    call this%gs_mg%init(this%dm_mg)
    call this%c_mg%init(this%gs_mg)
    call this%r_mg%init(this%dm_mg, 'hsmg_pnpn2_r_mg')
    call this%e_mg%init(this%dm_mg, 'hsmg_pnpn2_e_mg')
    call this%w_mg%init(this%dm_mg, 'hsmg_pnpn2_w_mg')

    call this%Xh_crs%init(GLL, lx_crs, lx_crs, lx_crs)
    call this%dm_crs%init(c_Xh%msh, this%Xh_crs)
    call this%gs_crs%init(this%dm_crs)
    call this%c_crs%init(this%gs_crs)
    call this%r_crs%init(this%dm_crs, 'hsmg_pnpn2_r_crs')
    call this%e_crs%init(this%dm_crs, 'hsmg_pnpn2_e_crs')

    call this%bc_mg%init_base(this%c_mg)
    call this%bc_crs%init_base(this%c_crs)
    do i = 1, bcs_prs_Xh%size()
      if (bcs_prs_Xh%strong(i)) then
        bc_i => bcs_prs_Xh%get(i)
        call this%bc_mg%mark_facets(bc_i%marked_facet)
        call this%bc_crs%mark_facets(bc_i%marked_facet)
      end if
    end do
    call this%bc_mg%finalize()
    call this%bc_crs%finalize()
    call this%bclst_mg%init()
    call this%bclst_crs%init()
    call this%bclst_mg%append(this%bc_mg)
    call this%bclst_crs%append(this%bc_crs)

    call this%schwarz_mg%init(this%Xh_mg, this%dm_mg, this%gs_mg, &
         this%bclst_mg, c_Xh%msh)
    call this%interp_top_mg%init(Yh, this%Xh_mg)
    call this%interp_mg_crs%init(this%Xh_mg, this%Xh_crs)
    call ax_helm_factory(this%ax_crs, full_formulation = .false.)
    call this%pc_crs%init(this%c_crs, this%dm_crs, this%gs_crs)
    call krylov_solver_factory(this%crs_solver, this%dm_crs%size(), 'cg', &
         KSP_MAX_ITER, M = this%pc_crs, monitor = .false.)
  end subroutine hsmg_pnpn2_init

  !> Release pressure-space HSMG storage.
  subroutine hsmg_pnpn2_free(this)
    class(hsmg_pnpn2_t), intent(inout) :: this

    call this%schwarz_Xh%free()
    call this%schwarz_mg%free()
    call this%interp_top_mg%free()
    call this%interp_mg_crs%free()
    call this%pc_crs%free()
    if (allocated(this%crs_solver)) then
      call this%crs_solver%free()
      deallocate(this%crs_solver)
    end if
    if (allocated(this%ax_crs)) then
      deallocate(this%ax_crs)
    end if
    call this%r_Yh%free()
    call this%w_Yh%free()
    call this%r_Xh%free()
    call this%z_Xh%free()
    call this%r_mg%free()
    call this%e_mg%free()
    call this%w_mg%free()
    call this%r_crs%free()
    call this%e_crs%free()
    call this%c_mg%free()
    call this%c_crs%free()
    call this%gs_mg%free()
    call this%gs_crs%free()
    call this%dm_mg%free()
    call this%dm_crs%free()
    call this%Xh_mg%free()
    call this%Xh_crs%free()
    call this%bc_mg%free()
    call this%bc_crs%free()
    call this%bclst_mg%free()
    call this%bclst_crs%free()
    call this%bclst_Xh%free()
    nullify(this%c_Xh)
    nullify(this%Yh)
    nullify(this%dm_Yh)
    nullify(this%dm_Xh)
    nullify(this%gs_Xh)
  end subroutine hsmg_pnpn2_free

  !> Apply the top Nek pressure-space Schwarz/FDM solve.
  !!
  !! This corresponds to Nek's `local_solves_fdm`: embed the GL pressure
  !! residual into the interior of the mesh-1 work field, perform the mesh-1
  !! overlapping Schwarz/FDM solve, and extract the interior back to GL.
  subroutine hsmg_pnpn2_solve(this, z, r, n)
    class(hsmg_pnpn2_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: z(n)
    real(kind=rp), intent(inout) :: r(n)
    integer :: e, i, j, k
    integer :: ix, iy, iz, nelv
    integer :: nx, ny, nz, nxyz_y
    integer :: id_x, id_y
    type(ksp_monitor_t) :: crs_info

    if (.not. associated(this%c_Xh)) then
      call neko_error('hsmg_pnpn2 used before initialization.')
    end if

    nx = this%c_Xh%Xh%lx
    ny = nx - 2
    nz = nx - 2
    nelv = this%c_Xh%msh%nelv
    if (this%c_Xh%msh%gdim .ne. 3) then
      nz = 1
    end if
    nxyz_y = ny * ny * nz

    if (n .ne. nxyz_y * nelv) then
      call neko_error('hsmg_pnpn2 received an unexpected pressure vector size.')
    end if

    call field_rzero(this%r_Xh)
    call field_rzero(this%z_Xh)

    do e = 1, nelv
      do k = 1, nz
        iz = k
        if (this%c_Xh%msh%gdim .eq. 3) iz = k + 1
        do j = 1, ny
          do i = 1, ny
            ix = i + 1
            iy = j + 1
            id_y = i + (j - 1) * ny + (k - 1) * ny * ny + (e - 1) * nxyz_y
            id_x = ix + (iy - 1) * nx + (iz - 1) * nx * nx + &
                 (e - 1) * nx * nx * this%c_Xh%Xh%lz
            this%r_Xh%x(id_x,1,1,1) = r(id_y)
          end do
        end do
      end do
    end do

    call this%schwarz_Xh%compute(this%z_Xh%x, this%r_Xh%x)
    call copy(this%r_Yh%x, r, n)

    do e = 1, nelv
      do k = 1, nz
        iz = k
        if (this%c_Xh%msh%gdim .eq. 3) iz = k + 1
        do j = 1, ny
          do i = 1, ny
            ix = i + 1
            iy = j + 1
            id_y = i + (j - 1) * ny + (k - 1) * ny * ny + (e - 1) * nxyz_y
            id_x = ix + (iy - 1) * nx + (iz - 1) * nx * nx + &
                 (e - 1) * nx * nx * this%c_Xh%Xh%lz
            z(id_y) = this%z_Xh%x(id_x,1,1,1)
          end do
        end do
      end do
    end do

    call this%interp_top_mg%map(this%r_mg%x, this%r_Yh%x, nelv, this%Xh_mg)
    call this%gs_mg%op(this%r_mg%x, this%dm_mg%size(), GS_OP_ADD)

    call this%schwarz_mg%compute(this%e_mg%x, this%r_mg%x)
    call col2(this%r_mg%x, this%c_mg%mult, this%dm_mg%size())

    call this%interp_mg_crs%map(this%r_crs%x, this%r_mg%x, nelv, this%Xh_crs)
    call this%bclst_crs%apply_scalar(this%r_crs%x, this%dm_crs%size())

    crs_info = this%crs_solver%solve(this%ax_crs, this%e_crs, this%r_crs%x, &
         this%dm_crs%size(), this%c_crs, this%bclst_crs, this%gs_crs, 10)
    call this%bclst_crs%apply_scalar(this%e_crs%x, this%dm_crs%size())

    call this%interp_mg_crs%map(this%w_mg%x, this%e_crs%x, nelv, this%Xh_mg)
    call add2(this%e_mg%x, this%w_mg%x, this%dm_mg%size())

    call this%interp_top_mg%map(this%w_Yh%x, this%e_mg%x, nelv, this%Yh)
    call add2(z, this%w_Yh%x, n)
  end subroutine hsmg_pnpn2_solve

  !> Refresh HSMG coefficients.
  subroutine hsmg_pnpn2_update(this)
    class(hsmg_pnpn2_t), intent(inout) :: this

    ! Coefficients for the top Schwarz/FDM solve are initialized from Xh.
  end subroutine hsmg_pnpn2_update

  !> Return Nek's middle HSMG point count for the current velocity GLL count.
  integer function hsmg_pnpn2_mid_lx(lx1) result(lx_mg)
    integer, intent(in) :: lx1
    integer, parameter :: mgn2(10) = [1, 2, 2, 2, 2, 3, 3, 5, 5, 5]
    integer :: mg_nx

    if (lx1 .le. 10) then
      mg_nx = mgn2(max(1, lx1))
    else
      mg_nx = 2 * ((lx1 - 2) / 4) + 1
    end if
    if (lx1 .eq. 8) then
      mg_nx = 3
    end if

    lx_mg = mg_nx + 1
  end function hsmg_pnpn2_mid_lx

end module hsmg_pnpn2
