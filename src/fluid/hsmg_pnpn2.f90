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
!!
!! This module contains the milestone pressure preconditioner used by the
!! unequal-order Pn/Pn-2 formulation when the user requests the `hsmg`
!! pressure preconditioner.
!!
!! The implementation deliberately mirrors the structure of Nek5000's H1
!! multigrid pressure preconditioner as far as the current Pn/Pn-2 split
!! allows:
!! - the incoming pressure residual lives on the discontinuous local `Y_h`
!!   pressure grid;
!! - the top local Schwarz/FDM stage is applied on an embedded `X_h`
!!   mesh-1 field;
!! - the remaining V-cycle works on continuous H1-like spaces on successively
!!   coarser GLL grids;
!! - the coarse solve is performed by [pnpn2_coarse_direct_t]
!!   (#pnpn2_coarse_direct::pnpn2_coarse_direct_t), which now assembles a
!!   dedicated sparse coarse operator on the unique coarse dofs and solves it
!!   with a dedicated coarse iteration.
!!
!! A reader should think of this module as the orchestration layer.  It does
!! not itself build element matrices or perform dense linear algebra.  Its job
!! is to:
!! 1. own the temporary fields for each level,
!! 2. construct the intermediate spaces and their boundary masks,
!! 3. apply the top Schwarz solve,
!! 4. perform the restriction/prolongation steps of the V-cycle, and
!! 5. hand the true coarse solve to the dedicated coarse solver object.
module hsmg_pnpn2
  use ax_product, only : ax_t, ax_helm_factory
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use coefs, only : coef_t
  use device_identity, only : device_ident_t
  use device_jacobi, only : device_jacobi_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_math, only : field_rzero
  use gather_scatter, only : gs_t, GS_OP_ADD
  use identity, only : ident_t
  use interpolation, only : interpolator_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use jacobi, only : jacobi_t
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, krylov_solver_factory
  use math, only : add2, col2, copy
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  use num_types, only : rp
  use fdm, only : fdm_t
  use pnpn2_coarse_direct, only : pnpn2_coarse_direct_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use profiler, only : profiler_start_region, profiler_end_region
  use schwarz, only : schwarz_t
  use space, only : GLL, space_t
  use sx_jacobi, only : sx_jacobi_t
  use tree_amg_multigrid, only : tamg_solver_t
  use utils, only : neko_error
  use zero_dirichlet, only : zero_dirichlet_t
  implicit none
  private

  !> Nek-style pressure-space HSMG preconditioner for Pn/Pn-2.
  !!
  !! The type groups together all state needed by one application of the
  !! pressure-space V-cycle:
  !! - pointers to the fine-grid data owned by the surrounding fluid scheme,
  !! - constructed multigrid levels (`Xh_mg` and `Xh_crs`),
  !! - level-specific coefficients and boundary masks,
  !! - work buffers for residuals and corrections, and
  !! - the specialised explicit coarse solver.
  !!
  !! The implementation does not expose the individual substeps to the outside.
  !! The public entry points are therefore only the standard preconditioner
  !! operations: `init`, `free`, `solve`, and `update`.
  type, public, extends(pc_t) :: hsmg_pnpn2_t
    !> Fine-grid coefficients on the embedded mesh-1 `X_h` space.
    type(coef_t), pointer :: c_Xh => null()
    !> Fine-grid pressure space `Y_h`, i.e. the local GL pressure grid.
    type(space_t), pointer :: Yh => null()
    !> Degree-of-freedom map for the local pressure grid.
    type(dofmap_t), pointer :: dm_Yh => null()
    !> Degree-of-freedom map for the embedded mesh-1 grid.
    type(dofmap_t), pointer :: dm_Xh => null()
    !> Gather-scatter handle for the embedded mesh-1 grid.
    type(gs_t), pointer :: gs_Xh => null()
    !> Strong pressure boundary conditions transferred to `X_h`.
    type(bc_list_t) :: bclst_Xh
    !> Top-level embedded mesh-1 FDM solver matching Nek `local_solves_fdm`.
    type(fdm_t) :: top_fdm_Xh
    !> Middle and coarse GLL spaces of the H1 hierarchy.
    type(space_t) :: Xh_mg, Xh_crs
    !> Degree-of-freedom maps associated with the middle and coarse spaces.
    type(dofmap_t) :: dm_mg, dm_crs
    !> Gather-scatter handles associated with the middle and coarse spaces.
    type(gs_t) :: gs_mg, gs_crs
    !> H1-form coefficients on the middle and coarse levels.
    type(coef_t) :: c_mg, c_crs
    !> Zero Dirichlet masks transferred to the middle and coarse levels.
    type(zero_dirichlet_t) :: bc_mg, bc_crs
    !> Boundary-condition lists used by the lower levels.
    type(bc_list_t) :: bclst_mg, bclst_crs
    !> Schwarz/FDM smoother on the middle level.
    type(schwarz_t) :: schwarz_mg
    !> Interpolation operators for the downward and upward legs.
    type(interpolator_t) :: interp_top_mg, interp_mg_crs
    !> Fine-grid work vectors on the local pressure grid.
    type(field_t) :: r_Yh, w_Yh
    !> Fine-grid work vectors on the embedded `X_h` grid.
    type(field_t) :: r_Xh
    type(field_t) :: z_Xh
    !> Nek `do_weight_op`-style multiplicity weights on the embedded pressure shell.
    type(field_t) :: top_wt_Xh
    !> Middle-level residual, correction, and prolongation work fields.
    type(field_t) :: r_mg, e_mg, w_mg
    !> Coarse-level residual and correction fields.
    type(field_t) :: r_crs, e_crs
    !> Dedicated sparse coarse solver for the Pn/Pn-2 H1 coarse grid.
    type(pnpn2_coarse_direct_t) :: crs_solver
    !> Generic coarse-grid Krylov solver reused from the standard HSMG stack.
    class(ksp_t), allocatable :: crs_ksp
    !> Generic coarse-grid preconditioner reused from the standard HSMG stack.
    class(pc_t), allocatable :: pc_crs_solver
    !> Optional TreeAMG coarse solver reused from the standard HSMG stack.
    type(tamg_solver_t), allocatable :: amg_solver
    !> Coarse-grid H1 operator used by the generic coarse solvers.
    class(ax_t), allocatable :: ax
    !> Number of iterations requested for the reusable coarse solver path.
    integer :: crs_niter = 10
    !> True once the selected coarse solver has been constructed successfully.
    logical :: crs_solver_initialized = .false.
    !> True when the dedicated Pn/Pn-2 coarse solver path is active.
    logical :: use_direct_coarse = .true.
   contains
     !> Initialise all levels of the pressure-space V-cycle.
     procedure, pass(this) :: init => hsmg_pnpn2_init
     !> Release all storage owned by the preconditioner.
     procedure, pass(this) :: free => hsmg_pnpn2_free
     !> Apply one pressure-space V-cycle.
     procedure, pass(this) :: solve => hsmg_pnpn2_solve
     !> Refresh level coefficients after the surrounding fluid update.
     procedure, pass(this) :: update => hsmg_pnpn2_update
     !> Force a coefficient object into the H1 form used by this preconditioner.
     procedure, pass(this), private :: set_h1_coeffs => &
          hsmg_pnpn2_set_h1_coeffs
     !> Build the top-level pressure-shell weights used after the local FDM solve.
     procedure, pass(this), private :: build_top_weight => &
          hsmg_pnpn2_build_top_weight
     !> Apply the top-level pressure-shell weights in-place.
     procedure, pass(this), private :: weight_top_pressure => &
          hsmg_pnpn2_weight_top_pressure
     !> Apply Nek's coarse masked-DSS weight field.
     procedure, pass(this), private :: weight_coarse_mask => &
          hsmg_pnpn2_weight_coarse_mask
     !> Apply Nek-style restriction weights only on level boundaries.
     procedure, pass(this), private :: weight_restriction_boundary => &
          hsmg_pnpn2_weight_restriction_boundary
     final :: hsmg_pnpn2_finalize
  end type hsmg_pnpn2_t

contains

  !> Initialise the pressure-space HSMG cycle on the Pn/Pn-2 pressure grid.
  !!
  !! @param c_Xh Fine-grid coefficient object on the embedded mesh-1 space.
  !! @param Yh Fine-grid local pressure space.
  !! @param dm_Yh Degree-of-freedom map of the local pressure space.
  !! @param dm_Xh Degree-of-freedom map of the embedded mesh-1 space.
  !! @param gs_Xh Gather-scatter handle on the embedded mesh-1 space.
  !! @param bcs_prs_Xh Strong pressure boundary conditions inherited from the
  !! surrounding fluid scheme.
  !! @param params Preconditioner parameter subtree from the case file.
  !!
  !! The routine builds the two lower HSMG levels explicitly.  It also transfers
  !! any strong pressure boundary conditions from the fine grid to the lower
  !! levels so that the restriction/prolongation steps and the coarse solve all
  !! see a consistent H1 boundary mask.
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
    logical :: coarse_null_space
    character(len=:), allocatable :: crs_solver_type, crs_pc_type
    logical :: crs_monitor
    integer :: crs_tamg_lvls, crs_tamg_itrs, crs_tamg_cheby_degree

    call this%free()

    this%c_Xh => c_Xh
    this%Yh => Yh
    this%dm_Yh => dm_Yh
    this%dm_Xh => dm_Xh
    this%gs_Xh => gs_Xh

    lx_crs = 2
    lx_mg = hsmg_pnpn2_mid_lx(c_Xh%Xh%lx)
    coarse_null_space = .true.
    call json_get_or_default(params, 'coarse_grid.solver', crs_solver_type, 'direct')
    call json_get_or_default(params, 'coarse_grid.preconditioner', crs_pc_type, &
         'jacobi')
    call json_get_or_default(params, 'coarse_grid.monitor', crs_monitor, .false.)
    call json_get_or_default(params, 'coarse_grid.iterations', this%crs_niter, 10)
    call json_get_or_default(params, 'coarse_grid.levels', crs_tamg_lvls, 3)
    call json_get_or_default(params, 'coarse_grid.cheby_degree', &
         crs_tamg_cheby_degree, 4)
    call json_get_or_default(params, 'coarse_grid.iterations', crs_tamg_itrs, 1)
    this%use_direct_coarse = (trim(crs_solver_type) .eq. 'direct')

    ! Collect the strong fine-grid pressure conditions that define the `X_h`
    ! solve.  The same logical information is later transferred to the middle
    ! and coarse levels.
    call this%bclst_Xh%init()
    do i = 1, bcs_prs_Xh%size()
      if (bcs_prs_Xh%strong(i)) then
        bc_i => bcs_prs_Xh%get(i)
        call this%bclst_Xh%append(bc_i)
        coarse_null_space = .false.
      end if
    end do

    ! Top level:
    ! - `r_Xh` and `z_Xh` are the embedded mesh-1 residual/correction fields;
    ! - `top_fdm_Xh` reproduces Nek's `local_solves_fdm` tensor solve on `X_h`.
    call this%r_Xh%init(dm_Xh, 'hsmg_pnpn2_r_Xh')
    call this%z_Xh%init(dm_Xh, 'hsmg_pnpn2_z_Xh')
    call this%top_wt_Xh%init(dm_Xh, 'hsmg_pnpn2_top_wt_Xh')
    call this%top_fdm_Xh%init_sem(c_Xh%Xh, dm_Xh, gs_Xh)
    call this%build_top_weight()

    ! Fine pressure-grid work buffers.  These remain local to each element and
    ! therefore use `dm_Yh`, not the continuous `X_h` dofmap.
    call this%r_Yh%init(dm_Yh, 'hsmg_pnpn2_r_Yh')
    call this%w_Yh%init(dm_Yh, 'hsmg_pnpn2_w_Yh')

    ! Middle H1 level.  This is the first true multigrid level below the top
    ! embedded `X_h` Schwarz solve.
    call this%Xh_mg%init(GLL, lx_mg, lx_mg, lx_mg)
    call this%dm_mg%init(c_Xh%msh, this%Xh_mg)
    call this%gs_mg%init(this%dm_mg)
    call this%c_mg%init(this%gs_mg)
    call this%set_h1_coeffs(this%c_mg)
    call this%r_mg%init(this%dm_mg, 'hsmg_pnpn2_r_mg')
    call this%e_mg%init(this%dm_mg, 'hsmg_pnpn2_e_mg')
    call this%w_mg%init(this%dm_mg, 'hsmg_pnpn2_w_mg')

    ! Coarse H1 level.  This level is not solved iteratively; instead it is
    ! handed to the explicit condensed coarse solver.
    call this%Xh_crs%init(GLL, lx_crs, lx_crs, lx_crs)
    call this%dm_crs%init(c_Xh%msh, this%Xh_crs)
    call this%gs_crs%init(this%dm_crs)
    call this%c_crs%init(this%gs_crs)
    call this%set_h1_coeffs(this%c_crs)
    call this%r_crs%init(this%dm_crs, 'hsmg_pnpn2_r_crs')
    call this%e_crs%init(this%dm_crs, 'hsmg_pnpn2_e_crs')

    ! Transfer the strong pressure facets to the lower H1 levels.  The middle
    ! and coarse levels use zero Dirichlet masks rather than the original fluid
    ! boundary objects.
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

    ! Finish the lower-level HSMG objects:
    ! - middle Schwarz smoother,
    ! - top-to-middle and middle-to-coarse interpolation operators,
    ! - coarse solver, either the dedicated Pn/Pn-2 CRS-like path or the
    !   alternative reusable HSMG path.
    call this%schwarz_mg%init(this%Xh_mg, this%dm_mg, this%gs_mg, &
         this%bclst_mg, c_Xh%msh)
    call this%interp_top_mg%init(Yh, this%Xh_mg)
    call this%interp_mg_crs%init(this%Xh_mg, this%Xh_crs)
    if (this%use_direct_coarse) then
      call this%crs_solver%init(this%c_Xh, this%dm_Xh, this%dm_crs, coarse_null_space)
    else if (trim(crs_solver_type) .eq. 'tamg') then
      call ax_helm_factory(this%ax, full_formulation = .false.)
      allocate(this%amg_solver)
      call this%amg_solver%init(this%ax, this%Xh_crs, this%c_crs, c_Xh%msh, &
           this%gs_crs, crs_tamg_lvls, this%bclst_crs, crs_tamg_itrs, &
           crs_tamg_cheby_degree)
    else
      call ax_helm_factory(this%ax, full_formulation = .false.)
      call precon_factory(this%pc_crs_solver, crs_pc_type)
      select type (pc => this%pc_crs_solver)
      type is (jacobi_t)
        call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
      type is (sx_jacobi_t)
        call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
      type is (device_jacobi_t)
        call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
      type is (ident_t)
      type is (device_ident_t)
      class default
        call neko_error('Unsupported PnPn-2 coarse-grid preconditioner.')
      end select
      call krylov_solver_factory(this%crs_ksp, this%dm_crs%size(), &
           trim(crs_solver_type), KSP_MAX_ITER, M = this%pc_crs_solver, &
           monitor = crs_monitor)
    end if
    this%crs_solver_initialized = .true.
  end subroutine hsmg_pnpn2_init

  !> Release pressure-space HSMG storage.
  !!
  !! The routine frees all owned work fields, interpolation operators, lower
  !! level coefficient objects, and the explicit coarse solver.  Pointers back
  !! to the surrounding fluid scheme are nullified but not deallocated, because
  !! ownership stays with the caller.
  subroutine hsmg_pnpn2_finalize(this)
    type(hsmg_pnpn2_t), intent(inout) :: this
    call this%free()
  end subroutine hsmg_pnpn2_finalize

  !> Release pressure-space HSMG storage (explicit).
  !!
  !! The routine frees all owned work fields, interpolation operators, lower
  !! level coefficient objects, and the explicit coarse solver.  Pointers back
  !! to the surrounding fluid scheme are nullified but not deallocated, because
  !! ownership stays with the caller.
  subroutine hsmg_pnpn2_free(this)
    class(hsmg_pnpn2_t), intent(inout) :: this

    call this%top_fdm_Xh%free()
    call this%schwarz_mg%free()
    call this%interp_top_mg%free()
    call this%interp_mg_crs%free()
    if (this%crs_solver_initialized) then
      if (this%use_direct_coarse) then
        call this%crs_solver%free()
      else
        if (allocated(this%crs_ksp)) then
          call this%crs_ksp%free()
          deallocate(this%crs_ksp)
        end if
        if (allocated(this%pc_crs_solver)) then
          call precon_destroy(this%pc_crs_solver)
          deallocate(this%pc_crs_solver)
        end if
        if (allocated(this%amg_solver)) then
          call this%amg_solver%free()
          deallocate(this%amg_solver)
        end if
        if (allocated(this%ax)) deallocate(this%ax)
      end if
      this%crs_solver_initialized = .false.
    end if
    call this%r_Yh%free()
    call this%w_Yh%free()
    call this%r_Xh%free()
    call this%z_Xh%free()
    call this%top_wt_Xh%free()
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

  !> Apply one pressure-space HSMG V-cycle.
  !!
  !! @param z Output correction on the local pressure grid.
  !! @param r Input residual on the local pressure grid.
  !! @param n Vector length of the local pressure grid.
  !!
  !! The algorithm is intentionally organised in the same order as Nek's H1
  !! multigrid cycle, with the important Pn/Pn-2 adaptation that the incoming
  !! vector is stored on `Y_h`:
  !! 1. embed the local `Y_h` residual into the interior of an `X_h` field;
  !! 2. apply the top Schwarz/FDM solve on that embedded field;
  !! 3. restrict the original residual to the middle H1 level;
  !! 4. apply the middle Schwarz/FDM solve;
  !! 5. weight the middle residual on the level boundary before coarse
  !!    restriction;
  !  ! 6. restrict to the coarse level and apply the dedicated coarse solve;
  !! 7. prolong the lower-level correction back through the hierarchy; and
  !! 8. combine the top `X_h` correction with the lower-level `Y_h`
  !!    contribution to form the final pressure correction.
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

    call profiler_start_region('HSMG_pnpn2_solve')
    call field_rzero(this%r_Xh)
    call field_rzero(this%z_Xh)

    ! Step 1:
    ! Embed the local `Y_h` residual into the interior points of an `X_h`
    ! field.  This is the current milestone analogue of Nek's fine-grid H1
    ! residual.  Boundary points remain zero because the local pressure space
    ! has no boundary nodes of its own.
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

    ! Step 2:
    ! Reproduce Nek's `local_solves_fdm` sequence on the embedded `X_h` field:
    ! extend the shell to the true faces, exchange those face values, remove the
    ! shell contribution, apply the lx1-sized tensor FDM, then fold the solved
    ! face values back to the shell before weighting.
    call profiler_start_region('HSMG_pnpn2_top_schwarz')
    call hsmg_pnpn2_pressure_face_ext(this%r_Xh%x, nx, nx, this%c_Xh%Xh%lz, nelv)
    call this%gs_Xh%op(this%r_Xh%x, this%dm_Xh%size(), GS_OP_ADD)
    call hsmg_pnpn2_pressure_face_add1si(this%r_Xh%x, -1.0_rp, nx, nx, &
         this%c_Xh%Xh%lz, nelv)
    call this%top_fdm_Xh%compute(this%z_Xh%x, this%r_Xh%x)
    call hsmg_pnpn2_pressure_shell_add_face(this%z_Xh%x, -1.0_rp, nx, nx, &
         this%c_Xh%Xh%lz, nelv)
    call this%gs_Xh%op(this%z_Xh%x, this%dm_Xh%size(), GS_OP_ADD)
    call hsmg_pnpn2_pressure_shell_add_face(this%z_Xh%x, 1.0_rp, nx, nx, &
         this%c_Xh%Xh%lz, nelv)
    call this%weight_top_pressure(this%z_Xh%x)
    call profiler_end_region('HSMG_pnpn2_top_schwarz')
    call copy(this%r_Yh%x, r, n)

    ! Step 3:
    ! Restrict the original local pressure residual to the middle H1 level and
    ! perform the middle-level gather-scatter.
    call profiler_start_region('HSMG_pnpn2_coarse_grid')
    call this%interp_top_mg%map(this%r_mg%x, this%r_Yh%x, nelv, this%Xh_mg)
    call this%gs_mg%op(this%r_mg%x, this%dm_mg%size(), GS_OP_ADD)

    ! Step 4:
    ! Apply the middle-level Schwarz/FDM smoother.
    call this%schwarz_mg%compute(this%e_mg%x, this%r_mg%x)

    ! Step 5:
    ! Apply Nek-style restriction weights only on the level boundary before the
    ! residual is transferred to the coarse grid.
    call this%weight_restriction_boundary(this%r_mg%x, this%c_mg)

    ! Step 6:
    ! Restrict to the coarse grid.  The dedicated coarse solver works on the
    ! unique coarse dofs internally, so this path does not do an additional
    ! coarse-grid gather-scatter beforehand.
    call this%interp_mg_crs%map(this%r_crs%x, this%r_mg%x, nelv, this%Xh_crs)
    if (this%use_direct_coarse) then
      call this%weight_coarse_mask(this%r_crs%x, this%c_crs)
      call profiler_start_region('HSMG_pnpn2_coarse_solve')
      call this%crs_solver%solve(this%e_crs%x, this%r_crs%x, this%c_crs)
      call profiler_end_region('HSMG_pnpn2_coarse_solve')
      call this%weight_coarse_mask(this%e_crs%x, this%c_crs)
    else
      call this%gs_crs%op(this%r_crs%x, this%dm_crs%size(), GS_OP_ADD)
      call this%bclst_crs%apply(this%r_crs)
      call field_rzero(this%e_crs)
      call profiler_start_region('HSMG_pnpn2_coarse_solve')
      if (allocated(this%amg_solver)) then
        call this%amg_solver%solve(this%e_crs%x, this%r_crs%x, this%dm_crs%size())
      else
        crs_info = this%crs_ksp%solve(this%ax, this%e_crs, this%r_crs%x, &
             this%dm_crs%size(), this%c_crs, this%bclst_crs, this%gs_crs, &
             this%crs_niter)
      end if
      call profiler_end_region('HSMG_pnpn2_coarse_solve')
      call this%bclst_crs%apply_scalar(this%e_crs%x, this%dm_crs%size())
    end if

    ! Step 7:
    ! Prolong the coarse correction back to the middle level and accumulate it
    ! with the middle Schwarz contribution.
    call this%interp_mg_crs%map(this%w_mg%x, this%e_crs%x, nelv, this%Xh_mg)
    call add2(this%e_mg%x, this%w_mg%x, this%dm_mg%size())

    ! Step 8:
    ! Prolong the lower-level correction to the local pressure grid.  The final
    ! result is the sum of:
    ! - the top `X_h` Schwarz correction extracted back to `Y_h`, and
    ! - the lower-level correction prolongated from the middle H1 level.
    call profiler_start_region('HSMG_pnpn2_prolongation')
    call this%interp_top_mg%map(this%w_Yh%x, this%e_mg%x, nelv, this%Yh)
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
    call add2(z, this%w_Yh%x, n)
    call profiler_end_region('HSMG_pnpn2_prolongation')
    call profiler_end_region('HSMG_pnpn2_coarse_grid')
    call profiler_end_region('HSMG_pnpn2_solve')
  end subroutine hsmg_pnpn2_solve

  !> Refresh HSMG coefficients.
  !!
  !! The current Pn/Pn-2 HSMG path treats the lower levels as pure H1 operators.
  !! Updating therefore means reapplying the H1 coefficient state to the middle
  !! and coarse coefficient objects after the surrounding fluid scheme has
  !! updated the fine-grid state.
  subroutine hsmg_pnpn2_update(this)
    class(hsmg_pnpn2_t), intent(inout) :: this
 
    call this%set_h1_coeffs(this%c_mg)
    call this%set_h1_coeffs(this%c_crs)
  end subroutine hsmg_pnpn2_update

  !> Force a coefficient object into the H1 form used by this preconditioner.
  !!
  !! @param coef Coefficient object to be rewritten in-place.
  !!
  !! HSMG pressure preconditioning on these levels uses the stiffness part only:
  !! `h1 = 1`, `h2 = 0`, and `ifh2 = .false.`.
  subroutine hsmg_pnpn2_set_h1_coeffs(this, coef)
    class(hsmg_pnpn2_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: n

    n = coef%dof%size()
    coef%h1(1:n,1,1,1) = 1.0_rp
    coef%h2(1:n,1,1,1) = 0.0_rp
    coef%ifh2 = .false.
  end subroutine hsmg_pnpn2_set_h1_coeffs

  !> Build Nek's `do_weight_op` multiplicity weights on the embedded pressure shell.
  !!
  !! Nek initialises these weights by marking the pressure-space boundary nodes
  !! inside the mesh-1 field, copying that shell to the true element faces,
  !! gather-scattering the face values, and then accumulating the resulting
  !! multiplicities back onto the shell.  The reciprocal of that count is then
  !! used after the top local FDM solve.
  !!
  !! This routine mirrors that setup on Neko's embedded `X_h` field.  The final
  !! field is `1` everywhere except on the pressure shell at indices
  !! `2`/`nx-1`, where it stores the reciprocal multiplicity.
  subroutine hsmg_pnpn2_build_top_weight(this)
    class(hsmg_pnpn2_t), intent(inout) :: this
    integer :: e, i, j, k
    integer :: nx, ny, nz, nelv

    nx = this%c_Xh%Xh%lx
    ny = this%c_Xh%Xh%ly
    nz = this%c_Xh%Xh%lz
    nelv = this%c_Xh%msh%nelv

    this%top_wt_Xh%x = 1.0_rp
    call field_rzero(this%r_Xh)

    ! Mark only the embedded pressure-space boundary shell.
    if (this%c_Xh%msh%gdim .eq. 2) then
      do e = 1, nelv
        do j = 2, ny - 1
          this%r_Xh%x(2,j,1,e) = 1.0_rp
          this%r_Xh%x(nx - 1,j,1,e) = 1.0_rp
        end do
        do i = 2, nx - 1
          this%r_Xh%x(i,2,1,e) = 1.0_rp
          this%r_Xh%x(i,ny - 1,1,e) = 1.0_rp
        end do
      end do
    else
      do e = 1, nelv
        do k = 2, nz - 1
          do j = 2, ny - 1
            this%r_Xh%x(2,j,k,e) = 1.0_rp
            this%r_Xh%x(nx - 1,j,k,e) = 1.0_rp
          end do
        end do
        do k = 2, nz - 1
          do i = 2, nx - 1
            this%r_Xh%x(i,2,k,e) = 1.0_rp
            this%r_Xh%x(i,ny - 1,k,e) = 1.0_rp
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            this%r_Xh%x(i,j,2,e) = 1.0_rp
            this%r_Xh%x(i,j,nz - 1,e) = 1.0_rp
          end do
        end do
      end do
    end if

    ! Mirror Nek's init_weight_op:
    ! 1. copy the shell to the true element faces,
    ! 2. gather-scatter on those shared face nodes,
    ! 3. subtract the local shell contribution from the face values, and
    ! 4. accumulate the resulting multiplicity back onto the shell.
    call hsmg_pnpn2_pressure_face_ext(this%r_Xh%x, nx, ny, nz, nelv)
    call this%gs_Xh%op(this%r_Xh%x, this%dm_Xh%size(), GS_OP_ADD)
    call hsmg_pnpn2_pressure_face_add1si(this%r_Xh%x, -1.0_rp, nx, ny, nz, nelv)
    call hsmg_pnpn2_pressure_shell_add_face(this%r_Xh%x, 1.0_rp, nx, ny, nz, nelv)

    ! Turn shell multiplicities into reciprocal weights and leave every other
    ! point as unity so a plain pointwise multiplication reproduces do_weight_op.
    if (this%c_Xh%msh%gdim .eq. 2) then
      do e = 1, nelv
        do j = 2, ny - 1
          this%top_wt_Xh%x(2,j,1,e) = 1.0_rp / this%r_Xh%x(2,j,1,e)
          this%top_wt_Xh%x(nx - 1,j,1,e) = 1.0_rp / this%r_Xh%x(nx - 1,j,1,e)
        end do
        do i = 2, nx - 1
          this%top_wt_Xh%x(i,2,1,e) = 1.0_rp / this%r_Xh%x(i,2,1,e)
          this%top_wt_Xh%x(i,ny - 1,1,e) = 1.0_rp / this%r_Xh%x(i,ny - 1,1,e)
        end do
      end do
    else
      do e = 1, nelv
        do k = 2, nz - 1
          do j = 2, ny - 1
            this%top_wt_Xh%x(2,j,k,e) = 1.0_rp / this%r_Xh%x(2,j,k,e)
            this%top_wt_Xh%x(nx - 1,j,k,e) = 1.0_rp / this%r_Xh%x(nx - 1,j,k,e)
          end do
        end do
        do k = 2, nz - 1
          do i = 2, nx - 1
            this%top_wt_Xh%x(i,2,k,e) = 1.0_rp / this%r_Xh%x(i,2,k,e)
            this%top_wt_Xh%x(i,ny - 1,k,e) = 1.0_rp / this%r_Xh%x(i,ny - 1,k,e)
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            this%top_wt_Xh%x(i,j,2,e) = 1.0_rp / this%r_Xh%x(i,j,2,e)
            this%top_wt_Xh%x(i,j,nz - 1,e) = 1.0_rp / this%r_Xh%x(i,j,nz - 1,e)
          end do
        end do
      end do
    end if

    call field_rzero(this%r_Xh)
  end subroutine hsmg_pnpn2_build_top_weight

  !> Apply Nek's top-level pressure-shell weights after the local FDM solve.
  !!
  !! @param u Embedded `X_h` correction to be weighted in-place.
  subroutine hsmg_pnpn2_weight_top_pressure(this, u)
    class(hsmg_pnpn2_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: u(:,:,:,:)

    u = u * this%top_wt_Xh%x
  end subroutine hsmg_pnpn2_weight_top_pressure

  !> Apply Nek's coarse `mg_mask` weighting in-place.
  !!
  !! @param u Coarse residual or correction field to be weighted in-place.
  subroutine hsmg_pnpn2_weight_coarse_mask(this, u, coef)
    class(hsmg_pnpn2_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: u(:,:,:,:)

    call this%bclst_crs%apply_scalar(u, this%dm_crs%size())
  end subroutine hsmg_pnpn2_weight_coarse_mask

  !> Apply Nek-style restriction weights on boundary points only.
  !!
  !! @param u Residual field to be weighted in-place.
  !! @param coef Coefficients that provide the point multiplicity field.
  !!
  !! Nek's `hsmg_do_wt` multiplies only the level-boundary points by reciprocal
  !! sharing counts before fine-to-coarse restriction.  Interior points are left
  !! unchanged.  On the Neko H1 levels the coefficient multiplicity field is
  !! already the level-local reciprocal DSS count, so restricting the action to
  !! boundary points reproduces the same data used by Nek's `mg_rstr_wt`.
  subroutine hsmg_pnpn2_weight_restriction_boundary(this, u, coef)
    class(hsmg_pnpn2_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: u(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv)
    integer :: e, i, j, k
    integer :: nx, ny, nz

    nx = coef%Xh%lx
    ny = coef%Xh%ly
    nz = coef%Xh%lz

    if (coef%msh%gdim .eq. 2) then
      do e = 1, coef%msh%nelv
        do j = 1, ny
          u(1,j,1,e) = u(1,j,1,e) * coef%mult(1,j,1,e)
          u(nx,j,1,e) = u(nx,j,1,e) * coef%mult(nx,j,1,e)
        end do
        do i = 2, nx - 1
          u(i,1,1,e) = u(i,1,1,e) * coef%mult(i,1,1,e)
          u(i,ny,1,e) = u(i,ny,1,e) * coef%mult(i,ny,1,e)
        end do
      end do
    else
      do e = 1, coef%msh%nelv
        do k = 1, nz
          do j = 1, ny
            u(1,j,k,e) = u(1,j,k,e) * coef%mult(1,j,k,e)
            u(nx,j,k,e) = u(nx,j,k,e) * coef%mult(nx,j,k,e)
          end do
        end do
        do k = 1, nz
          do i = 2, nx - 1
            u(i,1,k,e) = u(i,1,k,e) * coef%mult(i,1,k,e)
            u(i,ny,k,e) = u(i,ny,k,e) * coef%mult(i,ny,k,e)
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            u(i,j,1,e) = u(i,j,1,e) * coef%mult(i,j,1,e)
            u(i,j,nz,e) = u(i,j,nz,e) * coef%mult(i,j,nz,e)
          end do
        end do
      end do
    end if
  end subroutine hsmg_pnpn2_weight_restriction_boundary

  !> Copy the embedded pressure shell to the true element faces.
  !!
  !! This is the mesh-1 analogue of Nek's `dface_ext` used in `init_weight_op`.
  subroutine hsmg_pnpn2_pressure_face_ext(u, nx, ny, nz, nelv)
    integer, intent(in) :: nx, ny, nz, nelv
    real(kind=rp), intent(inout) :: u(nx, ny, nz, nelv)
    integer :: e, i, j, k

    if (nz .eq. 1) then
      do e = 1, nelv
        do i = 2, nx - 1
          u(i,1,1,e) = u(i,2,1,e)
          u(i,ny,1,e) = u(i,ny - 1,1,e)
        end do
        do j = 2, ny - 1
          u(1,j,1,e) = u(2,j,1,e)
          u(nx,j,1,e) = u(nx - 1,j,1,e)
        end do
      end do
    else
      do e = 1, nelv
        do k = 2, nz - 1
          do i = 2, nx - 1
            u(i,1,k,e) = u(i,2,k,e)
            u(i,ny,k,e) = u(i,ny - 1,k,e)
          end do
        end do
        do k = 2, nz - 1
          do j = 2, ny - 1
            u(1,j,k,e) = u(2,j,k,e)
            u(nx,j,k,e) = u(nx - 1,j,k,e)
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            u(i,j,1,e) = u(i,j,2,e)
            u(i,j,nz,e) = u(i,j,nz - 1,e)
          end do
        end do
      end do
    end if
  end subroutine hsmg_pnpn2_pressure_face_ext

  !> Add a scaled shell contribution to the true element faces.
  !!
  !! This mirrors Nek's `dface_add1si` for the embedded pressure shell.
  subroutine hsmg_pnpn2_pressure_face_add1si(u, c, nx, ny, nz, nelv)
    integer, intent(in) :: nx, ny, nz, nelv
    real(kind=rp), intent(in) :: c
    real(kind=rp), intent(inout) :: u(nx, ny, nz, nelv)
    integer :: e, i, j, k

    if (nz .eq. 1) then
      do e = 1, nelv
        do i = 2, nx - 1
          u(i,1,1,e) = u(i,1,1,e) + c * u(i,2,1,e)
          u(i,ny,1,e) = u(i,ny,1,e) + c * u(i,ny - 1,1,e)
        end do
        do j = 2, ny - 1
          u(1,j,1,e) = u(1,j,1,e) + c * u(2,j,1,e)
          u(nx,j,1,e) = u(nx,j,1,e) + c * u(nx - 1,j,1,e)
        end do
      end do
    else
      do e = 1, nelv
        do k = 2, nz - 1
          do i = 2, nx - 1
            u(i,1,k,e) = u(i,1,k,e) + c * u(i,2,k,e)
            u(i,ny,k,e) = u(i,ny,k,e) + c * u(i,ny - 1,k,e)
          end do
        end do
        do k = 2, nz - 1
          do j = 2, ny - 1
            u(1,j,k,e) = u(1,j,k,e) + c * u(2,j,k,e)
            u(nx,j,k,e) = u(nx,j,k,e) + c * u(nx - 1,j,k,e)
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            u(i,j,1,e) = u(i,j,1,e) + c * u(i,j,2,e)
            u(i,j,nz,e) = u(i,j,nz,e) + c * u(i,j,nz - 1,e)
          end do
        end do
      end do
    end if
  end subroutine hsmg_pnpn2_pressure_face_add1si

  !> Add a scaled face contribution back to the embedded pressure shell.
  !!
  !! This is the embedded-shell equivalent of Nek's `s_face_to_int`.
  subroutine hsmg_pnpn2_pressure_shell_add_face(u, c, nx, ny, nz, nelv)
    integer, intent(in) :: nx, ny, nz, nelv
    real(kind=rp), intent(in) :: c
    real(kind=rp), intent(inout) :: u(nx, ny, nz, nelv)
    integer :: e, i, j, k

    if (nz .eq. 1) then
      do e = 1, nelv
        do i = 2, nx - 1
          u(i,2,1,e) = u(i,2,1,e) + c * u(i,1,1,e)
          u(i,ny - 1,1,e) = u(i,ny - 1,1,e) + c * u(i,ny,1,e)
        end do
        do j = 2, ny - 1
          u(2,j,1,e) = u(2,j,1,e) + c * u(1,j,1,e)
          u(nx - 1,j,1,e) = u(nx - 1,j,1,e) + c * u(nx,j,1,e)
        end do
      end do
    else
      do e = 1, nelv
        do k = 2, nz - 1
          do i = 2, nx - 1
            u(i,2,k,e) = u(i,2,k,e) + c * u(i,1,k,e)
            u(i,ny - 1,k,e) = u(i,ny - 1,k,e) + c * u(i,ny,k,e)
          end do
        end do
        do k = 2, nz - 1
          do j = 2, ny - 1
            u(2,j,k,e) = u(2,j,k,e) + c * u(1,j,k,e)
            u(nx - 1,j,k,e) = u(nx - 1,j,k,e) + c * u(nx,j,k,e)
          end do
        end do
        do j = 2, ny - 1
          do i = 2, nx - 1
            u(i,j,2,e) = u(i,j,2,e) + c * u(i,j,1,e)
            u(i,j,nz - 1,e) = u(i,j,nz - 1,e) + c * u(i,j,nz,e)
          end do
        end do
      end do
    end if
  end subroutine hsmg_pnpn2_pressure_shell_add_face

  !> Return Nek's middle HSMG point count for the current velocity GLL count.
  !!
  !! @param lx1 Number of GLL points on the fine `X_h` grid.
  !! @return Number of GLL points used on the middle HSMG level.
  !!
  !! Nek tabulates the middle-level polynomial count `mg_nx`; Neko stores point
  !! counts, so this helper performs the same lookup and converts the result to
  !! the corresponding GLL point count.
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
