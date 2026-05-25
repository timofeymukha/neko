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
!> Implements type `fluid_pnpn2_t`.
module fluid_pnpn2
  use, intrinsic :: iso_fortran_env, only : error_unit
  use ax_product, only : ax_t, ax_helm_factory
  use advection, only : advection_t, advection_factory
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use checkpoint, only : chkp_t
  use coefs, only : coef_t
  use comm, only : NEKO_COMM
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_math, only : field_add2, field_cmult, field_rzero
  use field_series, only : field_series_t
  use fluid_aux, only : fluid_step_info
  use pnpn2_mixed_ops, only : pnpn2_mixed_ops_t
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use hsmg_pnpn2, only : hsmg_pnpn2_t
  use hsmg, only : hsmg_t
  use identity, only : ident_t
  use inflow, only : inflow_t
  use interpolation, only : interpolator_t
  use jacobi, only : jacobi_t
  use json_module, only : json_core, json_file, json_value
  use json_utils, only : json_get, json_get_or_default, json_get_or_lookup, &
       json_get_or_lookup_or_default, json_extract_item
  use krylov, only : ksp_monitor_t
  use logger, only : neko_log, LOG_SIZE
  use math, only : add2, add2s2, col2, cmult2, copy, glsc2, glsum, rzero
  use mathops, only : opadd2cm
  use mesh, only : mesh_t
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_LOGICAL, &
       MPI_LOR, MPI_MAX
  use mxm_wrapper, only : mxm
  use neko_config, only : NEKO_BCKND_DEVICE
  use no_slip, only : no_slip_t
  use num_types, only : i8, rp, xp
  use operators, only : ortho, rotate_cyc
  use pnpn2_prs_ax, only : pnpn2_prs_ax_t, pnpn2_prs_ax_clear, &
       pnpn2_prs_ax_init
  use phmg, only : phmg_t
  use precon, only : pc_t, precon_destroy, precon_factory
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use registry, only : neko_registry
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t, &
       rhs_maker_oifs_t, rhs_maker_bdf_fctry, rhs_maker_ext_fctry, &
       rhs_maker_oifs_fctry
  use scratch_registry, only : neko_scratch_registry
  use space, only : GL, GLL, space_t
  use tensor, only : tnsr3d
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t
  use utils, only : neko_error, neko_type_error
  use zero_dirichlet, only : zero_dirichlet_t
  implicit none
  private

  character(len=25), private :: FLUID_PNPN2_KNOWN_BCS(4) = &
       [character(len=25) :: "no_slip", "outflow", "velocity_inlet", &
       "velocity_value"]

  !> Pressure preconditioner wrapper using Nek-style \f$ Y_h \to X_h \to Y_h \f$ mapping.
  type, private, extends(pc_t) :: pnpn2_prs_precon_t
     !> Wrapped preconditioner acting on the velocity space.
     class(pc_t), allocatable :: pc_Xh
     !> Mixed-space operators providing Nek-style `MAP21E/MAP12` tensors.
     type(pnpn2_mixed_ops_t), pointer :: mixed_ops => null()
     !> Velocity coefficients reused as the mesh-1 H1 operator state.
     type(coef_t), pointer :: c_Xh => null()
     !> Velocity dofmap.
     type(dofmap_t), pointer :: dm_Xh => null()
     !> Velocity gather-scatter.
     type(gs_t), pointer :: gs_Xh => null()
     !> Elemental global point count on `X_h` for Nek-style mean removal.
     integer(kind=i8) :: glb_Xh_points = 0_i8
     !> Whether a strong pressure Dirichlet boundary exists.
     logical :: prs_dirichlet = .false.
     !> Velocity-space image of strong pressure boundaries.
     type(bc_list_t) :: bclst_Xh
     !> Velocity-space image of the pressure residual.
     type(field_t) :: r_Xh
     !> Velocity-space preconditioned correction.
     type(field_t) :: z_Xh
   contains
     !> Constructor.
     procedure, pass(this) :: init => pnpn2_prs_precon_init
     !> Destructor.
     procedure, pass(this) :: free => pnpn2_prs_precon_free
     !> Apply the wrapped pressure preconditioner.
     procedure, pass(this) :: solve => pnpn2_prs_precon_solve
     !> Update the wrapped pressure preconditioner.
     procedure, pass(this) :: update => pnpn2_prs_precon_update
    !> Set the dedicated H1 coefficients on `X_h`.
    procedure, pass(this) :: set_h1 => pnpn2_prs_precon_set_h1
    !> Apply Nek `MAP21E`.
    procedure, pass(this) :: map21e => pnpn2_prs_precon_map21e
    !> Apply Nek `MAP12`.
    procedure, pass(this) :: map12 => pnpn2_prs_precon_map12
  end type pnpn2_prs_precon_t

  !> Dedicated Nek-style local GMRES for the PnPn-2 pressure solve on \f$ Y_h \f$.
  type, private :: pnpn2_prs_gmres_t
     !> Restart length.
     integer :: lgmres = 30
     !> Maximum number of iterations.
     integer :: max_iter = 0
     !> Absolute convergence tolerance.
     real(kind=rp) :: abs_tol = 0.0_rp
     !> Enable residual logging.
     logical :: monitor = .false.
     !> Left split \f$ L = \sqrt{BM2INV} \f$.
     real(kind=rp), allocatable :: ml(:)
     !> Right split \f$ U = \sqrt{BM2} \f$.
     real(kind=rp), allocatable :: mu(:)
     !> Unscaled residual.
     real(kind=rp), allocatable :: r(:)
     !> Split residual \f$ \hat{r} = L r \f$.
     real(kind=rp), allocatable :: r_hat(:)
     !> Generic work vector on \f$ Y_h \f$.
     real(kind=rp), allocatable :: w(:)
     !> Split operator image \f$ \hat{w} = L A z \f$.
     real(kind=rp), allocatable :: w_hat(:)
     !> Solver-local right-hand side.
     real(kind=rp), allocatable :: rhs(:)
     !> Preconditioned Krylov directions on \f$ Y_h \f$.
     real(kind=rp), allocatable :: z(:,:)
     !> Arnoldi basis in the split pressure space.
     real(kind=rp), allocatable :: v(:,:)
     !> Hessenberg matrix.
     real(kind=xp), allocatable :: h(:,:)
     !> Givens-rotated residual coefficients.
     real(kind=xp), allocatable :: gamma(:)
     !> Cosines of Givens rotations.
     real(kind=xp), allocatable :: givens_c(:)
     !> Sines of Givens rotations.
     real(kind=xp), allocatable :: givens_s(:)
     !> Back-substitution coefficients.
     real(kind=xp), allocatable :: y(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => pnpn2_prs_gmres_init
     !> Destructor.
     procedure, pass(this) :: free => pnpn2_prs_gmres_free
     !> Solve the local pressure system with Nek-style split GMRES semantics.
     procedure, pass(this) :: solve => pnpn2_prs_gmres_solve
  end type pnpn2_prs_gmres_t

  !> Bare-bones CPU PnPn-2 pressure-correction scheme.
  type, public, extends(fluid_scheme_incompressible_t) :: fluid_pnpn2_t
     !> Pressure function space \f$ Y_h = P_{N-2} \f$.
     type(space_t) :: Yh
     !> Pressure dofmap on \f$ Y_h \f$.
     type(dofmap_t) :: dm_Yh
     !> Element-local no-op gather-scatter used only to satisfy generic APIs.
     type(gs_t) :: gs_prs
     !> Coefficients on \f$ Y_h \f$.
     type(coef_t) :: c_Yh
     !> Interpolator between \f$ X_h \f$ and \f$ Y_h \f$.
     type(interpolator_t) :: prs_interp
     !> Non-unique global number of lower-order pressure points.
     integer(kind=i8) :: glb_prs_points = 0_i8
     !> Solver-authoritative internal pressure on \f$ Y_h \f$.
     type(field_t) :: p_Yh
     !> Pressure right-hand side.
     type(field_t) :: p_res
     !> Pressure increment.
     type(field_t) :: dp
     !> Extrapolated pressure used by the predictor and pressure update.
     type(field_t) :: p_ext
     !> Velocity right-hand side.
     type(field_t) :: u_res, v_res, w_res
     !> Velocity increment.
     type(field_t) :: du, dv, dw
     !> Velocity Helmholtz operator.
     class(ax_t), allocatable :: Ax_vel
     !> Pressure projection operator.
     class(ax_t), allocatable :: Ax_prs
     !> Pressure projection history.
     type(projection_t) :: proj_prs
     !> Velocity projection history.
     type(projection_vel_t) :: proj_vel
     !> BDF history maker.
     class(rhs_maker_bdf_t), allocatable :: makebdf
     !> Pressure extrapolation history on \f$ Y_h \f$.
     type(field_series_t) :: plag
     !> Mixed-space pressure operators on \f$ X_h \leftrightarrow Y_h \f$.
     type(pnpn2_mixed_ops_t) :: mixed_ops
     !> Pressure-space boundary list for the linear solver.
     type(bc_list_t) :: bclst_dp
     !> Dummy strong-velocity marker for the residual.
     type(zero_dirichlet_t) :: bc_vel_res
     !> Velocity residual boundary list.
     type(bc_list_t) :: bclst_vel_res
     !> Dedicated local GMRES used for the pressure solve on \f$ Y_h \f$.
     type(pnpn2_prs_gmres_t) :: prs_gmres
     !> Dummy x-velocity increment Dirichlet marker.
     type(zero_dirichlet_t) :: bc_du
     !> X-velocity increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_du
     !> Dummy y-velocity increment Dirichlet marker.
     type(zero_dirichlet_t) :: bc_dv
     !> Y-velocity increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_dv
     !> Dummy z-velocity increment Dirichlet marker.
     type(zero_dirichlet_t) :: bc_dw
     !> Z-velocity increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_dw
     !> Whether a strong pressure Dirichlet boundary exists.
     logical :: prs_dirichlet = .false.
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_pnpn2_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_pnpn2_free
     !> Advance one time step.
     procedure, pass(this) :: step => fluid_pnpn2_step
     !> Restart from a checkpoint.
     procedure, pass(this) :: restart => fluid_pnpn2_restart
     !> Set up the supported boundary-condition configuration.
     procedure, pass(this) :: setup_bcs => fluid_pnpn2_setup_bcs
     !> Synchronize the public \f$ X_h \f$ pressure from the authoritative \f$ Y_h \f$ state.
     procedure, pass(this) :: sync_p_public => fluid_pnpn2_sync_p_public
     !> Synchronize the authoritative \f$ Y_h \f$ pressure from the public \f$ X_h \f$ field.
     procedure, pass(this) :: sync_p_from_public => fluid_pnpn2_sync_p_from_public
  end type fluid_pnpn2_t

  ! Keep the temporary milestone-1 advection state outside fluid_pnpn2_t.
  ! The current branch's pressure-preconditioner setup is sensitive to changes
  ! in the extended solver type layout.
  class(advection_t), allocatable, save :: pnpn2_adv
  logical, save :: pnpn2_oifs = .false.
  type(field_t), save :: pnpn2_abx1, pnpn2_aby1, pnpn2_abz1
  type(field_t), save :: pnpn2_abx2, pnpn2_aby2, pnpn2_abz2
  type(field_t), save :: pnpn2_advx, pnpn2_advy, pnpn2_advz
  class(rhs_maker_ext_t), allocatable, save :: pnpn2_makeabf
  class(rhs_maker_oifs_t), allocatable, save :: pnpn2_makeoifs

contains

  !> Initialize the wrapped pressure preconditioner on \f$ X_h \f$.
  subroutine pnpn2_prs_precon_init(this, precon_type, precon_params, mixed_ops, &
     c_Xh, dm_Xh, gs_Xh, bcs_prs_Xh, prs_dirichlet)
   class(pnpn2_prs_precon_t), intent(inout), target :: this
   character(len=*), intent(in) :: precon_type
   type(json_file), intent(inout) :: precon_params
   type(pnpn2_mixed_ops_t), target, intent(inout) :: mixed_ops
   type(coef_t), target, intent(inout) :: c_Xh
   type(dofmap_t), target, intent(inout) :: dm_Xh
   type(gs_t), target, intent(inout) :: gs_Xh
   type(bc_list_t), target, intent(inout) :: bcs_prs_Xh
   logical, intent(in) :: prs_dirichlet
   integer :: i

   call this%free()

   this%mixed_ops => mixed_ops
   this%c_Xh => c_Xh
   this%dm_Xh => dm_Xh
   this%gs_Xh => gs_Xh
   this%prs_dirichlet = prs_dirichlet
   this%glb_Xh_points = int(dm_Xh%msh%glb_nelv, i8) * int(dm_Xh%Xh%lxyz, i8)

   call this%bclst_Xh%init()
   do i = 1, bcs_prs_Xh%size()
     if (bcs_prs_Xh%strong(i)) then
       call this%bclst_Xh%append(bcs_prs_Xh%get(i))
     end if
   end do
   call this%r_Xh%init(dm_Xh, 'pnpn2_prs_r_Xh')
   call this%z_Xh%init(dm_Xh, 'pnpn2_prs_z_Xh')

   call precon_factory(this%pc_Xh, precon_type)

   select type (pcp => this%pc_Xh)
   type is (jacobi_t)
      call pcp%init(c_Xh, dm_Xh, gs_Xh)
   type is (hsmg_t)
      call pcp%init(c_Xh, this%bclst_Xh, precon_params)
   type is (phmg_t)
      call pcp%init(c_Xh, this%bclst_Xh, precon_params)
   type is (ident_t)
      continue
   class default
      call neko_error('Unsupported PnPn-2 pressure preconditioner on X_h.')
   end select
  end subroutine pnpn2_prs_precon_init

  !> Free the wrapped pressure preconditioner.
  subroutine pnpn2_prs_precon_free(this)
   class(pnpn2_prs_precon_t), intent(inout) :: this

   call this%bclst_Xh%free()
   call this%r_Xh%free()
   call this%z_Xh%free()

   if (allocated(this%pc_Xh)) then
     call precon_destroy(this%pc_Xh)
     deallocate(this%pc_Xh)
   end if

   nullify(this%mixed_ops)
   nullify(this%c_Xh)
   nullify(this%dm_Xh)
   nullify(this%gs_Xh)
   this%glb_Xh_points = 0_i8
   this%prs_dirichlet = .false.
  end subroutine pnpn2_prs_precon_free

  !> Apply the dedicated mesh-1 H1 coefficients used by the pressure preconditioner.
  subroutine pnpn2_prs_precon_set_h1(this)
   class(pnpn2_prs_precon_t), intent(inout) :: this
   integer :: i, n

   n = this%dm_Xh%size()
   do concurrent (i = 1:n)
    this%c_Xh%h1(i,1,1,1) = 1.0_rp
    this%c_Xh%h2(i,1,1,1) = 0.0_rp
   end do
   this%c_Xh%ifh2 = .false.
  end subroutine pnpn2_prs_precon_set_h1

  !> Apply Nek `MAP21E` from the local pressure mesh to the velocity mesh.
  subroutine pnpn2_prs_precon_map21e(this, y, x)
   class(pnpn2_prs_precon_t), intent(inout) :: this
   real(kind=rp), intent(inout) :: y(this%mixed_ops%Xh%lx, this%mixed_ops%Xh%ly, &
       this%mixed_ops%Xh%lz, this%dm_Xh%msh%nelv)
   real(kind=rp), intent(in) :: x(this%mixed_ops%Yh%lx, this%mixed_ops%Yh%ly, &
       this%mixed_ops%Yh%lz, this%dm_Xh%msh%nelv)

   call tnsr3d(y, this%mixed_ops%Xh%lx, x, this%mixed_ops%Yh%lx, &
       this%mixed_ops%i12t, this%mixed_ops%i12, this%mixed_ops%i12, &
       this%dm_Xh%msh%nelv)
  end subroutine pnpn2_prs_precon_map21e

  !> Apply Nek `MAP12` from the velocity mesh to the local pressure mesh.
  subroutine pnpn2_prs_precon_map12(this, y, x)
   class(pnpn2_prs_precon_t), intent(inout) :: this
   real(kind=rp), intent(inout) :: y(this%mixed_ops%Yh%lx, this%mixed_ops%Yh%ly, &
       this%mixed_ops%Yh%lz, this%dm_Xh%msh%nelv)
   real(kind=rp), intent(in) :: x(this%mixed_ops%Xh%lx, this%mixed_ops%Xh%ly, &
       this%mixed_ops%Xh%lz, this%dm_Xh%msh%nelv)

   call tnsr3d(y, this%mixed_ops%Yh%lx, x, this%mixed_ops%Xh%lx, &
       this%mixed_ops%i12, this%mixed_ops%i12t, this%mixed_ops%i12t, &
       this%dm_Xh%msh%nelv)
  end subroutine pnpn2_prs_precon_map12

  !> Apply `MAP21E`, a mesh-1 H1 preconditioner, and `MAP12`.
  subroutine pnpn2_prs_precon_solve(this, z, r, n)
   class(pnpn2_prs_precon_t), intent(inout) :: this
   integer, intent(in) :: n
   real(kind=rp), dimension(n), intent(inout) :: z
   real(kind=rp), dimension(n), intent(inout) :: r

   call this%map21e(this%r_Xh%x, r)
   if (.not. this%prs_dirichlet) then
     call ortho(this%r_Xh%x, this%glb_Xh_points, this%dm_Xh%size())
   end if
   call this%pc_Xh%solve(this%z_Xh%x, this%r_Xh%x, this%dm_Xh%size())
   call this%map12(z, this%z_Xh%x)
  end subroutine pnpn2_prs_precon_solve

  !> Update the wrapped \f$ X_h \f$ preconditioner.
  subroutine pnpn2_prs_precon_update(this)
   class(pnpn2_prs_precon_t), intent(inout) :: this

   call this%set_h1()
   call this%pc_Xh%update()
  end subroutine pnpn2_prs_precon_update

  !> Factory routine for supported PnPn-2 velocity boundary conditions.
  subroutine pnpn2_velocity_bc_factory(object, json, coef)
    class(bc_t), pointer, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    character(len=:), allocatable :: type
    character(len=:), allocatable :: default_name
    integer, allocatable :: zone_indices(:)
    character(len=64) :: buf
    integer :: i

    if (associated(object)) then
      call object%free()
      nullify(object)
    end if

    call json_get(json, 'type', type)

    select case (trim(type))
    case ('velocity_inlet', 'velocity_value')
      allocate(inflow_t :: object)
    case ('no_slip')
      allocate(no_slip_t :: object)
    case ('outflow')
      if (allocated(type)) then
        deallocate(type)
      end if
      return
    case default
      call neko_type_error('fluid_pnpn2 boundary conditions', type, &
           FLUID_PNPN2_KNOWN_BCS)
    end select

    call json_get_or_lookup(json, 'zone_indices', zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
      call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do

    write(buf, '("velocity_bc_", I0)') zone_indices(1)
    default_name = trim(buf)
    call json_get_or_default(json, 'name', object%name, default_name)
    object%zone_indices = zone_indices
    call object%finalize()

    if (allocated(default_name)) then
      deallocate(default_name)
    end if
    if (allocated(type)) then
      deallocate(type)
    end if
    if (allocated(zone_indices)) then
      deallocate(zone_indices)
    end if
  end subroutine pnpn2_velocity_bc_factory

  !> Factory routine for supported PnPn-2 pressure boundary conditions on `X_h`.
  subroutine pnpn2_pressure_bc_factory(object, json, coef)
    class(bc_t), pointer, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    character(len=:), allocatable :: type
    character(len=:), allocatable :: default_name
    integer, allocatable :: zone_indices(:)
    character(len=64) :: buf
    integer :: i

    if (associated(object)) then
      call object%free()
      nullify(object)
    end if

    call json_get(json, 'type', type)

    select case (trim(type))
    case ('outflow')
      allocate(zero_dirichlet_t :: object)
    case ('no_slip', 'velocity_inlet', 'velocity_value')
      if (allocated(type)) then
        deallocate(type)
      end if
      return
    case default
      call neko_type_error('fluid_pnpn2 boundary conditions', type, &
           FLUID_PNPN2_KNOWN_BCS)
    end select

    call json_get_or_lookup(json, 'zone_indices', zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
      call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do

    write(buf, '("pressure_bc_", I0)') zone_indices(1)
    default_name = trim(buf)
    call json_get_or_default(json, 'name', object%name, default_name)
    object%zone_indices = zone_indices
    call object%finalize()

    if (allocated(default_name)) then
      deallocate(default_name)
    end if
    if (allocated(type)) then
      deallocate(type)
    end if
    if (allocated(zone_indices)) then
      deallocate(zone_indices)
    end if
  end subroutine pnpn2_pressure_bc_factory

  !> Initialize the dedicated pressure GMRES work arrays and BM2 splits.
  subroutine pnpn2_prs_gmres_init(this, n, max_iter, abs_tol, monitor, b, binv)
   class(pnpn2_prs_gmres_t), intent(inout) :: this
   integer, intent(in) :: n
   integer, intent(in) :: max_iter
   real(kind=rp), intent(in) :: abs_tol
   logical, intent(in) :: monitor
   real(kind=rp), intent(in) :: b(n)
   real(kind=rp), intent(in) :: binv(n)
   integer :: i

   call this%free()

   this%max_iter = max_iter
   this%abs_tol = abs_tol
   this%monitor = monitor

   allocate(this%ml(n))
   allocate(this%mu(n))
   allocate(this%r(n))
   allocate(this%r_hat(n))
   allocate(this%w(n))
   allocate(this%w_hat(n))
   allocate(this%rhs(n))
   allocate(this%z(n, this%lgmres))
   allocate(this%v(n, this%lgmres + 1))
   allocate(this%h(this%lgmres + 1, this%lgmres))
   allocate(this%gamma(this%lgmres + 1))
   allocate(this%givens_c(this%lgmres))
   allocate(this%givens_s(this%lgmres))
   allocate(this%y(this%lgmres))

   do concurrent (i = 1:n)
     this%ml(i) = sqrt(binv(i))
     this%mu(i) = sqrt(b(i))
   end do
  end subroutine pnpn2_prs_gmres_init

  !> Free the dedicated pressure GMRES work arrays.
  subroutine pnpn2_prs_gmres_free(this)
   class(pnpn2_prs_gmres_t), intent(inout) :: this

   if (allocated(this%ml)) deallocate(this%ml)
   if (allocated(this%mu)) deallocate(this%mu)
   if (allocated(this%r)) deallocate(this%r)
   if (allocated(this%r_hat)) deallocate(this%r_hat)
   if (allocated(this%w)) deallocate(this%w)
   if (allocated(this%w_hat)) deallocate(this%w_hat)
   if (allocated(this%rhs)) deallocate(this%rhs)
   if (allocated(this%z)) deallocate(this%z)
   if (allocated(this%v)) deallocate(this%v)
   if (allocated(this%h)) deallocate(this%h)
   if (allocated(this%gamma)) deallocate(this%gamma)
   if (allocated(this%givens_c)) deallocate(this%givens_c)
   if (allocated(this%givens_s)) deallocate(this%givens_s)
   if (allocated(this%y)) deallocate(this%y)

   this%max_iter = 0
   this%abs_tol = 0.0_rp
   this%monitor = .false.
  end subroutine pnpn2_prs_gmres_free

  !> Solve the pressure system in the local \f$ lx2 \f$ split space.
  function pnpn2_prs_gmres_solve(this, Ax, M, x, f, n, coef, blst, &
      prs_dirichlet, glb_prs_points) result(ksp_results)
   class(pnpn2_prs_gmres_t), intent(inout) :: this
   class(ax_t), intent(in) :: Ax
   class(pc_t), intent(inout) :: M
   type(field_t), intent(inout) :: x
   integer, intent(in) :: n
   real(kind=rp), intent(in) :: f(n)
   type(coef_t), intent(inout) :: coef
   type(bc_list_t), intent(inout) :: blst
   logical, intent(in) :: prs_dirichlet
   integer(kind=i8), intent(in) :: glb_prs_points
   type(ksp_monitor_t) :: ksp_results
   character(len=LOG_SIZE) :: log_buf
   integer :: i, j, k
   integer :: iter, restart_len, m_used
   real(kind=rp) :: norm_fac, rnorm
   real(kind=xp) :: alpha, lr, temp
   logical :: converged

   converged = .false.
   iter = 0
   rnorm = 0.0_rp
   restart_len = min(this%lgmres, this%max_iter)

   if (restart_len .le. 0) then
     call neko_error('PnPn-2 pressure GMRES requires max_iter > 0.')
   end if

   norm_fac = 1.0_rp / sqrt(coef%volume)

   call copy(this%rhs, f, n)
   if (.not. prs_dirichlet) then
     call ortho(this%rhs, glb_prs_points, n)
   end if

   call rzero(x%x, n)

   do while ((.not. converged) .and. (iter .lt. this%max_iter))
     this%h = 0.0_xp
     this%gamma = 0.0_xp
     this%givens_c = 1.0_xp
     this%givens_s = 0.0_xp
     this%y = 0.0_xp
     m_used = 0

     call copy(this%r, this%rhs, n)
     if (iter .gt. 0) then
       call Ax%compute(this%w, x%x, coef, x%msh, x%Xh)
       call blst%apply(this%w, n)
       call add2s2(this%r, this%w, -1.0_rp, n)
       if (.not. prs_dirichlet) then
         call ortho(this%r, glb_prs_points, n)
       end if
     end if

     call copy(this%r_hat, this%r, n)
     call col2(this%r_hat, this%ml, n)
     this%gamma(1) = sqrt(real(glsc2(this%r_hat, this%r_hat, n), xp))

     if (iter .eq. 0) then
       ksp_results%res_start = real(this%gamma(1), rp) * norm_fac
     end if

     if (this%gamma(1) .eq. 0.0_xp) then
       converged = .true.
       exit
     end if

     call cmult2(this%v(1,1), this%r_hat, 1.0_rp / real(this%gamma(1), rp), n)

     do j = 1, restart_len
       m_used = j
       iter = iter + 1

       call copy(this%w, this%v(1,j), n)
       call col2(this%w, this%mu, n)

       call M%solve(this%z(1,j), this%w, n)
       if (.not. prs_dirichlet) then
         call ortho(this%z(1,j), glb_prs_points, n)
       end if

       call Ax%compute(this%w, this%z(1,j), coef, x%msh, x%Xh)
       call blst%apply(this%w, n)
       call copy(this%w_hat, this%w, n)
       call col2(this%w_hat, this%ml, n)

       do i = 1, j
         this%h(i,j) = real(glsc2(this%w_hat, this%v(1,i), n), xp)
         call add2s2(this%w_hat, this%v(1,i), -real(this%h(i,j), rp), n)
       end do

       alpha = sqrt(real(glsc2(this%w_hat, this%w_hat, n), xp))
       this%h(j + 1,j) = alpha

       do i = 1, j - 1
         temp = this%h(i,j)
         this%h(i,j) = this%givens_c(i) * temp + this%givens_s(i) * &
              this%h(i + 1,j)
         this%h(i + 1,j) = -this%givens_s(i) * temp + this%givens_c(i) * &
              this%h(i + 1,j)
       end do

       if (alpha .eq. 0.0_xp) then
         this%gamma(j + 1) = 0.0_xp
         rnorm = 0.0_rp
         converged = .true.
         exit
       end if

       lr = sqrt(this%h(j,j) * this%h(j,j) + this%h(j + 1,j) * this%h(j + 1,j))
       this%givens_c(j) = this%h(j,j) / lr
       this%givens_s(j) = this%h(j + 1,j) / lr
       this%h(j,j) = lr
       this%gamma(j + 1) = -this%givens_s(j) * this%gamma(j)
       this%gamma(j) = this%givens_c(j) * this%gamma(j)

       rnorm = abs(real(this%gamma(j + 1), rp)) * norm_fac
       if (this%monitor) then
         write(log_buf, '(A,I6,1X,ES13.6)') 'PnPn-2 pressure GMRES', iter, rnorm
         call neko_log%message(log_buf)
       end if

       if (rnorm .lt. this%abs_tol) then
         converged = .true.
         exit
       end if

       if (iter .ge. this%max_iter) exit

       call cmult2(this%v(1,j + 1), this%w_hat, 1.0_rp / real(alpha, rp), n)
     end do

     if (m_used .eq. 0) cycle

     do k = m_used, 1, -1
       temp = this%gamma(k)
       do i = m_used, k + 1, -1
         temp = temp - this%h(k,i) * this%y(i)
       end do
       this%y(k) = temp / this%h(k,k)
     end do

     do i = 1, m_used
       call add2s2(x%x, this%z(1,i), real(this%y(i), rp), n)
     end do

     if (.not. prs_dirichlet) then
       call ortho(x%x, glb_prs_points, n)
     end if
   end do

   if (.not. prs_dirichlet) then
     call ortho(x%x, glb_prs_points, n)
   end if

   ksp_results%iter = iter
   ksp_results%res_final = rnorm
   ksp_results%converged = converged .and. (rnorm .le. this%abs_tol)
  end function pnpn2_prs_gmres_solve

  !> Initialize lower-order GL pressure coefficients without any gather-scatter.
  subroutine pnpn2_pressure_coef_init(c, msh, Yh, dm)
    type(coef_t), intent(inout), target :: c
    type(mesh_t), target, intent(inout) :: msh
    type(space_t), target, intent(inout) :: Yh
    type(dofmap_t), target, intent(inout) :: dm
    integer :: e, i, lxy, lyz, lxyz, ntot

    call c%free()

    c%msh => msh
    c%Xh => Yh
    c%dof => dm

    allocate(c%G11(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%G22(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%G33(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%G12(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%G13(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%G23(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dxdr(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dydr(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dzdr(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dxds(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dyds(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dzds(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dxdt(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dydt(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dzdt(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%drdx(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%drdy(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%drdz(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dsdx(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dsdy(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dsdz(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dtdx(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dtdy(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%dtdz(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%jac(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%jacinv(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%B(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%Binv(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%h1(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%h2(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%mult(Yh%lx, Yh%ly, Yh%lz, msh%nelv))
    allocate(c%cyc_msk(0:0))

    c%Blag => c%B
    c%Blaglag => c%B

    lxy = Yh%lx * Yh%ly
    lyz = Yh%ly * Yh%lz
    lxyz = Yh%lxyz
    ntot = dm%size()

    do e = 1, msh%nelv
      call mxm(Yh%dx, Yh%lx, dm%x(1,1,1,e), Yh%lx, c%dxdr(1,1,1,e), lyz)
      call mxm(Yh%dx, Yh%lx, dm%y(1,1,1,e), Yh%lx, c%dydr(1,1,1,e), lyz)
      call mxm(Yh%dx, Yh%lx, dm%z(1,1,1,e), Yh%lx, c%dzdr(1,1,1,e), lyz)

      do i = 1, Yh%lz
        call mxm(dm%x(1,1,i,e), Yh%lx, Yh%dyt, Yh%ly, c%dxds(1,1,i,e), Yh%ly)
        call mxm(dm%y(1,1,i,e), Yh%lx, Yh%dyt, Yh%ly, c%dyds(1,1,i,e), Yh%ly)
        call mxm(dm%z(1,1,i,e), Yh%lx, Yh%dyt, Yh%ly, c%dzds(1,1,i,e), Yh%ly)
      end do

      if (msh%gdim .eq. 3) then
        call mxm(dm%x(1,1,1,e), lxy, Yh%dzt, Yh%lz, c%dxdt(1,1,1,e), Yh%lz)
        call mxm(dm%y(1,1,1,e), lxy, Yh%dzt, Yh%lz, c%dydt(1,1,1,e), Yh%lz)
        call mxm(dm%z(1,1,1,e), lxy, Yh%dzt, Yh%lz, c%dzdt(1,1,1,e), Yh%lz)
      else
        c%dxdt(:,:,:,e) = 0.0_rp
        c%dydt(:,:,:,e) = 0.0_rp
        c%dzdt(:,:,:,e) = 1.0_rp
      end if
    end do

    if (msh%gdim .eq. 2) then
      do i = 1, ntot
        c%jac(i,1,1,1) = c%dxdr(i,1,1,1) * c%dyds(i,1,1,1) - &
             c%dxds(i,1,1,1) * c%dydr(i,1,1,1)
        c%drdx(i,1,1,1) = c%dyds(i,1,1,1)
        c%drdy(i,1,1,1) = -c%dxds(i,1,1,1)
        c%dsdx(i,1,1,1) = -c%dydr(i,1,1,1)
        c%dsdy(i,1,1,1) = c%dxdr(i,1,1,1)
        c%drdz(i,1,1,1) = 0.0_rp
        c%dsdz(i,1,1,1) = 0.0_rp
        c%dtdx(i,1,1,1) = 0.0_rp
        c%dtdy(i,1,1,1) = 0.0_rp
        c%dtdz(i,1,1,1) = 1.0_rp
      end do
    else
      do i = 1, ntot
        c%jac(i,1,1,1) = &
             c%dxdr(i,1,1,1) * c%dyds(i,1,1,1) * c%dzdt(i,1,1,1) + &
             c%dxdt(i,1,1,1) * c%dydr(i,1,1,1) * c%dzds(i,1,1,1) + &
             c%dxds(i,1,1,1) * c%dydt(i,1,1,1) * c%dzdr(i,1,1,1) - &
             c%dxdr(i,1,1,1) * c%dydt(i,1,1,1) * c%dzds(i,1,1,1) - &
             c%dxds(i,1,1,1) * c%dydr(i,1,1,1) * c%dzdt(i,1,1,1) - &
             c%dxdt(i,1,1,1) * c%dyds(i,1,1,1) * c%dzdr(i,1,1,1)
        c%drdx(i,1,1,1) = c%dyds(i,1,1,1) * c%dzdt(i,1,1,1) - &
             c%dydt(i,1,1,1) * c%dzds(i,1,1,1)
        c%drdy(i,1,1,1) = c%dxdt(i,1,1,1) * c%dzds(i,1,1,1) - &
             c%dxds(i,1,1,1) * c%dzdt(i,1,1,1)
        c%drdz(i,1,1,1) = c%dxds(i,1,1,1) * c%dydt(i,1,1,1) - &
             c%dxdt(i,1,1,1) * c%dyds(i,1,1,1)
        c%dsdx(i,1,1,1) = c%dydt(i,1,1,1) * c%dzdr(i,1,1,1) - &
             c%dydr(i,1,1,1) * c%dzdt(i,1,1,1)
        c%dsdy(i,1,1,1) = c%dxdr(i,1,1,1) * c%dzdt(i,1,1,1) - &
             c%dxdt(i,1,1,1) * c%dzdr(i,1,1,1)
        c%dsdz(i,1,1,1) = c%dxdt(i,1,1,1) * c%dydr(i,1,1,1) - &
             c%dxdr(i,1,1,1) * c%dydt(i,1,1,1)
        c%dtdx(i,1,1,1) = c%dydr(i,1,1,1) * c%dzds(i,1,1,1) - &
             c%dyds(i,1,1,1) * c%dzdr(i,1,1,1)
        c%dtdy(i,1,1,1) = c%dxds(i,1,1,1) * c%dzdr(i,1,1,1) - &
             c%dxdr(i,1,1,1) * c%dzds(i,1,1,1)
        c%dtdz(i,1,1,1) = c%dxdr(i,1,1,1) * c%dyds(i,1,1,1) - &
             c%dxds(i,1,1,1) * c%dydr(i,1,1,1)
      end do
    end if

    do i = 1, ntot
      c%jacinv(i,1,1,1) = 1.0_rp / c%jac(i,1,1,1)
      c%drdx(i,1,1,1) = c%drdx(i,1,1,1) * c%jacinv(i,1,1,1)
      c%drdy(i,1,1,1) = c%drdy(i,1,1,1) * c%jacinv(i,1,1,1)
      c%drdz(i,1,1,1) = c%drdz(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dsdx(i,1,1,1) = c%dsdx(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dsdy(i,1,1,1) = c%dsdy(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dsdz(i,1,1,1) = c%dsdz(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dtdx(i,1,1,1) = c%dtdx(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dtdy(i,1,1,1) = c%dtdy(i,1,1,1) * c%jacinv(i,1,1,1)
      c%dtdz(i,1,1,1) = c%dtdz(i,1,1,1) * c%jacinv(i,1,1,1)
    end do

    do i = 1, ntot
      c%G11(i,1,1,1) = (c%drdx(i,1,1,1)**2 + c%drdy(i,1,1,1)**2 + &
           c%drdz(i,1,1,1)**2) * c%jacinv(i,1,1,1)
      c%G22(i,1,1,1) = (c%dsdx(i,1,1,1)**2 + c%dsdy(i,1,1,1)**2 + &
           c%dsdz(i,1,1,1)**2) * c%jacinv(i,1,1,1)
      c%G33(i,1,1,1) = (c%dtdx(i,1,1,1)**2 + c%dtdy(i,1,1,1)**2 + &
           c%dtdz(i,1,1,1)**2) * c%jacinv(i,1,1,1)
      c%G12(i,1,1,1) = (c%drdx(i,1,1,1) * c%dsdx(i,1,1,1) + &
           c%drdy(i,1,1,1) * c%dsdy(i,1,1,1) + &
           c%drdz(i,1,1,1) * c%dsdz(i,1,1,1)) * c%jacinv(i,1,1,1)
      c%G13(i,1,1,1) = (c%drdx(i,1,1,1) * c%dtdx(i,1,1,1) + &
           c%drdy(i,1,1,1) * c%dtdy(i,1,1,1) + &
           c%drdz(i,1,1,1) * c%dtdz(i,1,1,1)) * c%jacinv(i,1,1,1)
      c%G23(i,1,1,1) = (c%dsdx(i,1,1,1) * c%dtdx(i,1,1,1) + &
           c%dsdy(i,1,1,1) * c%dtdy(i,1,1,1) + &
           c%dsdz(i,1,1,1) * c%dtdz(i,1,1,1)) * c%jacinv(i,1,1,1)
    end do

    do e = 1, msh%nelv
      do i = 1, lxyz
        c%B(i,1,1,e) = c%jac(i,1,1,e) * Yh%w3(i,1,1)
        c%Binv(i,1,1,e) = 1.0_rp / c%B(i,1,1,e)
      end do
    end do

    c%h1 = 1.0_rp
    c%h2 = 1.0_rp
    c%mult = 1.0_rp
    c%ifh2 = .false.
    c%volume = glsum(c%B, ntot)
    c%cyc_msk(0) = 1
  end subroutine pnpn2_pressure_coef_init

  !> Initialize the milestone-1 PnPn-2 fluid scheme.
  subroutine fluid_pnpn2_init(this, msh, lx, params, user, chkp)
    class(fluid_pnpn2_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    character(len=17), parameter :: scheme = 'pnpn2 (Pn/Pn-2)'
    character(len=LOG_SIZE) :: log_buf
    integer :: lx2
    integer :: time_order
    integer :: solver_maxiter
    real(kind=rp) :: abs_tol
    logical :: monitor, advection
    character(len=:), allocatable :: solver_type, precon_type
    type(json_file) :: numerics_params, precon_params

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2 milestone-1 is currently CPU-only.')
    end if

    if (lx .le. 2) then
      call neko_error('pnpn2 requires polynomial order at least 3.')
    end if

    call this%init_base(msh, lx, params, scheme, user, .true.)

    lx2 = lx - 2
    if (msh%gdim .eq. 2) then
      call this%Yh%init(GL, lx2, lx2)
    else
      call this%Yh%init(GL, lx2, lx2, lx2)
    end if

    if (this%pr_projection_dim .gt. 0) then
      call neko_error('pnpn2 milestone-1 does not support pressure projection.')
    end if

    call this%dm_Yh%init(msh, this%Yh)
    call this%gs_prs%init_noop(this%dm_Yh)
    call pnpn2_pressure_coef_init(this%c_Yh, msh, this%Yh, this%dm_Yh)
    call this%prs_interp%init(this%Xh, this%Yh)

    call neko_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_registry%get_field('p')

    call json_get_or_lookup(params, 'case.numerics.time_order', time_order)
    allocate(this%ext_bdf)
    call this%ext_bdf%init(time_order)

    call ax_helm_factory(this%Ax_vel, full_formulation = .false.)
    allocate(pnpn2_prs_ax_t::this%Ax_prs)
    call this%mixed_ops%init(this%Xh, this%Yh, this%dm_Xh, this%dm_Yh, &
         this%c_Xh, this%c_Yh)
    call pnpn2_prs_ax_init(this%mixed_ops)

    call rhs_maker_bdf_fctry(this%makebdf)

    call this%p_Yh%init(this%dm_Yh, 'p_Yh')
    call this%p_res%init(this%dm_Yh, 'p_res')
    call this%dp%init(this%dm_Yh, 'dp')
    call this%p_ext%init(this%dm_Yh, 'p_ext')
    call this%plag%init(this%p_Yh, 3)
    call this%u_res%init(this%dm_Xh, 'u_res')
    call this%v_res%init(this%dm_Xh, 'v_res')
    call this%w_res%init(this%dm_Xh, 'w_res')
    call this%du%init(this%dm_Xh, 'du')
    call this%dv%init(this%dm_Xh, 'dv')
    call this%dw%init(this%dm_Xh, 'dw')

    call this%setup_bcs(user, params)

    call neko_log%section('Pressure solver')
    call json_get_or_lookup_or_default(params, &
         'case.fluid.pressure_solver.max_iterations', solver_maxiter, 800)
    call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
    call json_get(params, 'case.fluid.pressure_solver.preconditioner.type', &
         precon_type)
    call json_get(params, 'case.fluid.pressure_solver.preconditioner', &
         precon_params)
    call json_get_or_lookup(params, &
         'case.fluid.pressure_solver.absolute_tolerance', abs_tol)
    call json_get_or_default(params, 'case.fluid.pressure_solver.monitor', &
         monitor, .false.)
    call neko_log%message('Type       : (' // trim(solver_type) // ', ' // &
         trim(precon_type) // ')')
    write(log_buf, '(A,ES13.6)') 'Abs tol    :', abs_tol
    call neko_log%message(log_buf)
    if (trim(solver_type) .ne. 'gmres') then
      call neko_error('pnpn2 milestone-1 pressure solve currently requires gmres.')
    end if
    call this%solver_factory(this%ksp_prs, this%dm_Yh%size(), solver_type, &
         solver_maxiter, abs_tol, monitor)
    if (trim(precon_type) .eq. 'hsmg') then
      allocate(hsmg_pnpn2_t :: this%pc_prs)
      select type (pc_prs => this%pc_prs)
      type is (hsmg_pnpn2_t)
         call pc_prs%init(this%c_Xh, this%Yh, this%dm_Yh, this%dm_Xh, &
             this%gs_Xh, this%bcs_prs, precon_params)
      end select
    else
      allocate(pnpn2_prs_precon_t :: this%pc_prs)
      select type (pc_prs => this%pc_prs)
      type is (pnpn2_prs_precon_t)
         call pc_prs%init(precon_type, precon_params, this%mixed_ops, &
             this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_prs, &
             this%prs_dirichlet)
      end select
    end if
    call this%ksp_prs%set_pc(this%pc_prs)
    call this%prs_gmres%init(this%dm_Yh%size(), this%ksp_prs%max_iter, &
         this%ksp_prs%abs_tol, this%ksp_prs%monitor, this%c_Yh%B, &
         this%c_Yh%Binv)
    call neko_log%end_section()

    call json_get_or_default(params, 'case.numerics.oifs', pnpn2_oifs, .false.)
    call json_get_or_default(params, 'case.fluid.advection', advection, .true.)
    call json_get(params, 'case.numerics', numerics_params)
    call rhs_maker_ext_fctry(pnpn2_makeabf)
    call pnpn2_abx1%init(this%dm_Xh, 'abx1')
    call pnpn2_aby1%init(this%dm_Xh, 'aby1')
    call pnpn2_abz1%init(this%dm_Xh, 'abz1')
    call pnpn2_abx2%init(this%dm_Xh, 'abx2')
    call pnpn2_aby2%init(this%dm_Xh, 'aby2')
    call pnpn2_abz2%init(this%dm_Xh, 'abz2')
    if (pnpn2_oifs) then
      call rhs_maker_oifs_fctry(pnpn2_makeoifs)
      call pnpn2_advx%init(this%dm_Xh, 'advx')
      call pnpn2_advy%init(this%dm_Xh, 'advy')
      call pnpn2_advz%init(this%dm_Xh, 'advz')
    end if
    call advection_factory(pnpn2_adv, numerics_params, this%c_Xh, &
         this%ulag, this%vlag, this%wlag, chkp%dtlag, chkp%tlag, &
         this%ext_bdf, .not. advection)

    this%glb_prs_points = int(this%msh%glb_nelv, i8) * int(this%Yh%lxyz, i8)

    this%chkp => chkp
    call this%chkp%add_fluid(this%u, this%v, this%w, this%p)
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)

    call neko_log%message('Milestone-1 pnpn2: source terms are disabled.')
    call neko_log%end_section()
  end subroutine fluid_pnpn2_init

  !> Free the milestone-1 PnPn-2 scheme.
  subroutine fluid_pnpn2_free(this)
    class(fluid_pnpn2_t), intent(inout) :: this

    call pnpn2_prs_ax_clear()
    call this%prs_gmres%free()

    call this%scheme_free()

    if (allocated(this%ext_bdf)) then
      call this%ext_bdf%free()
      deallocate(this%ext_bdf)
    end if

    if (allocated(this%Ax_vel)) then
      deallocate(this%Ax_vel)
    end if
    if (allocated(this%Ax_prs)) then
      deallocate(this%Ax_prs)
    end if
    if (allocated(pnpn2_adv)) then
      call pnpn2_adv%free()
      deallocate(pnpn2_adv)
    end if
    if (allocated(pnpn2_makeabf)) then
      deallocate(pnpn2_makeabf)
    end if
    if (allocated(this%makebdf)) then
      deallocate(this%makebdf)
    end if
    if (allocated(pnpn2_makeoifs)) then
      deallocate(pnpn2_makeoifs)
    end if

    call this%plag%free()
    call this%mixed_ops%free()

    call this%p_Yh%free()
    call this%p_res%free()
    call this%dp%free()
    call this%p_ext%free()
    call this%u_res%free()
    call this%v_res%free()
    call this%w_res%free()
    call pnpn2_abx1%free()
    call pnpn2_aby1%free()
    call pnpn2_abz1%free()
    call pnpn2_abx2%free()
    call pnpn2_aby2%free()
    call pnpn2_abz2%free()
    call pnpn2_advx%free()
    call pnpn2_advy%free()
    call pnpn2_advz%free()
    call this%du%free()
    call this%dv%free()
    call this%dw%free()

    call this%bc_vel_res%free()
    call this%bc_du%free()
    call this%bc_dv%free()
    call this%bc_dw%free()

    call this%bclst_dp%free()
    call this%bclst_vel_res%free()
    call this%bclst_du%free()
    call this%bclst_dv%free()
    call this%bclst_dw%free()

    call this%prs_interp%free()
    call this%c_Yh%free()
    call this%gs_prs%free()
    call this%dm_Yh%free()
    call this%Yh%free()

    this%glb_prs_points = 0_i8
    this%prs_dirichlet = .false.
  end subroutine fluid_pnpn2_free

  !> Advance the bare-bones PnPn-2 milestone-1 scheme by one time step.
  subroutine fluid_pnpn2_step(this, time, dt_controller)
    class(fluid_pnpn2_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    type(field_t), pointer :: gx, gy, gz
    integer :: scratch_ids(3)
    integer :: n_x, n_y
    integer :: i
    real(kind=rp) :: rho_val, mu_val, a0, dt
    type(ksp_monitor_t) :: ksp_results(4)

    if (this%freeze) return

    n_x = this%dm_Xh%size()
    n_y = this%dm_Yh%size()
    rho_val = this%rho%x(1,1,1,1)
    mu_val = this%mu_tot%x(1,1,1,1)
    a0 = this%ext_bdf%diffusion_coeffs%x(1)
    dt = time%dt

    if (a0 .eq. 0.0_rp) then
      call neko_error('pnpn2 received zero leading BDF coefficient.')
    end if

    call this%update_material_properties(time)

    rho_val = this%rho%x(1,1,1,1)
    mu_val = this%mu_tot%x(1,1,1,1)

    call this%plag%update()

    call this%bc_apply_vel(time, strong = .true.)

    if (this%ext_bdf%ndiff .le. 2) then
      do concurrent (i = 1:n_y)
        this%p_ext%x(i,1,1,1) = this%plag%lf(1)%x(i,1,1,1)
      end do
    else
      do concurrent (i = 1:n_y)
        this%p_ext%x(i,1,1,1) = this%ext_bdf%advection_coeffs%x(1) * &
             this%plag%lf(1)%x(i,1,1,1) + &
             this%ext_bdf%advection_coeffs%x(2) * this%plag%lf(2)%x(i,1,1,1) + &
             this%ext_bdf%advection_coeffs%x(3) * this%plag%lf(3)%x(i,1,1,1)
      end do
    end if

    do concurrent (i = 1:n_x)
      this%c_Xh%h1(i,1,1,1) = mu_val
      this%c_Xh%h2(i,1,1,1) = rho_val * a0 / dt
    end do
    this%c_Xh%ifh2 = .true.

    do concurrent (i = 1:n_y)
      this%c_Yh%h1(i,1,1,1) = 1.0_rp / rho_val
    end do
    this%c_Yh%ifh2 = .false.

    call this%source_term%compute(time)
    call this%bcs_vel%apply_vector(this%f_x%x, this%f_y%x, this%f_z%x, n_x, &
         time, strong = .false.)

    if (pnpn2_oifs) then
      call pnpn2_adv%compute(this%u, this%v, this%w, pnpn2_advx, pnpn2_advy, &
           pnpn2_advz, this%Xh, this%c_Xh, n_x, dt)

      call pnpn2_makeabf%compute_fluid(pnpn2_abx1, pnpn2_aby1, pnpn2_abz1, &
           pnpn2_abx2, pnpn2_aby2, pnpn2_abz2, this%f_x%x, this%f_y%x, &
           this%f_z%x, rho_val, this%ext_bdf%advection_coeffs%x, n_x)

      call pnpn2_makeoifs%compute_fluid(pnpn2_advx%x, pnpn2_advy%x, &
           pnpn2_advz%x, this%f_x%x, this%f_y%x, this%f_z%x, rho_val, dt, n_x)
    else
      call pnpn2_adv%compute(this%u, this%v, this%w, this%f_x, this%f_y, &
           this%f_z, this%Xh, this%c_Xh, n_x)

      call pnpn2_makeabf%compute_fluid(pnpn2_abx1, pnpn2_aby1, pnpn2_abz1, &
           pnpn2_abx2, pnpn2_aby2, pnpn2_abz2, this%f_x%x, this%f_y%x, &
           this%f_z%x, rho_val, this%ext_bdf%advection_coeffs%x, n_x)

      call this%makebdf%compute_fluid(this%ulag, this%vlag, this%wlag, &
           this%f_x%x, this%f_y%x, this%f_z%x, this%u, this%v, this%w, &
           this%c_Xh%B, rho_val, dt, this%ext_bdf%diffusion_coeffs%x, &
           this%ext_bdf%ndiff, n_x, this%c_Xh%Blag, this%c_Xh%Blaglag)
    end if

    call this%ulag%update()
    call this%vlag%update()
    call this%wlag%update()

    call neko_scratch_registry%request_field(gx, scratch_ids(1), .false.)
    call neko_scratch_registry%request_field(gy, scratch_ids(2), .false.)
    call neko_scratch_registry%request_field(gz, scratch_ids(3), .false.)

    call this%Ax_vel%compute(this%u_res%x, this%u%x, this%c_Xh, this%msh, &
         this%Xh)
    call this%Ax_vel%compute(this%v_res%x, this%v%x, this%c_Xh, this%msh, &
         this%Xh)
    call this%Ax_vel%compute(this%w_res%x, this%w%x, this%c_Xh, this%msh, &
         this%Xh)
    call this%mixed_ops%opgradt(gx%x, gy%x, gz%x, this%p_ext%x)

    do concurrent (i = 1:n_x)
      this%u_res%x(i,1,1,1) = this%f_x%x(i,1,1,1) - this%u_res%x(i,1,1,1) + &
           gx%x(i,1,1,1)
      this%v_res%x(i,1,1,1) = this%f_y%x(i,1,1,1) - this%v_res%x(i,1,1,1) + &
           gy%x(i,1,1,1)
      this%w_res%x(i,1,1,1) = this%f_z%x(i,1,1,1) - this%w_res%x(i,1,1,1) + &
           gz%x(i,1,1,1)
    end do

    call rotate_cyc(this%u_res%x, this%v_res%x, this%w_res%x, 1, this%c_Xh)
    call this%gs_Xh%op(this%u_res, GS_OP_ADD)
    call this%gs_Xh%op(this%v_res, GS_OP_ADD)
    call this%gs_Xh%op(this%w_res, GS_OP_ADD)
    call rotate_cyc(this%u_res%x, this%v_res%x, this%w_res%x, 0, this%c_Xh)

    call this%bclst_vel_res%apply_vector(this%u_res%x, this%v_res%x, &
         this%w_res%x, n_x, time)

    call this%pc_vel%update()
    ksp_results(2:4) = this%ksp_vel%solve_coupled(this%Ax_vel, this%du, &
         this%dv, this%dw, this%u_res%x, this%v_res%x, this%w_res%x, n_x, &
         this%c_Xh, this%bclst_du, this%bclst_dv, this%bclst_dw, &
         this%gs_Xh, this%ksp_vel%max_iter)
    ksp_results(2)%name = 'X-Velocity'
    ksp_results(3)%name = 'Y-Velocity'
    ksp_results(4)%name = 'Z-Velocity'

    call opadd2cm(this%u%x, this%v%x, this%w%x, this%du%x, this%dv%x, &
         this%dw%x, 1.0_rp, n_x, this%msh%gdim)

    call this%mixed_ops%opdiv(this%p_res%x, this%u%x, this%v%x, this%w%x)
    call field_cmult(this%p_res, -1.0_rp, n_y)
    call this%bclst_dp%apply_scalar(this%p_res%x, n_y, time)
    if (.not. this%prs_dirichlet) then
     call ortho(this%p_res%x, this%glb_prs_points, n_y)
    end if

    call this%pc_prs%update()
    ksp_results(1) = this%prs_gmres%solve(this%Ax_prs, this%pc_prs, this%dp, &
         this%p_res%x, n_y, this%c_Yh, this%bclst_dp, this%prs_dirichlet, &
         this%glb_prs_points)
    ksp_results(1)%name = 'Pressure'

    if (.not. this%prs_dirichlet) then
      call ortho(this%dp%x, this%glb_prs_points, n_y)
    end if
    call field_cmult(this%dp, a0 / dt, n_y)

    do concurrent (i = 1:n_y)
      this%p_Yh%x(i,1,1,1) = this%p_ext%x(i,1,1,1) + this%dp%x(i,1,1,1)
    end do
    if (.not. this%prs_dirichlet) then
      call ortho(this%p_Yh%x, this%glb_prs_points, n_y)
    end if

    call this%mixed_ops%opgradt(gx%x, gy%x, gz%x, this%dp%x)
    call this%gs_Xh%op(gx, GS_OP_ADD)
    call this%gs_Xh%op(gy, GS_OP_ADD)
    call this%gs_Xh%op(gz, GS_OP_ADD)
    call col2(gx%x, this%c_Xh%Binv, n_x)
    call col2(gy%x, this%c_Xh%Binv, n_x)
    call col2(gz%x, this%c_Xh%Binv, n_x)
    call opadd2cm(this%u%x, this%v%x, this%w%x, gx%x, gy%x, gz%x, dt / a0, &
         n_x, this%msh%gdim)

    call this%bc_apply_vel(time, strong = .true.)
    call this%sync_p_public()
    call this%bc_apply_prs(time)

    call neko_scratch_registry%relinquish_field(scratch_ids)

    call fluid_step_info(time, ksp_results, .false., this%strict_convergence, &
         this%allow_stabilization, 1)
  end subroutine fluid_pnpn2_step

  !> Restart the milestone-1 PnPn-2 scheme.
  subroutine fluid_pnpn2_restart(this, chkp)
    class(fluid_pnpn2_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2 restart is currently CPU-only.')
    end if

    call this%sync_p_from_public()
    call this%plag%set(this%p_Yh)
    call pnpn2_prs_ax_init(this%mixed_ops)
  end subroutine fluid_pnpn2_restart

  !> Set up the supported boundary-condition configuration for milestone 1.
  subroutine fluid_pnpn2_setup_bcs(this, user, params)
    class(fluid_pnpn2_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: global_zone_size, i, ierr, j, n_bcs, zone_size
    logical :: found
    logical :: has_pressure_bc
    logical, allocatable :: marked_zones(:)
    integer, allocatable :: zone_indices(:)
    class(bc_t), pointer :: bc_i
    type(json_core) :: core
    type(json_file) :: bc_subdict
    type(json_value), pointer :: bc_object

    call this%bclst_vel_res%init()
    call this%bclst_du%init()
    call this%bclst_dv%init()
    call this%bclst_dw%init()
    call this%bclst_dp%init()

    call this%bc_vel_res%init_from_components(this%c_Xh)
    call this%bc_du%init_from_components(this%c_Xh)
    call this%bc_dv%init_from_components(this%c_Xh)
    call this%bc_dw%init_from_components(this%c_Xh)
    has_pressure_bc = .false.

    if (params%valid_path('case.fluid.boundary_conditions')) then
      call params%info('case.fluid.boundary_conditions', n_children = n_bcs)
      call params%get_core(core)
      call params%get('case.fluid.boundary_conditions', bc_object, found)

      call this%bcs_vel%init(n_bcs)
      call this%bcs_prs%init(n_bcs)

      allocate(marked_zones(size(this%msh%labeled_zones)))
      marked_zones = .false.

      do i = 1, n_bcs
        call json_extract_item(core, bc_object, i, bc_subdict)
        call json_get_or_lookup(bc_subdict, 'zone_indices', zone_indices)

        do j = 1, size(zone_indices)
          zone_size = this%msh%labeled_zones(zone_indices(j))%size
          call MPI_Allreduce(zone_size, global_zone_size, 1, MPI_INTEGER, &
               MPI_MAX, NEKO_COMM, ierr)

          if (global_zone_size .eq. 0) then
            write(error_unit, '(A, A, I0, A, A, I0, A)') "*** ERROR ***: ", &
                 "Zone index ", zone_indices(j), &
                 " is invalid as this zone has 0 size, meaning it ", &
                 "is not in the mesh. Check fluid boundary condition ", i, "."
            error stop
          end if

          if (marked_zones(zone_indices(j))) then
            write(error_unit, '(A, A, I0, A, A, A, A)') "*** ERROR ***: ", &
                 "Zone with index ", zone_indices(j), &
                 " has already been assigned a boundary condition. ", &
                 "Please check your boundary_conditions entry for the ", &
                 "fluid and make sure that each zone index appears only ", &
                 "in a single boundary condition."
            error stop
          else
            marked_zones(zone_indices(j)) = .true.
          end if
        end do

        bc_i => null()
        call pnpn2_velocity_bc_factory(bc_i, bc_subdict, this%c_Xh)
        if (associated(bc_i)) then
          if (bc_i%strong) then
            call this%bc_vel_res%mark_facets(bc_i%marked_facet)
            call this%bc_du%mark_facets(bc_i%marked_facet)
            call this%bc_dv%mark_facets(bc_i%marked_facet)
            call this%bc_dw%mark_facets(bc_i%marked_facet)
          end if
          call this%bcs_vel%append(bc_i)
        end if

        bc_i => null()
        call pnpn2_pressure_bc_factory(bc_i, bc_subdict, this%c_Xh)
        if (associated(bc_i)) then
          has_pressure_bc = .true.
          call this%bcs_prs%append(bc_i)
        end if

        if (allocated(zone_indices)) then
          deallocate(zone_indices)
        end if
      end do

      do i = 1, size(this%msh%labeled_zones)
        if ((this%msh%labeled_zones(i)%size .gt. 0) .and. &
             (.not. marked_zones(i))) then
          write(error_unit, '(A, A, I0)') "*** ERROR ***: ", &
               "No fluid boundary condition assigned to zone ", i
          error stop
        end if
      end do
    else
      do i = 1, size(this%msh%labeled_zones)
        if (this%msh%labeled_zones(i)%size .gt. 0) then
          call neko_error('No boundary_conditions entry in the case file!')
        end if
      end do

      call this%bcs_vel%init()
      call this%bcs_prs%init()
    end if

    call this%bc_vel_res%finalize()
    call this%bc_du%finalize()
    call this%bc_dv%finalize()
    call this%bc_dw%finalize()
    call this%bclst_vel_res%append(this%bc_vel_res)
    call this%bclst_du%append(this%bc_du)
    call this%bclst_dv%append(this%bc_dv)
    call this%bclst_dw%append(this%bc_dw)

    this%prs_dirichlet = has_pressure_bc
    call MPI_Allreduce(MPI_IN_PLACE, this%prs_dirichlet, 1, MPI_LOGICAL, &
         MPI_LOR, NEKO_COMM)

    if (allocated(marked_zones)) then
      deallocate(marked_zones)
    end if
    if (allocated(zone_indices)) then
      deallocate(zone_indices)
    end if
  end subroutine fluid_pnpn2_setup_bcs

  !> Refresh the public \f$ X_h \f$ pressure view from the authoritative \f$ Y_h \f$ state.
  subroutine fluid_pnpn2_sync_p_public(this)
    class(fluid_pnpn2_t), intent(inout) :: this

    call this%prs_interp%map(this%p%x, this%p_Yh%x, this%msh%nelv, this%Xh)
  end subroutine fluid_pnpn2_sync_p_public

  !> Refresh the authoritative \f$ Y_h \f$ pressure from the public \f$ X_h \f$ field.
  subroutine fluid_pnpn2_sync_p_from_public(this)
    class(fluid_pnpn2_t), intent(inout) :: this

    call this%prs_interp%map(this%p_Yh%x, this%p%x, this%msh%nelv, this%Yh)
  end subroutine fluid_pnpn2_sync_p_from_public

end module fluid_pnpn2
