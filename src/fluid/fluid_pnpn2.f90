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
  use ax_product, only : ax_t, ax_helm_factory
  use bc_list, only : bc_list_t
  use checkpoint, only : chkp_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_math, only : field_add2, field_cmult, field_rzero
  use field_series, only : field_series_t
  use fluid_aux, only : fluid_step_info
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use interpolation, only : interpolator_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use krylov, only : ksp_monitor_t
  use logger, only : neko_log, LOG_SIZE
  use math, only : add2, col2
  use mathops, only : opadd2cm
  use mesh, only : mesh_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : i8, rp
  use operators, only : cdtp, opgrad, ortho, rotate_cyc
  use pnpn2_prs_ax, only : pnpn2_prs_ax_t, pnpn2_prs_ax_clear, &
       pnpn2_prs_ax_init
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use registry, only : neko_registry
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_bdf_fctry
  use scratch_registry, only : neko_scratch_registry
  use space, only : GLL, space_t
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t
  use utils, only : neko_error
  implicit none
  private

  !> Bare-bones CPU PnPn-2 pressure-correction scheme.
  type, public, extends(fluid_scheme_incompressible_t) :: fluid_pnpn2_t
     !> Pressure function space \f$ Y_h = P_{N-2} \f$.
     type(space_t) :: Yh
     !> Pressure dofmap on \f$ Y_h \f$.
     type(dofmap_t) :: dm_Yh
     !> Pressure gather-scatter on \f$ Y_h \f$.
     type(gs_t) :: gs_Yh
     !> Coefficients on \f$ Y_h \f$.
     type(coef_t) :: c_Yh
     !> Interpolator between \f$ X_h \f$ and \f$ Y_h \f$.
     type(interpolator_t) :: prs_interp
     !> Non-unique global number of lower-order pressure points.
     integer(kind=i8) :: glb_prs_points = 0_i8
     !> Pressure right-hand side.
     type(field_t) :: p_res
     !> Pressure increment.
     type(field_t) :: dp
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
     !> Empty pressure-increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_dp
     !> Empty x-velocity increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_du
     !> Empty y-velocity increment boundary list for the linear solver.
     type(bc_list_t) :: bclst_dv
     !> Empty z-velocity increment boundary list for the linear solver.
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
  end type fluid_pnpn2_t

contains

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
    logical :: monitor
    character(len=:), allocatable :: solver_type, precon_type
    type(json_file) :: precon_params

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2 milestone-1 is currently CPU-only.')
    end if

    if (lx .le. 2) then
      call neko_error('pnpn2 requires polynomial order at least 3.')
    end if

    call this%init_base(msh, lx, params, scheme, user, .true.)

    if (params%valid_path('case.fluid.source_terms')) then
      call neko_error('pnpn2 milestone-1 does not support source terms.')
    end if

    lx2 = lx - 2
    if (msh%gdim .eq. 2) then
      call this%Yh%init(GLL, lx2, lx2)
    else
      call this%Yh%init(GLL, lx2, lx2, lx2)
    end if

    call this%dm_Yh%init(msh, this%Yh)
    call this%gs_Yh%init(this%dm_Yh)
    call this%c_Yh%init(this%gs_Yh)
    call this%prs_interp%init(this%Xh, this%Yh)

    call neko_registry%add_field(this%dm_Yh, 'p')
    this%p => neko_registry%get_field('p')

    call json_get_or_lookup(params, 'case.numerics.time_order', time_order)
    allocate(this%ext_bdf)
    call this%ext_bdf%init(time_order)

    call ax_helm_factory(this%Ax_vel, full_formulation = .false.)
    allocate(pnpn2_prs_ax_t::this%Ax_prs)
    call pnpn2_prs_ax_init(this%Xh, this%Yh, this%c_Xh, this%prs_interp)

    call rhs_maker_bdf_fctry(this%makebdf)

    call this%p_res%init(this%dm_Yh, 'p_res')
    call this%dp%init(this%dm_Yh, 'dp')
    call this%u_res%init(this%dm_Xh, 'u_res')
    call this%v_res%init(this%dm_Xh, 'v_res')
    call this%w_res%init(this%dm_Xh, 'w_res')
    call this%du%init(this%dm_Xh, 'du')
    call this%dv%init(this%dm_Xh, 'dv')
    call this%dw%init(this%dm_Xh, 'dw')

    call this%setup_bcs(user, params)

    call this%proj_prs%init(this%dm_Yh%size(), this%pr_projection_dim, &
         this%pr_projection_activ_step, &
         this%pr_projection_reorthogonalize_basis)
    call this%proj_vel%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)

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
    call this%solver_factory(this%ksp_prs, this%dm_Yh%size(), solver_type, &
         solver_maxiter, abs_tol, monitor)
    call this%precon_factory_(this%pc_prs, this%ksp_prs, this%c_Yh, &
         this%dm_Yh, this%gs_Yh, this%bcs_prs, precon_type, precon_params)
    call neko_log%end_section()

    this%glb_prs_points = int(this%msh%glb_nelv, i8) * int(this%Yh%lxyz, i8)

    this%chkp => chkp
    call this%chkp%add_fluid(this%u, this%v, this%w, this%p)
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)

    call neko_log%message('Milestone-1 pnpn2: advection, source terms, and ' // &
         'boundary-condition terms are disabled.')
    call neko_log%end_section()
  end subroutine fluid_pnpn2_init

  !> Free the milestone-1 PnPn-2 scheme.
  subroutine fluid_pnpn2_free(this)
    class(fluid_pnpn2_t), intent(inout) :: this

    call pnpn2_prs_ax_clear()

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
    if (allocated(this%makebdf)) then
      deallocate(this%makebdf)
    end if

    call this%proj_prs%free()
    call this%proj_vel%free()

    call this%p_res%free()
    call this%dp%free()
    call this%u_res%free()
    call this%v_res%free()
    call this%w_res%free()
    call this%du%free()
    call this%dv%free()
    call this%dw%free()

    call this%bclst_dp%free()
    call this%bclst_du%free()
    call this%bclst_dv%free()
    call this%bclst_dw%free()

    call this%prs_interp%free()
    call this%c_Yh%free()
    call this%gs_Yh%free()
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
    type(field_t), pointer :: p_Xh, gx, gy, gz, div_Xh
    integer :: scratch_ids(5)
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

    call this%ulag%update()
    call this%vlag%update()
    call this%wlag%update()

    do concurrent (i = 1:n_x)
      this%c_Xh%h1(i,1,1,1) = mu_val
      this%c_Xh%h2(i,1,1,1) = rho_val * a0 / dt
    end do
    this%c_Xh%ifh2 = .true.

    do concurrent (i = 1:n_y)
      this%c_Yh%h1(i,1,1,1) = 1.0_rp / rho_val
    end do
    this%c_Yh%ifh2 = .false.

    call field_rzero(this%f_x)
    call field_rzero(this%f_y)
    call field_rzero(this%f_z)

    call this%makebdf%compute_fluid(this%ulag, this%vlag, this%wlag, &
         this%f_x%x, this%f_y%x, this%f_z%x, this%u, this%v, this%w, &
         this%c_Xh%B, rho_val, dt, this%ext_bdf%diffusion_coeffs%x, &
         this%ext_bdf%ndiff, n_x, this%c_Xh%Blag, this%c_Xh%Blaglag)

    call neko_scratch_registry%request_field(p_Xh, scratch_ids(1), .false.)
    call neko_scratch_registry%request_field(gx, scratch_ids(2), .false.)
    call neko_scratch_registry%request_field(gy, scratch_ids(3), .false.)
    call neko_scratch_registry%request_field(gz, scratch_ids(4), .false.)
    call neko_scratch_registry%request_field(div_Xh, scratch_ids(5), .false.)

    call this%prs_interp%map(p_Xh%x, this%p%x, this%msh%nelv, this%Xh)

    call this%Ax_vel%compute(this%u_res%x, this%u%x, this%c_Xh, this%msh, &
         this%Xh)
    call this%Ax_vel%compute(this%v_res%x, this%v%x, this%c_Xh, this%msh, &
         this%Xh)
    call this%Ax_vel%compute(this%w_res%x, this%w%x, this%c_Xh, this%msh, &
         this%Xh)
    call opgrad(gx%x, gy%x, gz%x, p_Xh%x, this%c_Xh)

    do concurrent (i = 1:n_x)
      this%u_res%x(i,1,1,1) = this%f_x%x(i,1,1,1) - this%u_res%x(i,1,1,1) - &
           gx%x(i,1,1,1)
      this%v_res%x(i,1,1,1) = this%f_y%x(i,1,1,1) - this%v_res%x(i,1,1,1) - &
           gy%x(i,1,1,1)
      this%w_res%x(i,1,1,1) = this%f_z%x(i,1,1,1) - this%w_res%x(i,1,1,1) - &
           gz%x(i,1,1,1)
    end do

    call rotate_cyc(this%u_res%x, this%v_res%x, this%w_res%x, 1, this%c_Xh)
    call this%gs_Xh%op(this%u_res, GS_OP_ADD)
    call this%gs_Xh%op(this%v_res, GS_OP_ADD)
    call this%gs_Xh%op(this%w_res, GS_OP_ADD)
    call rotate_cyc(this%u_res%x, this%v_res%x, this%w_res%x, 0, this%c_Xh)

    call this%proj_vel%pre_solving(this%u_res%x, this%v_res%x, this%w_res%x, &
         time%tstep, this%c_Xh, n_x, dt_controller, 'Velocity')

    call this%pc_vel%update()
    ksp_results(2:4) = this%ksp_vel%solve_coupled(this%Ax_vel, this%du, &
         this%dv, this%dw, this%u_res%x, this%v_res%x, this%w_res%x, n_x, &
         this%c_Xh, this%bclst_du, this%bclst_dv, this%bclst_dw, &
         this%gs_Xh, this%ksp_vel%max_iter)
    ksp_results(2)%name = 'X-Velocity'
    ksp_results(3)%name = 'Y-Velocity'
    ksp_results(4)%name = 'Z-Velocity'

    call this%proj_vel%post_solving(this%du%x, this%dv%x, this%dw%x, &
         this%Ax_vel, this%c_Xh, this%bclst_du, this%bclst_dv, this%bclst_dw, &
         this%gs_Xh, n_x, time%tstep, dt_controller)

    call opadd2cm(this%u%x, this%v%x, this%w%x, this%du%x, this%dv%x, &
         this%dw%x, 1.0_rp, n_x, this%msh%gdim)

    call cdtp(div_Xh%x, this%u%x, this%c_Xh%drdx, this%c_Xh%dsdx, &
         this%c_Xh%dtdx, this%c_Xh)
    call cdtp(gx%x, this%v%x, this%c_Xh%drdy, this%c_Xh%dsdy, &
         this%c_Xh%dtdy, this%c_Xh)
    call cdtp(gy%x, this%w%x, this%c_Xh%drdz, this%c_Xh%dsdz, &
         this%c_Xh%dtdz, this%c_Xh)
    call add2(div_Xh%x, gx%x, n_x)
    call add2(div_Xh%x, gy%x, n_x)

    call this%prs_interp%map(this%p_res%x, div_Xh%x, this%msh%nelv, this%Yh)
    call field_cmult(this%p_res, -(a0 / dt), n_y)
    if (.not. this%prs_dirichlet) then
      call ortho(this%p_res%x, this%glb_prs_points, n_y)
    end if
    call this%gs_Yh%op(this%p_res, GS_OP_ADD)

    call this%proj_prs%pre_solving(this%p_res%x, time%tstep, this%c_Yh, n_y, &
         dt_controller, Ax = this%Ax_prs, gs_h = this%gs_Yh, &
         bclst = this%bclst_dp, string = 'Pressure')

    call this%pc_prs%update()
    ksp_results(1) = this%ksp_prs%solve(this%Ax_prs, this%dp, this%p_res%x, &
         n_y, this%c_Yh, this%bclst_dp, this%gs_Yh)
    ksp_results(1)%name = 'Pressure'

    call this%proj_prs%post_solving(this%dp%x, this%Ax_prs, this%c_Yh, &
         this%bclst_dp, this%gs_Yh, n_y, time%tstep, dt_controller)

    call field_add2(this%p, this%dp, n_y)
    if (.not. this%prs_dirichlet) then
      call ortho(this%p%x, this%glb_prs_points, n_y)
    end if

    call this%prs_interp%map(p_Xh%x, this%dp%x, this%msh%nelv, this%Xh)
    call opgrad(gx%x, gy%x, gz%x, p_Xh%x, this%c_Xh)
    call col2(gx%x, this%c_Xh%Binv, n_x)
    call col2(gy%x, this%c_Xh%Binv, n_x)
    call col2(gz%x, this%c_Xh%Binv, n_x)
    call opadd2cm(this%u%x, this%v%x, this%w%x, gx%x, gy%x, gz%x, dt / a0, &
         n_x, this%msh%gdim)

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

    call pnpn2_prs_ax_init(this%Xh, this%Yh, this%c_Xh, this%prs_interp)
  end subroutine fluid_pnpn2_restart

  !> Set up the supported boundary-condition configuration for milestone 1.
  subroutine fluid_pnpn2_setup_bcs(this, user, params)
    class(fluid_pnpn2_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: i

    if (params%valid_path('case.fluid.boundary_conditions')) then
      call neko_error('pnpn2 milestone-1 does not support boundary conditions.')
    end if

    do i = 1, size(this%msh%labeled_zones)
      if (this%msh%labeled_zones(i)%size .gt. 0) then
        call neko_error('pnpn2 milestone-1 currently supports periodic cases ' // &
             'only.')
      end if
    end do

    call this%bcs_vel%init()
    call this%bcs_prs%init()
    call this%bclst_dp%init()
    call this%bclst_du%init()
    call this%bclst_dv%init()
    call this%bclst_dw%init()

    this%prs_dirichlet = .false.
  end subroutine fluid_pnpn2_setup_bcs

end module fluid_pnpn2
