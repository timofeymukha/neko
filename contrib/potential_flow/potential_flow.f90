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

!> Initialize an incompressible flow field from a velocity potential.
program potential_flow
  use neko
  use bc, only : bc_t
  use neumann, only : neumann_t
  use inflow, only : inflow_t
  use no_slip, only : no_slip_t
  use symmetry, only : symmetry_t
  use blasius, only : blasius_t
  use expression_dirichlet_vector, only : expression_dirichlet_vector_t
  use json_module, only : json_core, json_value
  use json_utils, only : json_get_or_lookup_or_default
  use registry, only : neko_const_registry
  use math, only : NEKO_EPS
  use num_types, only : i8
  use operators, only : rotate_cyc
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use phmg, only : phmg_t
  use precon, only : pc_t, precon_allocator, precon_destroy
  use field_math, only : field_glsum
  use vector_math, only : vector_vdot3, vector_masked_gather_copy_0, &
       vector_face_masked_gather_copy_0
  use opr_device, only : device_ortho
  use device, only : DEVICE_TO_HOST, device_sync
  use device_math, only : device_absval, device_add2, device_cfill, &
       device_col2, device_copy, device_glmax, device_glsc3, device_glsum
  use fld_file, only : fld_file_t
  use mpi_f08
  implicit none

  type :: potential_boundary_t
     character(len=:), allocatable :: type_name
     integer, allocatable :: zone_indices(:)
     logical :: neumann = .false.
     logical :: dirichlet = .false.
     type(neumann_t) :: flux_bc
     type(zero_dirichlet_t) :: potential_bc
     class(bc_t), allocatable :: velocity_bc
  end type potential_boundary_t

  character(len=NEKO_FNAME_LEN) :: config_fname, case_fname
  character(len=NEKO_FNAME_LEN) :: mesh_fname
  character(len=NEKO_FNAME_LEN) :: output_fname
  character(len=:), allocatable :: case_setting, mesh_setting
  character(len=:), allocatable :: output_setting, scheme
  character(len=:), allocatable :: solver_type, precon_type
  character(len=:), allocatable :: boundary_path
  character(len=LOG_SIZE) :: log_buf
  character(len=80) :: suffix
  type(json_file) :: config, params, time_params, precon_params
  type(mesh_t) :: msh
  type(file_t) :: mesh_file, output_file
  type(space_t) :: Xh
  type(dofmap_t), target :: dm
  type(gs_t), target :: gs
  type(coef_t), target :: coef
  type(field_t), target :: phi, rhs, p, u, v, w, div_u, work
  type(field_list_t) :: output_fields
  type(time_state_t) :: time
  type(bc_list_t) :: potential_bcs, flux_bcs, velocity_bcs
  type(potential_boundary_t), allocatable, target :: boundaries(:)
  class(ax_t), allocatable :: Ax
  class(ksp_t), allocatable :: solver
  class(pc_t), allocatable, target :: preconditioner
  type(ksp_monitor_t) :: monitor
  integer :: argc, polynomial_order, lx, n
  integer :: i, max_iter
  integer(kind=i8) :: glb_n_points
  real(kind=rp) :: tolerance, relative_tolerance
  real(kind=rp) :: net_flux, abs_flux
  real(kind=rp) :: div_l2, div_linf
  logical :: cyclic, has_dirichlet, solver_monitor
  logical :: has_relative_tolerance

  call neko_init()
  argc = command_argument_count()

  if (argc .ne. 1) then
     call usage()
     call neko_finalize()
     stop
  end if

  call get_command_argument(1, config_fname)
  call neko_job_info()
  call config%load_file(filename = trim(config_fname))

  call json_get(config, "case_file", case_setting)
  case_fname = relative_to_file(config_fname, case_setting)
  call params%load_file(filename = trim(case_fname))
  call load_case_constants(params)

  call json_get(params, "case.fluid.scheme", scheme)
  if (trim(scheme) .ne. "pnpn") then
     call neko_error("potential_flow supports the incompressible pnpn " // &
          "scheme only")
  end if

  if (config%valid_path("mesh_file")) then
     call json_get(config, "mesh_file", mesh_setting)
     mesh_fname = relative_to_file(config_fname, mesh_setting)
  else
     call json_get(params, "case.mesh_file", mesh_setting)
     mesh_fname = relative_to_file(case_fname, mesh_setting)
  end if

  if (config%valid_path("polynomial_order")) then
     call json_get_or_lookup(config, "polynomial_order", polynomial_order)
  else
     call json_get_or_lookup(params, "case.numerics.polynomial_order", &
          polynomial_order)
  end if
  lx = polynomial_order + 1

  call json_get(params, "case.time", time_params)
  call time%init(time_params)
  if (config%valid_path("evaluation_time")) then
     call json_get_or_lookup(config, "evaluation_time", time%t)
  end if

  call json_get_or_default(config, "output_filename", output_setting, &
       "potential_flow.fld")
  output_fname = relative_to_file(config_fname, output_setting)
  call filename_suffix(output_fname, suffix)
  if (trim(suffix) .ne. "fld") then
     call neko_error("The potential_flow output must have a .fld suffix")
  end if

  call neko_log%section("Potential flow initialization")
  call neko_log%message("Configuration file: " // trim(config_fname))
  call neko_log%message("Case file: " // trim(case_fname))
  call neko_log%message("Mesh file: " // trim(mesh_fname))

  call mesh_file%init(trim(mesh_fname))
  call mesh_file%read(msh)
  call mesh_file%free()

  call Xh%init(GLL, lx, lx, lx)
  call dm%init(msh, Xh)
  call gs%init(dm)
  call coef%init(gs)

  cyclic = .false.
  call json_get_or_default(params, "case.fluid.cyclic", cyclic, .false.)
  if (config%valid_path("cyclic")) then
     call json_get(config, "cyclic", cyclic)
  end if
  coef%cyclic = cyclic
  call coef%generate_cyclic_bc()

  call phi%init(dm, "potential")
  call rhs%init(dm, "potential_rhs")
  call p%init(dm, "p")
  call u%init(dm, "u")
  call v%init(dm, "v")
  call w%init(dm, "w")
  call div_u%init(dm, "div_u")
  call work%init(dm, "work")
  n = dm%size()
  glb_n_points = int(msh%glb_nelv, i8) * int(Xh%lxyz, i8)

  if (config%valid_path("boundary_conditions")) then
     boundary_path = "boundary_conditions"
     call neko_log%message("Boundary conditions: configuration override")
     call setup_boundaries(config, boundary_path, coef, time, boundaries, &
          potential_bcs, flux_bcs, velocity_bcs, has_dirichlet)
  else
     boundary_path = "case.fluid.boundary_conditions"
     call neko_log%message("Boundary conditions: referenced case")
     call setup_boundaries(params, boundary_path, coef, time, boundaries, &
          potential_bcs, flux_bcs, velocity_bcs, has_dirichlet)
  end if

  rhs = 0.0_rp
  call flux_bcs%apply_scalar(rhs%x, n, time, strong = .false.)
  net_flux = field_glsum(rhs, n)
  abs_flux = global_absolute_sum(rhs, work, n)

  write(log_buf, '(A,ES13.5)') &
       "Integrated prescribed Neumann flux: ", net_flux
  call neko_log%message(log_buf)

  if (.not. has_dirichlet) then
     tolerance = max(100.0_rp * NEKO_EPS, 1.0e-10_rp) * &
          max(1.0_rp, abs_flux)
     if (abs(net_flux) .gt. tolerance) then
        write(log_buf, '(A,ES13.5,A,ES13.5)') &
             "Incompatible pure-Neumann problem: net flux ", net_flux, &
             ", tolerance ", tolerance
        call neko_error(trim(log_buf))
     end if
     call orthogonalize_field(rhs, glb_n_points, n)
  end if

  call gs%op(rhs, GS_OP_ADD)
  if (NEKO_BCKND_DEVICE .eq. 1) call device_sync()
  call potential_bcs%apply_scalar(rhs%x, n)

  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_cfill(coef%h1_d, 1.0_rp, n)
     call device_cfill(coef%h2_d, 0.0_rp, n)
  else
     coef%h1 = 1.0_rp
     coef%h2 = 0.0_rp
  end if
  coef%ifh2 = .false.
  call ax_helm_factory(Ax, full_formulation = .false.)

  call json_get(config, "solver.type", solver_type)
  call json_get(config, "solver.preconditioner.type", precon_type)
  call json_get(config, "solver.preconditioner", precon_params)
  call json_get_or_lookup(config, "solver.absolute_tolerance", tolerance)
  has_relative_tolerance = config%valid_path("solver.relative_tolerance")
  if (has_relative_tolerance) then
     call json_get_or_lookup(config, "solver.relative_tolerance", &
          relative_tolerance)
  end if
  call json_get_or_lookup_or_default(config, "solver.max_iterations", &
       max_iter, 5000)
  call json_get_or_default(config, "solver.monitor", solver_monitor, .false.)

  call initialize_preconditioner(preconditioner, precon_type, coef, &
       potential_bcs, precon_params, dm, gs)
  call krylov_solver_allocator(solver, solver_type)
  if (has_relative_tolerance) then
     call solver%init(n, max_iter, M = preconditioner, &
          rel_tol = relative_tolerance, abs_tol = tolerance, &
          monitor = solver_monitor)
  else
     call solver%init(n, max_iter, M = preconditioner, &
          abs_tol = tolerance, monitor = solver_monitor)
  end if

  call neko_log%message("Potential solver: " // trim(solver_type) // &
       " with " // trim(precon_type))
  write(log_buf, '(A,ES13.5)') "Potential solver tolerance: ", tolerance
  call neko_log%message(log_buf)
  if (has_relative_tolerance) then
     write(log_buf, '(A,ES13.5)') "Potential relative tolerance: ", &
          relative_tolerance
     call neko_log%message(log_buf)
  end if
  write(log_buf, '(A,I0)') "Potential solver maximum iterations: ", max_iter
  call neko_log%message(log_buf)

  phi = 0.0_rp
  monitor = solver%solve(Ax, phi, rhs%x, n, coef, potential_bcs, gs)
  write(log_buf, '(A,I0,A,ES13.5)') "Potential solve: ", monitor%iter, &
       " iterations, residual ", monitor%res_final
  call neko_log%message(log_buf)
  if (.not. monitor%converged) then
     call neko_warning("Potential solve did not converge; " // &
          "writing the current iterate")
  end if
  if (.not. has_dirichlet) then
     call orthogonalize_field(phi, glb_n_points, n)
  end if

  call opgrad(u%x, v%x, w%x, phi%x, coef)
  call rotate_cyc(u, v, w, 1, coef)
  call gs%op(u%x, v%x, w%x, n, GS_OP_ADD)
  if (NEKO_BCKND_DEVICE .eq. 1) call device_sync()
  call rotate_cyc(u, v, w, 0, coef)
  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_col2(u%x_d, coef%Binv_d, n)
     call device_col2(v%x_d, coef%Binv_d, n)
     call device_col2(w%x_d, coef%Binv_d, n)
  else
     call col2(u%x, coef%Binv, n)
     call col2(v%x, coef%Binv, n)
     call col2(w%x, coef%Binv, n)
  end if

  call divergence_norms(div_l2, div_linf, div_u, work, u, v, w, &
       coef, n)
  write(log_buf, '(A,ES13.5,A,ES13.5)') &
       "Projected velocity divergence: L2 ", div_l2, ", Linf ", div_linf
  call neko_log%message(log_buf)

  call velocity_bcs%apply_vector_field(u, v, w, time = time, strong = .true.)
  call divergence_norms(div_l2, div_linf, div_u, work, u, v, w, &
       coef, n)
  write(log_buf, '(A,ES13.5,A,ES13.5)') &
       "After velocity BCs: divergence L2 ", div_l2, ", Linf ", div_linf
  call neko_log%message(log_buf)

  p = 0.0_rp
  call output_fields%init(4)
  call output_fields%assign_to_field(1, p)
  call output_fields%assign_to_field(2, u)
  call output_fields%assign_to_field(3, v)
  call output_fields%assign_to_field(4, w)
  call output_file%init(trim(output_fname), precision = rp)
  select type (fld => output_file%file_type)
  type is (fld_file_t)
     fld%skip_pressure = .false.
     fld%skip_velocity = .false.
  end select
  if (NEKO_BCKND_DEVICE .eq. 1) then
     call output_fields%copy_from(DEVICE_TO_HOST, .true.)
  end if
  call output_file%write(output_fields, time%t)
  call output_file%free()
  call output_fields%free()

  call neko_log%message("Wrote initialized p, u, v, w fields to " // &
       trim(output_fname))
  call neko_log%end_section()

  call potential_bcs%free()
  call flux_bcs%free()
  call velocity_bcs%free()
  do i = 1, size(boundaries)
     if (boundaries(i)%neumann) call boundaries(i)%flux_bc%free()
     if (boundaries(i)%dirichlet) call boundaries(i)%potential_bc%free()
     if (allocated(boundaries(i)%velocity_bc)) then
        call boundaries(i)%velocity_bc%free()
        deallocate(boundaries(i)%velocity_bc)
     end if
  end do
  deallocate(boundaries)

  call solver%free()
  deallocate(solver)
  call precon_destroy(preconditioner)
  deallocate(preconditioner)
  deallocate(Ax)
  call work%free()
  call div_u%free()
  call w%free()
  call v%free()
  call u%free()
  call p%free()
  call rhs%free()
  call phi%free()
  call coef%free()
  call gs%free()
  call dm%free()
  call Xh%free()
  call msh%free()
  call neko_finalize()

contains

  !> Print command-line usage.
  subroutine usage()
    if (pe_rank .eq. 0) then
       write(*,*) "Usage: potential_flow <configuration.json>"
    end if
  end subroutine usage

  !> Resolve a path relative to the file in which it was configured.
  function relative_to_file(source_file, path) result(resolved)
    character(len=*), intent(in) :: source_file
    character(len=*), intent(in) :: path
    character(len=NEKO_FNAME_LEN) :: resolved
    integer :: slash

    if (len_trim(path) .gt. 0 .and. path(1:1) .eq. "/") then
       resolved = trim(path)
       return
    end if

    slash = scan(trim(source_file), "/", back = .true.)
    if (slash .gt. 0) then
       resolved = source_file(1:slash) // trim(path)
    else
       resolved = trim(path)
    end if
  end function relative_to_file

  !> Populate the constant registry used by json_get_or_lookup.
  subroutine load_case_constants(case_params)
    type(json_file), intent(inout) :: case_params
    type(json_file) :: constant
    type(vector_t), pointer :: vector_value
    character(len=:), allocatable :: name
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: real_value
    integer :: integer_value
    integer :: i, n_constants, value_type
    logical :: found

    if (.not. case_params%valid_path("case.constants")) return

    call case_params%info("case.constants", n_children = n_constants)
    do i = 1, n_constants
       call json_extract_item(case_params, "case.constants", i, constant)
       call json_get(constant, "name", name)
       call constant%info("value", found = found, var_type = value_type)
       if (.not. found) then
          call neko_error("Missing value for case constant " // name)
       end if

       select case (value_type)
       case (5)
          call json_get(constant, "value", integer_value)
          call neko_const_registry%add_integer_scalar(integer_value, name)
       case (6)
          call json_get(constant, "value", real_value)
          call neko_const_registry%add_real_scalar(real_value, name)
       case (3)
          call json_get(constant, "value", values)
          call neko_const_registry%add_vector(size(values), name)
          vector_value => neko_const_registry%get_vector(name)
          vector_value%x = values
       case default
          call neko_error("Unsupported value in case.constants for " // name)
       end select
    end do
  end subroutine load_case_constants

  !> Construct potential and velocity boundary conditions from JSON.
  subroutine setup_boundaries(boundary_params, path, c, time_state, boundary, &
       phi_list, flux_list, velocity_list, any_dirichlet)
    type(json_file), intent(inout) :: boundary_params
    character(len=*), intent(in) :: path
    type(coef_t), target, intent(in) :: c
    type(time_state_t), intent(in) :: time_state
    type(potential_boundary_t), allocatable, target, intent(out) :: boundary(:)
    type(bc_list_t), intent(inout) :: phi_list
    type(bc_list_t), intent(inout) :: flux_list
    type(bc_list_t), intent(inout) :: velocity_list
    logical, intent(out) :: any_dirichlet
    type(json_core) :: core
    type(json_value), pointer :: bc_array
    type(json_file) :: bc_json
    logical, allocatable :: marked_zones(:)
    logical :: found
    integer :: i, j, n_bcs, zone

    n_bcs = 0
    if (boundary_params%valid_path(path)) then
       call boundary_params%info(path, n_children = n_bcs)
    end if

    allocate(boundary(n_bcs))
    allocate(marked_zones(size(c%msh%labeled_zones)))
    marked_zones = .false.
    call phi_list%init(max(1, n_bcs))
    call flux_list%init(max(1, n_bcs))
    call velocity_list%init(max(1, n_bcs))
    any_dirichlet = .false.

    if (n_bcs .eq. 0) then
       call neko_warning("No fluid boundary conditions found; the " // &
            "potential-flow velocity will be zero")
       deallocate(marked_zones)
       return
    end if

    call boundary_params%get_core(core)
    call boundary_params%get(path, bc_array, found)
    if (.not. found) then
       call neko_error("Unable to read potential-flow boundary conditions")
    end if

    do i = 1, n_bcs
       call json_extract_item(core, bc_array, i, bc_json)
       call json_get(bc_json, "type", boundary(i)%type_name)
       call json_get_or_lookup(bc_json, "zone_indices", &
            boundary(i)%zone_indices)

       do j = 1, size(boundary(i)%zone_indices)
          zone = boundary(i)%zone_indices(j)
          if (zone .lt. 1 .or. zone .gt. size(marked_zones)) then
             call neko_error("Invalid boundary zone index in potential_flow")
          end if
          if (marked_zones(zone)) then
             call neko_error("A boundary zone is assigned more than once")
          end if
          marked_zones(zone) = .true.
       end do

       call configure_boundary(boundary(i), bc_json)

       if (boundary(i)%neumann) then
          call boundary(i)%flux_bc%init_from_components(c, 0.0_rp)
          call mark_zones(boundary(i)%flux_bc, c, &
               boundary(i)%zone_indices)
          call boundary(i)%flux_bc%finalize(.true.)
          call flux_list%append(boundary(i)%flux_bc)
       end if

       if (boundary(i)%dirichlet) then
          call boundary(i)%potential_bc%init_from_components(c)
          call mark_zones(boundary(i)%potential_bc, c, &
               boundary(i)%zone_indices)
          call boundary(i)%potential_bc%finalize()
          call phi_list%append(boundary(i)%potential_bc)
          any_dirichlet = .true.
       end if

       if (allocated(boundary(i)%velocity_bc)) then
          call boundary(i)%velocity_bc%init(c, bc_json)
          call mark_zones(boundary(i)%velocity_bc, c, &
               boundary(i)%zone_indices)
          call boundary(i)%velocity_bc%finalize()
          call velocity_list%append(boundary(i)%velocity_bc)
          call set_normal_flux(boundary(i), c, time_state)
       end if
    end do

    deallocate(marked_zones)
  end subroutine setup_boundaries

  !> Select how a Neko velocity boundary maps to a potential boundary.
  subroutine configure_boundary(boundary, json)
    type(potential_boundary_t), intent(inout) :: boundary
    type(json_file), intent(inout) :: json
    logical :: moving

    select case (trim(boundary%type_name))
    case ("velocity_value")
       boundary%neumann = .true.
       allocate(inflow_t :: boundary%velocity_bc)
    case ("expression_velocity")
       boundary%neumann = .true.
       allocate(expression_dirichlet_vector_t :: boundary%velocity_bc)
    case ("blasius_profile")
       boundary%neumann = .true.
       allocate(blasius_t :: boundary%velocity_bc)
    case ("no_slip")
       moving = .false.
       call json_get_or_default(json, "moving", moving, .false.)
       if (moving) then
          call neko_error("Moving no_slip walls are not supported")
       end if
       boundary%neumann = .true.
       allocate(no_slip_t :: boundary%velocity_bc)
    case ("symmetry")
       boundary%neumann = .true.
       allocate(symmetry_t :: boundary%velocity_bc)
    case ("shear_stress", "wall_model")
       boundary%neumann = .true.
    case ("outflow", "normal_outflow", "outflow+dong", &
         "normal_outflow+dong", "outflow+user", &
         "normal_outflow+user", "expression_pressure", "user_pressure")
       boundary%dirichlet = .true.
    case ("user_velocity", "overset_interface")
       call neko_error("Boundary type " // trim(boundary%type_name) // &
            " requires user data and is not supported by potential_flow")
    case default
       call neko_error("Unsupported potential_flow boundary type: " // &
            trim(boundary%type_name))
    end select
  end subroutine configure_boundary

  !> Mark all zones assigned to a boundary condition.
  subroutine mark_zones(bc_object, c, zones)
    class(bc_t), intent(inout) :: bc_object
    type(coef_t), intent(in) :: c
    integer, intent(in) :: zones(:)
    integer :: i

    do i = 1, size(zones)
       call bc_object%mark_zone(c%msh%labeled_zones(zones(i)))
    end do
    bc_object%zone_indices = zones
  end subroutine mark_zones

  !> Initialize a Neko preconditioner from the standalone solver settings.
  subroutine initialize_preconditioner(pc, type_name, c, bcs, pc_params, &
       dof, gather_scatter)
    class(pc_t), allocatable, target, intent(inout) :: pc
    character(len=*), intent(in) :: type_name
    type(coef_t), target, intent(in) :: c
    type(bc_list_t), target, intent(inout) :: bcs
    type(json_file), intent(inout) :: pc_params
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gather_scatter

    call precon_allocator(pc, type_name)
    select case (trim(type_name))
    case ("jacobi")
       select type (jacobi => pc)
       type is (jacobi_t)
          call jacobi%init(c, dof, gather_scatter)
       type is (sx_jacobi_t)
          call jacobi%init(c, dof, gather_scatter)
       type is (device_jacobi_t)
          call jacobi%init(c, dof, gather_scatter)
       class default
          call neko_error("Unable to initialize the Jacobi preconditioner")
       end select
    case ("hsmg")
       select type (multigrid => pc)
       type is (hsmg_t)
          call multigrid%init(c, bcs, pc_params)
       class default
          call neko_error("Unable to initialize the HSMG preconditioner")
       end select
    case ("phmg")
       select type (multigrid => pc)
       type is (phmg_t)
          call multigrid%init(c, bcs, pc_params)
       class default
          call neko_error("Unable to initialize the PHMG preconditioner")
       end select
    case ("ident")
       continue
    case default
       call neko_error("Unsupported potential_flow preconditioner: " // &
            trim(type_name))
    end select
  end subroutine initialize_preconditioner

  !> Remove the constant null-space component on the active backend.
  subroutine orthogonalize_field(x, global_size, local_size)
    type(field_t), intent(inout) :: x
    integer(kind=i8), intent(in) :: global_size
    integer, intent(in) :: local_size

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_ortho(x%x_d, global_size, local_size)
    else
       call ortho(x%x, global_size, local_size)
    end if
  end subroutine orthogonalize_field

  !> Evaluate a prescribed velocity and set its normal potential flux.
  subroutine set_normal_flux(boundary, c, time_state)
    type(potential_boundary_t), intent(inout) :: boundary
    type(coef_t), target, intent(in) :: c
    type(time_state_t), intent(in) :: time_state
    type(field_t) :: bx, by, bz
    type(vector_t) :: bx_surface, by_surface, bz_surface, normal_flux
    type(vector_t) :: normal_x, normal_y, normal_z
    integer :: n_boundary, n_local

    if (.not. boundary%neumann) return

    call bx%init(c%dof, "boundary_u")
    call by%init(c%dof, "boundary_v")
    call bz%init(c%dof, "boundary_w")
    bx = 0.0_rp
    by = 0.0_rp
    bz = 0.0_rp

    call boundary%velocity_bc%apply_vector_generic(bx, by, bz, &
         time = time_state, strong = .true.)

    n_local = c%dof%size()
    n_boundary = boundary%flux_bc%msk(0)
    call bx_surface%init(n_boundary)
    call by_surface%init(n_boundary)
    call bz_surface%init(n_boundary)
    call normal_x%init(n_boundary)
    call normal_y%init(n_boundary)
    call normal_z%init(n_boundary)
    call normal_flux%init(n_boundary)

    if (n_boundary .gt. 0) then
       call vector_masked_gather_copy_0(bx_surface, bx%x, &
            boundary%flux_bc%msk, n_local, n_boundary)
       call vector_masked_gather_copy_0(by_surface, by%x, &
            boundary%flux_bc%msk, n_local, n_boundary)
       call vector_masked_gather_copy_0(bz_surface, bz%x, &
            boundary%flux_bc%msk, n_local, n_boundary)
       call vector_face_masked_gather_copy_0(normal_x, c%nx, &
            boundary%flux_bc%msk, boundary%flux_bc%facet, &
            c%Xh%lx, c%Xh%ly, c%Xh%lz, n_boundary)
       call vector_face_masked_gather_copy_0(normal_y, c%ny, &
            boundary%flux_bc%msk, boundary%flux_bc%facet, &
            c%Xh%lx, c%Xh%ly, c%Xh%lz, n_boundary)
       call vector_face_masked_gather_copy_0(normal_z, c%nz, &
            boundary%flux_bc%msk, boundary%flux_bc%facet, &
            c%Xh%lx, c%Xh%ly, c%Xh%lz, n_boundary)
       call vector_vdot3(normal_flux, bx_surface, by_surface, bz_surface, &
            normal_x, normal_y, normal_z)
    end if
    call boundary%flux_bc%set_flux(normal_flux, 1)

    if (NEKO_BCKND_DEVICE .eq. 1) call device_sync()
    call normal_flux%free()
    call normal_z%free()
    call normal_y%free()
    call normal_x%free()
    call bz_surface%free()
    call by_surface%free()
    call bx_surface%free()
    call bz%free()
    call by%free()
    call bx%free()
  end subroutine set_normal_flux

  !> Compute global L2 and maximum norms of velocity divergence.
  subroutine divergence_norms(l2_norm, max_norm, divergence, scratch, &
       ux, uy, uz, c, n_local)
    real(kind=rp), intent(out) :: l2_norm, max_norm
    type(field_t), intent(inout) :: divergence, scratch
    type(field_t), intent(in) :: ux, uy, uz
    type(coef_t), intent(in) :: c
    integer, intent(in) :: n_local
    real(kind=rp) :: local_max
    integer :: ierr

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call dudxyz(divergence%x_d, ux%x_d, c%drdx_d, c%dsdx_d, &
            c%dtdx_d, c)
       call dudxyz(scratch%x_d, uy%x_d, c%drdy_d, c%dsdy_d, c%dtdy_d, c)
       call device_add2(divergence%x_d, scratch%x_d, n_local)
       call dudxyz(scratch%x_d, uz%x_d, c%drdz_d, c%dsdz_d, c%dtdz_d, c)
       call device_add2(divergence%x_d, scratch%x_d, n_local)

       l2_norm = sqrt(device_glsc3(divergence%x_d, c%B_d, &
            divergence%x_d, n_local) / c%volume)
       call device_copy(scratch%x_d, divergence%x_d, n_local)
       call device_absval(scratch%x_d, n_local)
       max_norm = device_glmax(scratch%x_d, n_local)
    else
       call dudxyz(divergence%x, ux%x, c%drdx, c%dsdx, c%dtdx, c)
       call dudxyz(scratch%x, uy%x, c%drdy, c%dsdy, c%dtdy, c)
       call add2(divergence%x, scratch%x, n_local)
       call dudxyz(scratch%x, uz%x, c%drdz, c%dsdz, c%dtdz, c)
       call add2(divergence%x, scratch%x, n_local)

       l2_norm = sqrt(glsc3(divergence%x, c%B, divergence%x, n_local) / &
            c%volume)
       local_max = maxval(abs(divergence%x))
       call MPI_Allreduce(local_max, max_norm, 1, MPI_REAL_PRECISION, &
            MPI_MAX, NEKO_COMM, ierr)
    end if
  end subroutine divergence_norms

  !> Compute the global one-norm of a field.
  function global_absolute_sum(x, scratch, n_local) result(value)
    integer, intent(in) :: n_local
    type(field_t), intent(in) :: x
    type(field_t), intent(inout) :: scratch
    real(kind=rp) :: value, local_value
    integer :: ierr

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(scratch%x_d, x%x_d, n_local)
       call device_absval(scratch%x_d, n_local)
       value = device_glsum(scratch%x_d, n_local)
    else
       local_value = sum(abs(x%x))
       call MPI_Allreduce(local_value, value, 1, MPI_REAL_PRECISION, &
            MPI_SUM, NEKO_COMM, ierr)
    end if
  end function global_absolute_sum

end program potential_flow
