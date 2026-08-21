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
  use registry, only : neko_const_registry
  use math, only : NEKO_EPS
  use num_types, only : i8
  use operators, only : rotate_cyc
  use jacobi, only : jacobi_t
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

  character(len=NEKO_FNAME_LEN) :: case_fname
  character(len=NEKO_FNAME_LEN) :: mesh_fname
  character(len=NEKO_FNAME_LEN) :: output_fname
  character(len=:), allocatable :: mesh_from_case
  character(len=:), allocatable :: scheme
  character(len=LOG_SIZE) :: log_buf
  character(len=80) :: suffix
  type(json_file) :: params, time_params
  type(mesh_t) :: msh
  type(file_t) :: mesh_file, output_file
  type(space_t) :: Xh
  type(dofmap_t), target :: dm
  type(gs_t), target :: gs
  type(coef_t), target :: coef
  type(jacobi_t), target :: preconditioner
  type(field_t), target :: phi, rhs, p, u, v, w, div_u, work
  type(field_list_t) :: output_fields
  type(time_state_t) :: time
  type(bc_list_t) :: potential_bcs, flux_bcs, velocity_bcs
  type(potential_boundary_t), allocatable, target :: boundaries(:)
  class(ax_t), allocatable :: Ax
  class(ksp_t), allocatable :: solver
  type(ksp_monitor_t) :: monitor
  integer :: argc, polynomial_order, lx, n
  integer :: i, max_iter
  integer(kind=i8) :: glb_n_points
  real(kind=rp) :: tolerance
  real(kind=rp) :: net_flux, abs_flux
  real(kind=rp) :: div_l2, div_linf
  logical :: cyclic, has_dirichlet

  call neko_init()
  argc = command_argument_count()

  if (argc .lt. 1 .or. argc .gt. 2) then
     call usage()
     call neko_finalize()
     stop
  end if

  if (NEKO_BCKND_DEVICE .eq. 1) then
     call neko_error("potential_flow currently supports CPU backends only")
  end if

  call get_command_argument(1, case_fname)
  output_fname = "potential_flow.fld"
  if (argc .eq. 2) call get_command_argument(2, output_fname)

  call filename_suffix(output_fname, suffix)
  if (trim(suffix) .ne. "fld") then
     call neko_error("The potential_flow output must have a .fld suffix")
  end if

  call neko_job_info()
  call params%load_file(filename = trim(case_fname))
  call load_case_constants(params)

  call json_get(params, "case.fluid.scheme", scheme)
  if (trim(scheme) .ne. "pnpn") then
     call neko_error("potential_flow supports the incompressible pnpn " // &
          "scheme only")
  end if

  call json_get(params, "case.mesh_file", mesh_from_case)
  mesh_fname = relative_to_case(case_fname, mesh_from_case)
  call json_get_or_lookup(params, "case.numerics.polynomial_order", &
       polynomial_order)
  lx = polynomial_order + 1

  call json_get(params, "case.time", time_params)
  call time%init(time_params)

  call neko_log%section("Potential flow initialization")
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

  call setup_boundaries(params, coef, time, boundaries, potential_bcs, &
       flux_bcs, velocity_bcs, has_dirichlet)

  rhs%x = 0.0_rp
  call flux_bcs%apply_scalar(rhs%x, n, time, strong = .false.)
  net_flux = glsum(rhs%x, n)
  abs_flux = global_absolute_sum(rhs%x, n)

  write(log_buf, '(A,ES13.5)') "Net prescribed boundary flux: ", net_flux
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
     call ortho(rhs%x, glb_n_points, n)
  end if

  call gs%op(rhs, GS_OP_ADD)
  call potential_bcs%apply_scalar(rhs%x, n)

  coef%h1 = 1.0_rp
  coef%h2 = 0.0_rp
  coef%ifh2 = .false.
  call ax_helm_factory(Ax, full_formulation = .false.)

  max_iter = 5000
  tolerance = 1.0e-8_rp
  call preconditioner%init(coef, dm, gs)
  call krylov_solver_factory(solver, n, "cg", max_iter, &
       abstol = tolerance, M = preconditioner, monitor = .false.)

  phi%x = 0.0_rp
  monitor = solver%solve(Ax, phi, rhs%x, n, coef, potential_bcs, gs)
  write(log_buf, '(A,I0,A,ES13.5)') "Potential solve: ", monitor%iter, &
       " iterations, residual ", monitor%res_final
  call neko_log%message(log_buf)
  if (.not. monitor%converged) then
     call neko_error("Potential solve did not converge")
  end if
  if (.not. has_dirichlet) call ortho(phi%x, glb_n_points, n)

  call opgrad(u%x, v%x, w%x, phi%x, coef)
  call rotate_cyc(u, v, w, 1, coef)
  call gs%op(u%x, v%x, w%x, n, GS_OP_ADD)
  call rotate_cyc(u, v, w, 0, coef)
  call col2(u%x, coef%Binv, n)
  call col2(v%x, coef%Binv, n)
  call col2(w%x, coef%Binv, n)

  call divergence_norms(div_l2, div_linf, div_u, work, u, v, w, &
       coef, n)
  write(log_buf, '(A,ES13.5,A,ES13.5)') &
       "Projected velocity divergence: L2 ", div_l2, ", Linf ", div_linf
  call neko_log%message(log_buf)

  call velocity_bcs%apply_vector(u%x, v%x, w%x, n, time, &
       strong = .true.)
  call divergence_norms(div_l2, div_linf, div_u, work, u, v, w, &
       coef, n)
  write(log_buf, '(A,ES13.5,A,ES13.5)') &
       "After velocity BCs: divergence L2 ", div_l2, ", Linf ", div_linf
  call neko_log%message(log_buf)

  p%x = 0.0_rp
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
  call preconditioner%free()
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
       write(*,*) "Usage: potential_flow <case file> [output.fld]"
    end if
  end subroutine usage

  !> Resolve a mesh path relative to the case-file directory.
  function relative_to_case(case_file, path) result(resolved)
    character(len=*), intent(in) :: case_file
    character(len=*), intent(in) :: path
    character(len=NEKO_FNAME_LEN) :: resolved
    integer :: slash

    if (len_trim(path) .gt. 0 .and. path(1:1) .eq. "/") then
       resolved = trim(path)
       return
    end if

    slash = scan(trim(case_file), "/", back = .true.)
    if (slash .gt. 0) then
       resolved = case_file(1:slash) // trim(path)
    else
       resolved = trim(path)
    end if
  end function relative_to_case

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

  !> Construct potential and velocity boundary conditions from the case.
  subroutine setup_boundaries(case_params, c, time_state, boundary, &
       phi_list, flux_list, velocity_list, any_dirichlet)
    type(json_file), intent(inout) :: case_params
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
    if (case_params%valid_path("case.fluid.boundary_conditions")) then
       call case_params%info("case.fluid.boundary_conditions", &
            n_children = n_bcs)
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

    call case_params%get_core(core)
    call case_params%get("case.fluid.boundary_conditions", bc_array, found)
    if (.not. found) then
       call neko_error("Unable to read fluid boundary conditions")
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

  !> Evaluate a prescribed velocity and set its normal potential flux.
  subroutine set_normal_flux(boundary, c, time_state)
    type(potential_boundary_t), intent(inout) :: boundary
    type(coef_t), intent(in) :: c
    type(time_state_t), intent(in) :: time_state
    type(field_t) :: bx, by, bz
    type(vector_t) :: normal_flux
    real(kind=rp) :: normal(3)
    integer :: i, k, facet, idx(4), n_local

    if (.not. boundary%neumann) return

    call bx%init(c%dof, "boundary_u")
    call by%init(c%dof, "boundary_v")
    call bz%init(c%dof, "boundary_w")
    bx%x = 0.0_rp
    by%x = 0.0_rp
    bz%x = 0.0_rp

    n_local = c%dof%size()
    !$omp parallel
    call boundary%velocity_bc%apply_vector(bx%x, by%x, bz%x, n_local, &
         time_state, strong = .true.)
    !$omp end parallel

    call normal_flux%init(boundary%flux_bc%msk(0))
    do i = 1, boundary%flux_bc%msk(0)
       k = boundary%flux_bc%msk(i)
       facet = boundary%flux_bc%facet(i)
       idx = nonlinear_index(k, c%Xh%lx, c%Xh%ly, c%Xh%lz)
       normal = c%get_normal(idx(1), idx(2), idx(3), idx(4), facet)
       normal_flux%x(i) = bx%x(k, 1, 1, 1) * normal(1) + &
            by%x(k, 1, 1, 1) * normal(2) + &
            bz%x(k, 1, 1, 1) * normal(3)
    end do
    call boundary%flux_bc%set_flux(normal_flux, 1)

    call normal_flux%free()
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
  end subroutine divergence_norms

  !> Compute the global one-norm of an array.
  function global_absolute_sum(x, n_local) result(value)
    integer, intent(in) :: n_local
    real(kind=rp), intent(in) :: x(n_local)
    real(kind=rp) :: value, local_value
    integer :: ierr

    local_value = sum(abs(x))
    call MPI_Allreduce(local_value, value, 1, MPI_REAL_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
  end function global_absolute_sum

end program potential_flow
