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
!> Implements `volume_operation_t`.
module volume_operation
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use time_state, only : time_state_t
  use coefs, only : coef_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use vector, only : vector_t
  use vector_math, only : vector_glsc2, vector_glsum, vector_glmin, &
       vector_glmax
  use device_math, only : device_glsc2, device_glsum, &
        device_glmin, device_glmax, device_masked_gather_copy_aligned
  use math, only : glsc2, glsum, glmin, glmax, masked_gather_copy
  use mask, only : mask_t
  use time_based_controller, only : time_based_controller_t
  use point_zone, only : point_zone_t
  use point_zone_registry, only : neko_point_zone_registry
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM, MPI_INTEGER
  implicit none
  private

  !> A simulation component for volume reductions over all GLL nodes or a point
  !! zone selection.
  type, public, extends(simulation_component_t) :: volume_operation_t
     !> Field to sample.
     type(field_t), pointer :: field => null()
     !> Volume coefficients.
     type(coef_t), pointer :: coef => null()
     !> Name of the selected point zone when a selector is used.
     character(len=:), allocatable :: point_zone_name
     !> Name of the input field.
     character(len=:), allocatable :: field_name
     !> Requested operations.
     character(len=16), allocatable :: operations(:)
     !> Whether a point zone selector is active.
     logical :: use_point_zone = .false.
     !> Whether to compute the integral.
     logical :: compute_integral = .false.
     !> Whether to compute the average.
     logical :: compute_average = .false.
     !> Whether to compute the minimum.
     logical :: compute_min = .false.
     !> Whether to compute the maximum.
     logical :: compute_max = .false.
     !> Whether to write results to the log.
     logical :: log = .true.
     !> Whether CSV output is enabled.
     logical :: csv_output_enabled = .false.
     !> Optional CSV output file.
     type(file_t) :: csv_output
     !> Reusable row buffer for CSV output.
     type(vector_t) :: csv_row
     !> Private point-zone mask for masked gathers.
     type(mask_t), private :: point_zone_mask_
     !> Precomputed point-zone quadrature weights.
     type(vector_t), private :: point_zone_weights_
     !> Gathered field values on the point-zone selection.
     type(vector_t), private :: point_zone_values_
     !> Most recently computed integral.
     real(kind=rp) :: integral = 0.0_rp
     !> Most recently computed average.
     real(kind=rp) :: average = 0.0_rp
     !> Most recently computed minimum.
     real(kind=rp) :: minimum = huge(0.0_rp)
     !> Most recently computed maximum.
     real(kind=rp) :: maximum = -huge(0.0_rp)
   contains
     !> Construct the component from a case-file JSON object.
     procedure, pass(this) :: init => volume_operation_init_from_json
     !> Generic constructor from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Construct from explicit time-based controllers.
     procedure, pass(this) :: init_from_controllers => &
          volume_operation_init_from_controllers
     !> Construct from time-based controller properties.
     procedure, pass(this) :: init_from_controllers_properties => &
          volume_operation_init_from_controllers_properties
     !> Common constructor used by all public constructors.
     procedure, private, pass(this) :: init_common => &
          volume_operation_init_common
     !> Free the component.
     procedure, pass(this) :: free => volume_operation_free
     !> Compute the requested volume operations.
     procedure, pass(this) :: compute_ => volume_operation_compute
  end type volume_operation_t

contains

  !> Construct from JSON.
  !! @param json JSON object describing the simcomp.
  !! @param case Simulation case.
  subroutine volume_operation_init_from_json(this, json, case)
    class(volume_operation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    character(len=:), allocatable :: field_name
    character(len=16), allocatable :: operations(:)
    logical :: log
    character(len=:), allocatable :: output_filename
    character(len=:), allocatable :: point_zone_name

    call this%free()

    call json_get_or_default(json, "name", name, "volume_operation")
    call json_get_or_default(json, "log", log, .true.)
    call this%init_base(json, case)

    call json_get(json, "field_name", field_name)
    call json_get(json, "operations", operations)

    if (json%valid_path("point_zone")) then
       call json_get(json, "point_zone", point_zone_name)
    end if

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", output_filename)
       if (allocated(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log, point_zone_name, output_filename)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log, output_filename = output_filename)
       end if
    else
       if (allocated(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log, point_zone_name)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log)
       end if
    end if

    if (allocated(name)) deallocate(name)
    if (allocated(field_name)) deallocate(field_name)
    if (allocated(operations)) deallocate(operations)
    if (allocated(output_filename)) deallocate(output_filename)
    if (allocated(point_zone_name)) deallocate(point_zone_name)
  end subroutine volume_operation_init_from_json

  !> Common constructor shared by all public constructors.
  !! @param name Unique simcomp name.
  !! @param coef SEM coefficients.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Whether to emit tabular log output.
  !! @param point_zone_name Optional point zone selector.
  !! @param output_filename Optional CSV output filename.
  subroutine volume_operation_init_common(this, name, coef, field_name, &
       operations, log, point_zone_name, output_filename)
    class(volume_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(coef_t), intent(inout), target :: coef
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in) :: log
    character(len=*), intent(in), optional :: point_zone_name
    character(len=*), intent(in), optional :: output_filename
    character(len=:), allocatable :: csv_header
    character(len=LOG_SIZE) :: log_buf
    class(point_zone_t), pointer :: point_zone
    integer :: i
    integer :: ierr
    integer :: n_pts
    integer :: n_pts_global

    point_zone => null()

    this%name = name
    this%log = log
    this%use_point_zone = .false.
    this%compute_integral = .false.
    this%compute_average = .false.
    this%compute_min = .false.
    this%compute_max = .false.
    this%csv_output_enabled = .false.
    this%integral = 0.0_rp
    this%average = 0.0_rp
    this%minimum = huge(0.0_rp)
    this%maximum = -huge(0.0_rp)

    this%coef => coef
    this%field => neko_registry%get_field_by_name(field_name)

    if (size(operations) .eq. 0) then
       call neko_error("volume_operation requires at least one operation")
    end if

    this%field_name = field_name
    allocate(character(len=16) :: this%operations(size(operations)))
    this%operations = operations

    do i = 1, size(this%operations)
       select case (trim(this%operations(i)))
       case ("integral")
          this%compute_integral = .true.
       case ("average")
          this%compute_average = .true.
       case ("min")
          this%compute_min = .true.
       case ("max")
          this%compute_max = .true.
       case default
          call neko_error("volume_operation supports only operations = " // &
               "[integral, average, min, max]")
       end select
    end do

    if (present(point_zone_name)) then
       if (len_trim(point_zone_name) .eq. 0) then
          call neko_error("volume_operation point_zone must not be empty")
       end if

       point_zone => &
            neko_point_zone_registry%get_point_zone(trim(point_zone_name))
       n_pts = point_zone%mask%size()

       n_pts_global = n_pts
       call MPI_Allreduce(MPI_IN_PLACE, n_pts_global, 1, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM, ierr)
       if (n_pts_global .eq. 0) then
          call neko_error("volume_operation point_zone selects no GLL nodes")
       end if

       this%use_point_zone = .true.
       this%point_zone_name = trim(point_zone_name)
       call this%point_zone_mask_%init(point_zone%mask)
       n_pts = this%point_zone_mask_%size()

       call this%point_zone_weights_%init(n_pts)
       call this%point_zone_values_%init(n_pts)
       if (n_pts .gt. 0) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_masked_gather_copy_aligned( &
                  this%point_zone_weights_%x_d, this%coef%B_d, &
                  this%point_zone_mask_%get_d(), this%coef%dof%size(), n_pts)
          else
             call masked_gather_copy(this%point_zone_weights_%x, this%coef%B, &
                  this%point_zone_mask_%get(), this%coef%dof%size(), n_pts)
          end if
       end if
    end if
    n_pts = this%coef%dof%size()

    if (present(output_filename)) then
       csv_header = "tstep,time"
       do i = 1, size(this%operations)
           csv_header = trim(csv_header) // "," // trim(this%operations(i))
       end do
       call this%csv_output%init(trim(output_filename), header = &
             trim(csv_header), overwrite = .true.)
       call this%csv_row%init(2 + size(this%operations))
       this%csv_output_enabled = .true.
    end if

    call neko_log%section("Volume operation")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Field: ", trim(this%field_name)
    call neko_log%message(log_buf)
    call neko_log%message("Operations:")
    do i = 1, size(this%operations)
       write(log_buf, '(A,A)') "  ", trim(this%operations(i))
       call neko_log%message(log_buf)
    end do
    if (this%use_point_zone) then
       write(log_buf, '(A,A)') "Point zone: ", trim(this%point_zone_name)
    else
       write(log_buf, '(A)') "Selection: all GLL nodes"
    end if
    call neko_log%message(log_buf)
    if (this%use_point_zone) then
       n_pts = this%point_zone_mask_%size()
    end if
    n_pts_global = n_pts
    call MPI_Allreduce(MPI_IN_PLACE, n_pts_global, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)
    write(log_buf, '(A,I0)') "Selected volume quadrature points: ", &
         n_pts_global
    call neko_log%message(log_buf)
    call neko_log%end_section()

    if (allocated(csv_header)) deallocate(csv_header)
  end subroutine volume_operation_init_common

  !> Construct from explicit time-based controllers.
  !! @param name Unique simcomp name.
  !! @param case Simulation case owning the simcomp.
  !! @param order Execution priority.
  !! @param preprocess_controller Controller for preprocessing.
  !! @param compute_controller Controller for computation.
  !! @param output_controller Controller for output.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Optional flag controlling log output.
  !! @param point_zone_name Optional point zone selector.
  !! @param output_filename Optional CSV output filename.
  subroutine volume_operation_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, field_name, &
       operations, log, point_zone_name, output_filename)
    class(volume_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in), optional :: log
    character(len=*), intent(in), optional :: point_zone_name
    character(len=*), intent(in), optional :: output_filename
    logical :: log_enabled

    call this%free()

    log_enabled = .true.
    if (present(log)) log_enabled = log

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)

    if (present(output_filename)) then
       if (present(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, point_zone_name, output_filename)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, output_filename = output_filename)
       end if
    else
       if (present(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, point_zone_name)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled)
       end if
    end if
  end subroutine volume_operation_init_from_controllers

  !> Construct from time-based controller properties.
  !! @param name Unique simcomp name.
  !! @param case Simulation case owning the simcomp.
  !! @param order Execution priority.
  !! @param preprocess_control Control mode for preprocessing.
  !! @param preprocess_value Control value for preprocessing.
  !! @param compute_control Control mode for computation.
  !! @param compute_value Control value for computation.
  !! @param output_control Control mode for output.
  !! @param output_value Control value for output.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Optional flag controlling log output.
  !! @param point_zone_name Optional point zone selector.
  !! @param output_filename Optional CSV output filename.
  subroutine volume_operation_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, field_name, operations, &
       log, point_zone_name, output_filename)
    class(volume_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in), optional :: log
    character(len=*), intent(in), optional :: point_zone_name
    character(len=*), intent(in), optional :: output_filename
    logical :: log_enabled

    call this%free()

    log_enabled = .true.
    if (present(log)) log_enabled = log

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)

    if (present(output_filename)) then
       if (present(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, point_zone_name, output_filename)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, output_filename = output_filename)
       end if
    else
       if (present(point_zone_name)) then
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled, point_zone_name)
       else
          call this%init_common(name, case%fluid%c_Xh, field_name, operations, &
               log_enabled)
       end if
    end if
  end subroutine volume_operation_init_from_controllers_properties

  !> Free all resources owned by the component.
  subroutine volume_operation_free(this)
    class(volume_operation_t), intent(inout) :: this

    call this%csv_output%free()
    call this%csv_row%free()
    call this%point_zone_mask_%free()
    call this%point_zone_weights_%free()
    call this%point_zone_values_%free()
    if (allocated(this%point_zone_name)) deallocate(this%point_zone_name)
    if (allocated(this%field_name)) deallocate(this%field_name)
    if (allocated(this%operations)) deallocate(this%operations)

    nullify(this%field)
    nullify(this%coef)

    call this%free_base()
  end subroutine volume_operation_free

  !> Compute and optionally output the requested volume operations.
  !! @param time Current simulation time state.
  subroutine volume_operation_compute(this, time)
    class(volume_operation_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: n_pts
    integer :: i
    real(kind=rp) :: volume
    character(len=18) :: value_buf
    character(len=12) :: step_str
    character(len=:), allocatable :: section_title, header_line, value_line
    integer :: output_col

    n_pts = this%coef%dof%size()
    if (this%use_point_zone) then
       n_pts = this%point_zone_mask_%size()
    end if

    this%integral = 0.0_rp
    this%average = 0.0_rp
    this%minimum = huge(0.0_rp)
    this%maximum = -huge(0.0_rp)
    volume = 0.0_rp

    if (this%use_point_zone) then
       if (n_pts .gt. 0) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_masked_gather_copy_aligned( &
                  this%point_zone_values_%x_d, this%field%x_d, &
                  this%point_zone_mask_%get_d(), this%field%size(), n_pts)
          else
             call masked_gather_copy(this%point_zone_values_%x, this%field%x, &
                  this%point_zone_mask_%get(), this%field%size(), n_pts)
          end if
       end if

       if (this%compute_integral .or. this%compute_average) then
          this%integral = vector_glsc2(this%point_zone_values_, &
               this%point_zone_weights_, n_pts)
       end if

       if (this%compute_average) then
          volume = vector_glsum(this%point_zone_weights_, n_pts)
          if (volume .gt. 0.0_rp) then
             this%average = this%integral / volume
          else
             this%average = 0.0_rp
          end if
       end if

       if (this%compute_min) then
          this%minimum = vector_glmin(this%point_zone_values_, n_pts)
       end if

       if (this%compute_max) then
          this%maximum = vector_glmax(this%point_zone_values_, n_pts)
       end if
    else
       if (this%compute_integral .or. this%compute_average) then
           if (NEKO_BCKND_DEVICE .eq. 1) then
              this%integral = device_glsc2(this%field%x_d, this%coef%B_d, n_pts)
           else
              this%integral = glsc2(this%field%x, this%coef%B, n_pts)
           end if
       end if

       if (this%compute_average) then
           if (NEKO_BCKND_DEVICE .eq. 1) then
              volume = device_glsum(this%coef%B_d, n_pts)
           else
              volume = glsum(this%coef%B, n_pts)
           end if
          if (volume .gt. 0.0_rp) then
             this%average = this%integral / volume
          else
             this%average = 0.0_rp
          end if
       end if

       if (this%compute_min) then
           if (NEKO_BCKND_DEVICE .eq. 1) then
              this%minimum = device_glmin(this%field%x_d, n_pts)
           else
              this%minimum = glmin(this%field%x, n_pts)
           end if
       end if

       if (this%compute_max) then
           if (NEKO_BCKND_DEVICE .eq. 1) then
              this%maximum = device_glmax(this%field%x_d, n_pts)
           else
              this%maximum = glmax(this%field%x, n_pts)
           end if
       end if
    end if

    if (this%log) then
       section_title = trim(this%name)
       call neko_log%section(section_title)

       header_line = repeat(' ', 12) // ' |'
       write(step_str, '(I12)') time%tstep
       step_str = adjustl(step_str)
       value_line = step_str // ' |'

       do i = 1, size(this%operations)
          select case (trim(this%operations(i)))
          case ('integral')
             header_line = header_line // left_pad('Integral:', 18)
             write(value_buf, '(ES18.9)') this%integral
             value_line = value_line // value_buf
          case ('average')
             header_line = header_line // left_pad('Average:', 18)
             write(value_buf, '(ES18.9)') this%average
             value_line = value_line // value_buf
          case ('min')
             header_line = header_line // left_pad('Min:', 18)
             write(value_buf, '(ES18.9)') this%minimum
             value_line = value_line // value_buf
          case ('max')
             header_line = header_line // left_pad('Max:', 18)
             write(value_buf, '(ES18.9)') this%maximum
             value_line = value_line // value_buf
          end select
       end do

       call neko_log%message(header_line)
       call neko_log%message(value_line)
       call neko_log%end_section()
    end if

    if (this%csv_output_enabled) then
       if (this%output_controller%check(time)) then
          this%csv_row%x(1) = real(time%tstep, rp)
          this%csv_row%x(2) = time%t
          output_col = 3
          do i = 1, size(this%operations)
             select case (trim(this%operations(i)))
             case ('integral')
                this%csv_row%x(output_col) = this%integral
             case ('average')
                this%csv_row%x(output_col) = this%average
             case ('min')
                this%csv_row%x(output_col) = this%minimum
             case ('max')
                this%csv_row%x(output_col) = this%maximum
             end select
             output_col = output_col + 1
          end do
          call this%csv_output%write(this%csv_row)
          call this%output_controller%register_execution()
       end if
    end if

    if (allocated(section_title)) deallocate(section_title)
    if (allocated(header_line)) deallocate(header_line)
    if (allocated(value_line)) deallocate(value_line)
  end subroutine volume_operation_compute

  !> Left-pad a string to a fixed width.
  !! @param text Input string.
  !! @param width Width of the padded result.
  pure function left_pad(text, width) result(padded)
    character(len=*), intent(in) :: text
    integer, intent(in) :: width
    character(len=:), allocatable :: padded
    integer :: pad_width

    pad_width = max(0, width - len_trim(text))
    padded = repeat(' ', pad_width) // trim(text)
  end function left_pad

end module volume_operation
