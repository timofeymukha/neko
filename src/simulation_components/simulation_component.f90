! Copyright (c) 2023, The Neko Authors
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
!> Simulation components are objects that encapsulate functionality that can be
!! fit to a particular compute pattern. In particular, each simcomp is executed
!! up to two times per timestep at a given frequency to perfrom some
!! computations, and also has the possibility to perform output at another
!! user-selected frequency. A lot can fit into this pattern, from simple
!! computation of a field, to LES models, to entire solvers.
!! @note
!! The canonical way to abbreviate simulation_component is simcomp.
module simulation_component
  use num_types, only : rp
  use json_module, only : json_file, json_core, json_value
  use case, only : case_t
  use time_based_controller, only : time_based_controller_t
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use utils, only : neko_error
  use logger, only : neko_log
  implicit none
  private

  !> Base abstract class for simulation components.
  type, abstract, public :: simulation_component_t
     !> Pointer to the simulation case.
     type(case_t), pointer :: case
     !> Controller for when to run `preprocess`.
     type(time_based_controller_t) :: preprocess_controller
     !> Controller for when to run `compute`.
     type(time_based_controller_t) :: compute_controller
     !> Controller for when to do output.
     type(time_based_controller_t) :: output_controller
     !> The execution order, lowest excutes first.
     integer :: order
     !> The unique name of the simcomp.
     character(len=20) :: name
     !> The simcomp_executor_t this simcomp belongs to.
     type(simcomp_executor_t), pointer, private :: executor_ => null()
   contains
     !> Constructor for the simulation_component_t (base) class from JSON.
     procedure, pass(this) :: init_base => simulation_component_init_base
     !> Constructor for the simulation_component_t (base) class from components.
     procedure, pass(this) :: init_base_from_components => &
          simulation_component_init_base_from_components
     !> Destructor for the simulation_component_t (base) class.
     procedure, pass(this) :: free_base => simulation_component_free_base
     !> Wrapper for calling `set_counter` for the time based controllers.
     !! Serves as the public interface.
     procedure, pass(this) :: restart => simulation_component_restart_wrapper
     !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: preprocess => &
          simulation_component_preprocess_wrapper
     !> Wrapper for calling `compute_` based on the `compute_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: compute => simulation_component_compute_wrapper
     !> The common constructor using a JSON dictionary.
     procedure(simulation_component_init), pass(this), deferred :: init
     !> Destructor.
     procedure(simulation_component_free), pass(this), deferred :: free
     !> The preprocessing function to be executed during the run.
     procedure, pass(this) :: preprocess_
     !> The main function to be executed during the run.
     procedure, pass(this) :: compute_
     !> The restart function to be called upon restarting simulation
     procedure, pass(this) :: restart_
  end type simulation_component_t

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: simulation_component_alloc_t
     class(simulation_component_t), allocatable :: simcomp
  end type simulation_component_alloc_t

  !> Singleton type that serves as a driver for the simulation components.
  !! Stores all the components in the case and provides an interface matching
  !! that of a single simcomp, which executes the corresponding routines for
  !! each stored simcomp.
  !! The execution order is based on the order property of each simcomp.
  !! By default, the order is by the order of apparence in the case file.
  type, public :: simcomp_executor_t
     !> The simcomps.
     class(simulation_component_alloc_t), allocatable :: simcomps(:)
     !> Number of simcomps
     integer :: n_simcomps
     !> The case
     type(case_t), pointer :: case
   contains
     !> Constructor.
     procedure, pass(this) :: init => simcomp_executor_init
     !> Destructor.
     procedure, pass(this) :: free => simcomp_executor_free
     !> Appending a new simcomp to the executor.
     procedure, pass(this) :: add_user_simcomp => simcomp_executor_add
     !> Execute preprocess_ for all simcomps.
     procedure, pass(this) :: preprocess => simcomp_executor_preprocess
     !> Execute compute_ for all simcomps.
     procedure, pass(this) :: compute => simcomp_executor_compute
     !> Execute restart for all simcomps.
     procedure, pass(this) :: restart=> simcomp_executor_restart
     !> Finalize the initialization.
     procedure, pass(this) :: finalize => simcomp_executor_finalize
  end type simcomp_executor_t

  !> Global variable for the simulation component driver.
  type(simcomp_executor_t), target, public :: neko_simcomps

  !
  ! simulation_component_t procedures
  !

  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param json The JSON with properties.
     !! @param case The case_t object.
     subroutine simulation_component_init(this, json, case)
       import simulation_component_t, json_file, case_t
       class(simulation_component_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_init
  end interface

  abstract interface
     !> Destructor.
     subroutine simulation_component_free(this)
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine simulation_component_free
  end interface

  interface
     !> Simulation component factory.
     !! Both constructs and initializes the object.
     !! @param object The object to be created and initialized.
     !! @param json JSON object initializing the simulation component.
     !! @param case The simulation case.
     module subroutine simulation_component_factory(object, json, case)
       class(simulation_component_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_factory
  end interface

  interface
     module subroutine communicate(name)
       character(len=20) :: name
     end subroutine communicate
  end interface

  public :: simulation_component_factory, communicate

contains

  !
  ! simulation_component_t procedures
  !

  !> Constructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_init_base(this, json, case)
    class(simulation_component_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: preprocess_control, compute_control, &
         output_control
    real(kind=rp) :: preprocess_value, compute_value, output_value
    integer :: order
    character(len=:), allocatable :: name

    ! We default to preprocess every time-step
    call json_get_or_default(json, "preprocess_control", preprocess_control, &
                             "tsteps")
    call json_get_or_default(json, "preprocess_value", preprocess_value, 1.0_rp)

    ! We default to compute every time-step
    call json_get_or_default(json, "compute_control", compute_control, &
                             "tsteps")
    call json_get_or_default(json, "compute_value", compute_value, 1.0_rp)

    ! We default to output whenever we execute
    call json_get_or_default(json, "output_control", output_control, &
                             compute_control)
    call json_get_or_default(json, "output_value", output_value, &
                             compute_value)


    if (output_control == "global") then
       call json_get(case%params, 'case.fluid.output_control', output_control)
       call json_get(case%params, 'case.fluid.output_value', output_value)
    end if

    call json_get_or_default(json, "order", order, -1)

    call json_get_or_default(json, "name", name, "simcomp")

    write(*,*) "Hi from simcomp json"
    call this%init_base_from_components(case, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, order, name)

  end subroutine simulation_component_init_base

  subroutine simulation_component_init_base_from_components(this, case, &
       preprocess_control, preprocess_value, compute_control, compute_value, &
       output_control, output_value, order, name)
    class(simulation_component_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    integer, intent(in) :: order
    character(len=*), intent(inout) :: name

    write(*,*) "Hi"
    this%case => case
    this%order = order
    this%executor_ => neko_simcomps
    write(*,*) "H2"

    call this%preprocess_controller%init(case%end_time, preprocess_control, &
         preprocess_value)
    call this%compute_controller%init(case%end_time, compute_control, &
         compute_value)
    call this%output_controller%init(case%end_time, output_control, &
         output_value)

    write(*,*) "Getting name"
    call make_unique_name(this%executor_, name)
    this%name = name
    write(*,*) "My name is ", name

  end subroutine simulation_component_init_base_from_components

  recursive subroutine make_unique_name(executor, name)
    class(simcomp_executor_t), intent(in) :: executor
    character(len=*), intent(inout) :: name
    integer :: i
    logical :: name_exists
    character(len=:), allocatable :: new_name

    ! Check if the name already exists in the name_list
    name_exists = .false.
    do i = 1, executor%n_simcomps
       if (trim(name) == trim(executor%simcomps(i)%simcomp%name)) then
           name_exists = .true.
          exit
       end if
    end do

    write(*,*) name, name_exists

    if (.not. name_exists) return  ! Name is unique, return

    ! If the name exists, add an underscore and recursively check again
    new_name = trim(name) // '_'
    call make_unique_name(executor, new_name)
  end subroutine make_unique_name

  !> Destructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_free_base(this)
    class(simulation_component_t), intent(inout) :: this

    nullify(this%case)
    nullify(this%executor_)
  end subroutine simulation_component_free_base

  !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
  !! Serves as the public interface.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine simulation_component_preprocess_wrapper(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (this%preprocess_controller%check(t, tstep)) then
       call this%preprocess_(t, tstep)
       call this%preprocess_controller%register_execution()
    end if
  end subroutine simulation_component_preprocess_wrapper

  !> Wrapper for calling `compute_` based on the `compute_controller`.
  !! Serves as the public interface.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine simulation_component_compute_wrapper(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (this%compute_controller%check(t, tstep)) then
       call this%compute_(t, tstep)
       call this%compute_controller%register_execution()
    end if
  end subroutine simulation_component_compute_wrapper

  !> Wrapper for calling `set_counter_` based for the controllers.
  !! Serves as the public interface.
  !! @param t The time value.
  subroutine simulation_component_restart_wrapper(this, t)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%compute_controller%set_counter(t)
    call this%output_controller%set_counter(t)
    call this%restart_(t) 

  end subroutine simulation_component_restart_wrapper

  !> Dummy restart function.
  !! @param t The time value.
  subroutine restart_(this, t)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    ! Do nothing
  end subroutine restart_

  !> Dummy preprocessing function.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine preprocess_(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Do nothing
  end subroutine preprocess_

  !> Dummy compute function.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine compute_(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Do nothing
  end subroutine compute_

  !
  ! simcomp_exectutor_t procedures
  !

  !> Constructor.
  !! @param case The case.
  !! @param simcomp_root The root name of the simulation components in the case.
  !! If not provided, the default is 'case.simulation_components'.
  subroutine simcomp_executor_init(this, case, simcomp_root)
    class(simcomp_executor_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: case
    character(len=*), optional, intent(in) :: simcomp_root
    integer :: n_simcomps, i
    type(json_core) :: core
    type(json_value), pointer :: simcomp_object
    type(json_file) :: comp_subdict
    logical :: found, is_user, has_user
    ! Help array for finding minimal values
    logical, allocatable :: mask(:)
    ! The order value for each simcomp in order of appearance in the case file.
    integer, allocatable :: read_order(:), order(:)
    ! Location of the min value
    integer :: loc(1)
    integer :: max_order
    character(len=:), allocatable :: root_name, comp_type

    call this%free()
    this%case => case

    ! Get the root name of the simulation components if specified
    if (present(simcomp_root)) then
       root_name = simcomp_root
    else
       root_name = 'case.simulation_components'
    end if

    ! Get the core json object and the simulation components object
    call case%params%get_core(core)
    call case%params%get(root_name, simcomp_object, found)
    if (.not. found) return
    call neko_log%section('Initialize simcomp')

    ! Set the number of simcomps and allocate the arrays
    call case%params%info(root_name, n_children = n_simcomps)
    this%n_simcomps = n_simcomps
    allocate(this%simcomps(n_simcomps))
    allocate(order(n_simcomps))
    allocate(read_order(n_simcomps))
    allocate(mask(n_simcomps), source = .true.)

    ! We need a separate loop to figure out the order, so that we can
    ! apply the order to the initialization as well.
    max_order = 0
    has_user = .false.
    do i = 1, n_simcomps
       ! Create a new json containing just the subdict for this simcomp
       call json_extract_item(core, simcomp_object, i, comp_subdict)

       call json_get_or_default(comp_subdict, "is_user", is_user, .false.)
       has_user = has_user .or. is_user

       call json_get_or_default(comp_subdict, "order", read_order(i), -1)
       if (read_order(i) .gt. max_order) then
          max_order = read_order(i)
       end if
    end do

    ! If the order was not specified, we use the order of appearance in the
    ! case file.
    do i = 1, n_simcomps
       if (read_order(i) == -1) then
          max_order = max_order + 1
          read_order(i) = max_order
       end if
    end do

    ! Figure out the execution order using a poor man's argsort.
    ! Searches for the location of the min value, each time masking out the
    ! found location prior to the next search.
    do i = 1, n_simcomps
       loc = minloc(read_order, mask = mask)
       order(i) = loc(1)
       mask(loc) = .false.
    end do

    ! Init in the determined order.
    do i = 1, n_simcomps
       call json_extract_item(core, simcomp_object, order(i), comp_subdict)

       ! Log the component type if it is not a user component
       call json_get(comp_subdict, "type", comp_type)
       call json_get_or_default(comp_subdict, "is_user", is_user, .false.)
       if (.not. is_user) call neko_log%message('- ' // trim(comp_type))

       call simulation_component_factory(this%simcomps(i)%simcomp, &
            comp_subdict, case)
    end do

    if (has_user) then
       call neko_log%message('Initialize user simcomp')

       comp_subdict = json_file(simcomp_object)
       call case%usr%init_user_simcomp(comp_subdict)
    end if

    ! Cleanup
    call neko_simcomps%finalize()
    deallocate(order)
    deallocate(read_order)
    deallocate(mask)

    call neko_log%end_section()
  end subroutine simcomp_executor_init

  !> Destructor.
  subroutine simcomp_executor_free(this)
    class(simcomp_executor_t), intent(inout) :: this
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%free
       end do
       deallocate(this%simcomps)
    end if
  end subroutine simcomp_executor_free

  !> Appending a new simcomp to the executor.
  !! @param new_object The simcomp to append.
  !! @param settings The settings for the simcomp.
  subroutine simcomp_executor_add(this, object, settings)
    class(simcomp_executor_t), intent(inout) :: this
    class(simulation_component_t), intent(in) :: object
    type(json_file), intent(inout) :: settings

    class(simulation_component_alloc_t), allocatable :: tmp_simcomps(:)
    integer :: i, position

    ! Find the first empty position
    position = 0
    do i = 1, this%n_simcomps
       if (.not. allocated(this%simcomps(i)%simcomp)) then
          position = i
          exit
       end if
    end do

    ! If no empty position was found, append to the end
    if (position == 0) then
       call move_alloc(this%simcomps, tmp_simcomps)
       allocate(this%simcomps(this%n_simcomps + 1))

       if (allocated(tmp_simcomps)) then
          do i = 1, this%n_simcomps
             call move_alloc(tmp_simcomps(i)%simcomp, this%simcomps(i)%simcomp)
          end do
          deallocate(tmp_simcomps)
       end if

       this%n_simcomps = this%n_simcomps + 1
       position = this%n_simcomps
    end if

    this%simcomps(position)%simcomp = object
    call this%simcomps(position)%simcomp%init(settings, this%case)

    if (allocated(tmp_simcomps)) then
       deallocate(tmp_simcomps)
    end if

  end subroutine simcomp_executor_add

  !> Finalize the initialization.
  !! Sorts the simcomps based on the order property.
  !! Additionally we check that the order is unique, contiguous, starts at 1 and
  !! within bounds.
  subroutine simcomp_executor_finalize(this)
    class(simcomp_executor_t), intent(inout) :: this
    integer :: i, order, max_order
    logical :: order_found, previous_found

    class(simulation_component_alloc_t), allocatable :: tmp_simcomps(:)
    integer, allocatable :: order_list(:)

    ! Check that all components are initialized
    do i = 1, this%n_simcomps
       if (.not. allocated(this%simcomps(i)%simcomp)) then
          call neko_error("Simulation component not initialized.")
       end if
    end do

    ! Check that the order is unique and contiguous
    previous_found = .true.
    do order = 1, this%n_simcomps
       order_found = .false.
       do i = 1, this%n_simcomps
          if (this%simcomps(i)%simcomp%order == order .and. order_found) then
             call neko_error("Simulation component order must be unique.")
          else if (this%simcomps(i)%simcomp%order == order) then
             order_found = .true.
          end if
       end do
       if (order_found .and. .not. previous_found) then
          call neko_error("Simulation component order must be contiguous &
               &starting at 1.")
       end if
       previous_found = order_found
    end do

    allocate(order_list(this%n_simcomps))
    order_list = 0
    max_order = 0
    do i = 1, this%n_simcomps
       order_list(i) = this%simcomps(i)%simcomp%order
       if (order_list(i) .gt. max_order) then
          max_order = order_list(i)
       end if
    end do

    do i = 1, this%n_simcomps
       if (order_list(i) .eq. -1) then
          order_list(i) = max_order + 1
          max_order = max_order + 1
       end if
    end do

    ! Check that the order is within bounds
    do i = 1, this%n_simcomps
       if (order_list(i) .gt. this%n_simcomps) then
          deallocate(order_list)
          call neko_error("Simulation component order is out of bounds.")
       end if
    end do

    ! Reorder the simcomps based on the order specified
    call move_alloc(this%simcomps, tmp_simcomps)
    allocate(this%simcomps(this%n_simcomps))
    do i = 1, this%n_simcomps
       order = order_list(i)
       call move_alloc(tmp_simcomps(i)%simcomp, this%simcomps(order)%simcomp)
    end do

    if (allocated(tmp_simcomps)) then
       deallocate(tmp_simcomps)
    end if
    if (allocated(order_list)) then
       deallocate(order_list)
    end if

  end subroutine simcomp_executor_finalize
  !> Execute preprocess_ for all simcomps.
  !! @param t The time value.
  !! @param tstep The timestep number.
  subroutine simcomp_executor_preprocess(this, t, tstep)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, size(this%simcomps)
          call this%simcomps(i)%simcomp%preprocess(t, tstep)
       end do
    end if

  end subroutine simcomp_executor_preprocess

  !> Execute compute_ for all simcomps.
  !! @param t The time value.
  !! @param tstep The timestep number.
  subroutine simcomp_executor_compute(this, t, tstep)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%compute(t, tstep)
       end do
    end if

  end subroutine simcomp_executor_compute

  !> Execute restart for all simcomps.
  !! @param t The time value.
  subroutine simcomp_executor_restart(this, t)
    class(simcomp_executor_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (allocated(this%simcomps)) then
       do i = 1, this%n_simcomps
          call this%simcomps(i)%simcomp%restart(t)
       end do
    end if

  end subroutine simcomp_executor_restart
end module simulation_component
