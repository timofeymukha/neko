module user
  use neko
  implicit none

  type(file_t) :: error_file
  type(vector_t) :: linf_error
  character(len=:), allocatable :: error_filename

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%compute => compute
    user%finalize => finalize
    user%source_term => source_term
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    if (allocated(error_filename)) then
      deallocate(error_filename)
    end if

    call json_get_or_default(params, "case.error_file", &
         error_filename, "error.csv")
  end subroutine startup

  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: rhs_u, rhs_v
    real(kind=rp) :: x, y
    integer :: i

    if (scheme_name .ne. "fluid") return

    rhs_u => rhs%get_by_index(1)
    rhs_v => rhs%get_by_index(2)

    do i = 1, rhs_u%size()
      x = rhs_u%dof%x(i,1,1,1)
      y = rhs_u%dof%y(i,1,1,1)

      rhs_u%x(i,1,1,1) = -pi * sin(time%t) * sin(pi * x) * sin(pi * y) - &
           2.0_rp * pi**3 * cos(pi * x)**2 * sin(time%t) * sin(2.0_rp * pi * y) + &
           pi * cos(time%t) * sin(pi * x)**2 * sin(2.0_rp * pi * y) + &
           6.0_rp * pi**3 * sin(time%t) * sin(pi * x)**2 * sin(2.0_rp * pi * y)

      rhs_v%x(i,1,1,1) = pi * cos(pi * x) * cos(pi * y) * sin(time%t) + &
           2.0_rp * pi**3 * cos(pi * y)**2 * sin(time%t) * sin(2.0_rp * pi * x) - &
           pi * cos(time%t) * sin(pi * y)**2 * sin(2.0_rp * pi * x) - &
           6.0_rp * pi**3 * sin(time%t) * sin(pi * y)**2 * sin(2.0_rp * pi * x)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(rhs_u%x, rhs_u%x_d, rhs_u%size(), &
           HOST_TO_DEVICE, sync=.false.)
      call device_memcpy(rhs_v%x, rhs_v%x_d, rhs_v%size(), &
           HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine source_term

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    u => neko_registry%get_field("u")

    call neko_registry%add_field(u%dof, "error_u", ignore_existing=.true.)
    call neko_registry%add_field(u%dof, "error_v", ignore_existing=.true.)
    call neko_registry%add_field(u%dof, "error_p", ignore_existing=.true.)

    call error_file%init(trim(error_filename))
    call error_file%set_header("# t,error_u,error_v,error_p")
    call linf_error%init(3)

    call compute(time)
  end subroutine initialize

  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u, v, p
    type(field_t), pointer :: error_u, error_v, error_p
    real(kind=rp) :: x, y
    real(kind=rp) :: true_u, true_v, true_p
    integer :: i

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    p => neko_registry%get_field("p")

    error_u => neko_registry%get_field("error_u")
    error_v => neko_registry%get_field("error_v")
    error_p => neko_registry%get_field("error_p")

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call u%copy_from(DEVICE_TO_HOST, sync=.false.)
      call v%copy_from(DEVICE_TO_HOST, sync=.false.)
      call p%copy_from(DEVICE_TO_HOST, sync=.true.)
    end if

    do i = 1, u%size()
      x = u%dof%x(i,1,1,1)
      y = u%dof%y(i,1,1,1)

      call true_solution(true_u, true_v, true_p, x, y, time%t)

      error_u%x(i,1,1,1) = abs(u%x(i,1,1,1) - true_u)
      error_v%x(i,1,1,1) = abs(v%x(i,1,1,1) - true_v)
      error_p%x(i,1,1,1) = abs(p%x(i,1,1,1) - true_p)
    end do

    linf_error%x(1) = glmax(error_u%x, u%size())
    linf_error%x(2) = glmax(error_v%x, u%size())
    linf_error%x(3) = glmax(error_p%x, u%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(error_u%x, error_u%x_d, error_u%size(), &
           HOST_TO_DEVICE, sync=.false.)
      call device_memcpy(error_v%x, error_v%x_d, error_v%size(), &
           HOST_TO_DEVICE, sync=.false.)
      call device_memcpy(error_p%x, error_p%x_d, error_p%size(), &
           HOST_TO_DEVICE, sync=.false.)
    end if

    call error_file%write(linf_error, time%t)
  end subroutine compute

  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    call file_free(error_file)
    call linf_error%free()

    if (allocated(error_filename)) then
      deallocate(error_filename)
    end if
  end subroutine finalize

  subroutine true_solution(u, v, p, x, y, t)
    real(kind=rp), intent(inout) :: u, v, p
    real(kind=rp), intent(in) :: x, y, t

    u = pi * sin(t) * sin(2.0_rp * pi * y) * sin(pi * x)**2
    v = -pi * sin(t) * sin(2.0_rp * pi * x) * sin(pi * y)**2
    p = sin(t) * sin(pi * y) * cos(pi * x)
  end subroutine true_solution

end module user
