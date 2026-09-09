! Copyright (c) 2024-2025, The Neko Authors
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
!> Defines a list of `vector_bc_t`.
module vector_bc_list
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use field, only : field_t
  use device, only : device_get_ptr, glb_cmd_queue
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use vector_bc, only : vector_bc_t, vector_bc_ptr_t
  use time_state, only : time_state_t
  !$ use omp_lib
  implicit none
  private

  !> A list of boundary conditions acting on a vector field.
  !! Follows the standard interface of lists.
  type, public, extends(bc_list_t) :: vector_bc_list_t
     ! The items of the list.
     type(vector_bc_ptr_t), allocatable, private :: items(:)
     !> Number of items in the list that are themselves allocated.
     integer, private :: size_ = 0
     !> Capacity, i.e. the size of the items list. Some items may themselves be
     !! unallocated.
     integer, private :: capacity = 0
   contains
     !> Constructor.
     procedure, pass(this) :: init => vector_bc_list_init
     !> Destructor.
     procedure, pass(this) :: free => vector_bc_list_free

     !> Append an item to the end of the list.
     procedure, pass(this) :: append => vector_bc_list_append
     !> Get the item at the given index.
     procedure, pass(this) :: get => vector_bc_list_get
     !> Get the item at the given index, with its valence.
     procedure, pass(this) :: get_vector => vector_bc_list_get_vector
     !> Return the number of items in the list.
     procedure, pass(this) :: size => vector_bc_list_size

     !> Apply all boundary conditions in the list.
     generic :: apply => apply_vector, apply_vector_device, apply_vector_field
     !> Apply the boundary conditions to a vector array.
     procedure, pass(this) :: apply_vector => vector_bc_list_apply_array
     !> Apply the boundary conditions to a vector device array.
     procedure, pass(this) :: apply_vector_device => &
          vector_bc_list_apply_device
     !> Apply the boundary conditions to a vector field.
     procedure, pass(this) :: apply_vector_field => vector_bc_list_apply_field
  end type vector_bc_list_t

contains

  !> Constructor.
  !! @param capacity The size of the list to allocate.
  subroutine vector_bc_list_init(this, capacity)
    class(vector_bc_list_t), intent(inout), target :: this
    integer, optional :: capacity
    integer :: n

    call this%free()

    n = 1
    if (present(capacity)) n = capacity

    allocate(this%items(n))

    this%size_ = 0
    this%capacity = n

  end subroutine vector_bc_list_init

  !> Destructor.
  !! @note This will only nullify all pointers, not deallocate any
  !! conditions pointed to by the list
  subroutine vector_bc_list_free(this)
    class(vector_bc_list_t), intent(inout) :: this
    integer :: i

    if (allocated(this%items)) then
       do i = 1, this%size_
          this%items(i)%ptr => null()
       end do

       deallocate(this%items)
    end if

    this%size_ = 0
    this%capacity = 0
  end subroutine vector_bc_list_free

  !> Append a condition to the end of the list.
  !! @param bc The boundary condition to add.
  !! @details Will add the object to the list, even if the mask has zero size.
  subroutine vector_bc_list_append(this, bc)
    class(vector_bc_list_t), intent(inout) :: this
    class(vector_bc_t), intent(inout), target :: bc
    type(vector_bc_ptr_t), allocatable :: tmp(:)

    if (this%size_ .ge. this%capacity) then
       this%capacity = max(this%capacity * 2, 1)
       allocate(tmp(this%capacity))
       tmp(1:this%size_) = this%items(1:this%size_)
       call move_alloc(tmp, this%items)
    end if

    this%size_ = this%size_ + 1
    this%items(this%size_)%ptr => bc

  end subroutine vector_bc_list_append

  !> Get the item at the given index.
  !! @param i The index of the item to get.
  !! @return The item at the given index.
  function vector_bc_list_get(this, i) result(bc)
    class(vector_bc_list_t), intent(in) :: this
    class(bc_t), pointer :: bc
    integer, intent(in) :: i

    if (i .lt. 1 .or. i .gt. this%size_) then
       call neko_error("Index out of bounds in vector_bc_list_get")
    end if

    bc => this%items(i)%ptr

  end function vector_bc_list_get

  !> Get the item at the given index, keeping its valence in the declared
  !! type. `get` returns a `bc_t` pointer, which is what the mask-level code
  !! wants; use this where the caller needs to apply the condition or pass it
  !! on to something that does.
  !! @param i The index of the item to get.
  !! @return The item at the given index.
  function vector_bc_list_get_vector(this, i) result(bc)
    class(vector_bc_list_t), intent(in) :: this
    class(vector_bc_t), pointer :: bc
    integer, intent(in) :: i

    if (i .lt. 1 .or. i .gt. this%size_) then
       call neko_error("Index out of bounds in vector_bc_list_get_vector")
    end if

    bc => this%items(i)%ptr

  end function vector_bc_list_get_vector

  !> Return the number of items in the list.
  pure function vector_bc_list_size(this) result(size)
    class(vector_bc_list_t), intent(in) :: this
    integer :: size

    size = this%size_
  end function vector_bc_list_size

  !> Apply a list of boundary conditions to a vector array.
  !! @param x The x comp of the field for which to apply the bcs.
  !! @param y The y comp of the field for which to apply the bcs.
  !! @param z The z comp of the field for which to apply the bcs.
  !! @param n The size of x, y, z.
  !! @param time Current time state.
  !! @param strong Filter for strong or weak boundary conditions. Default is to
  !! apply the whole list.
  !! @param strm Device stream
  subroutine vector_bc_list_apply_array(this, x, y, z, n, time, strong, strm)
    class(vector_bc_list_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout), optional :: strm
    logical :: strong_
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then

       x_d = device_get_ptr(x)
       y_d = device_get_ptr(y)
       z_d = device_get_ptr(z)

       call this%apply_vector_device(x_d, y_d, z_d, time = time, &
            strong = strong, strm = strm)
    else
       ! Resolve strong into a concrete, always-present local before opening
       ! the parallel region. CCE's outlined region prologue dereferences a
       ! null descriptor for any absent optional dummy referenced inside the
       ! region, and it does not treat an unallocated allocatable as not
       ! present, so optionals must not be forwarded into the region. time
       ! cannot be given a meaningful concrete default, so we branch on
       ! present(time) outside the region instead.
       strong_ = .true.
       if (present(strong)) strong_ = strong

       if (present(time)) then
          !$omp parallel
          do i = 1, this%size_
             call this%items(i)%ptr%apply_vector(x, y, z, n, time = time, &
                  strong = strong_)
          end do
          !$omp end parallel
       else
          !$omp parallel
          do i = 1, this%size_
             call this%items(i)%ptr%apply_vector(x, y, z, n, strong = strong_)
          end do
          !$omp end parallel
       end if
    end if

  end subroutine vector_bc_list_apply_array

  !> Apply a list of boundary conditions to a vector array on the device.
  !! @param x_d The x comp of the field for which to apply the bcs.
  !! @param y_d The y comp of the field for which to apply the bcs.
  !! @param z_d The z comp of the field for which to apply the bcs.
  !! @param time Current time state.
  !! @param strong Filter for strong or weak boundary conditions. Default is to
  !! apply the whole list.
  !! @param strm Device stream
  subroutine vector_bc_list_apply_device(this, x_d, y_d, z_d, time, strong, &
       strm)
    class(vector_bc_list_t), intent(inout) :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout), optional :: strm
    type(c_ptr) :: strm_
    integer :: i

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    do i = 1, this%size_
       call this%items(i)%ptr%apply_vector_dev(x_d, y_d, z_d, time = time, &
            strong = strong, strm = strm_)
    end do

  end subroutine vector_bc_list_apply_device

  !> Apply a list of boundary conditions to a vector field.
  !! @param x The x comp of the field for which to apply the bcs.
  !! @param y The y comp of the field for which to apply the bcs.
  !! @param z The z comp of the field for which to apply the bcs.
  !! @param time Current time state.
  !! @param strong Filter for strong or weak boundary conditions. Default is to
  !! apply the whole list.
  !! @param strm Device stream
  subroutine vector_bc_list_apply_field(this, x, y, z, time, strong, strm)
    class(vector_bc_list_t), intent(inout) :: this
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout), optional :: strm
    logical :: strong_
    type(c_ptr) :: strm_
    integer :: i

    ! Resolve strong and strm into concrete, always-present locals before
    ! opening the parallel region. CCE's outlined region prologue dereferences
    ! a null descriptor for any absent optional dummy referenced inside the
    ! region, and it does not treat an unallocated allocatable as not present,
    ! so optionals must not be forwarded into the region. time cannot be given
    ! a meaningful concrete default, so we branch on present(time) outside the
    ! region instead.
    strong_ = .true.
    if (present(strong)) strong_ = strong
    strm_ = glb_cmd_queue
    if (present(strm)) strm_ = strm

    if (present(time)) then
       !$omp parallel if (.not. omp_in_parallel())
       do i = 1, this%size_
          call this%items(i)%ptr%apply_vector_generic(x, y, z, time = time, &
               strong = strong_, strm = strm_)
       end do
       !$omp end parallel
    else
       !$omp parallel if (.not. omp_in_parallel())
       do i = 1, this%size_
          call this%items(i)%ptr%apply_vector_generic(x, y, z, &
               strong = strong_, strm = strm_)
       end do
       !$omp end parallel
    end if

  end subroutine vector_bc_list_apply_field

end module vector_bc_list
