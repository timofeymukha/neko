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
!> Defines the valence-independent interface of a list of `bc_t`.
module bc_list
  use utils, only : neko_error
  use bc, only : bc_t
  implicit none
  private

  !> Abstract base for a list of boundary conditions.
  !! @details This type holds the queries that do not depend on whether the
  !! conditions in the list act on a scalar or on a vector field: iteration,
  !! lookup by name or zone, and the mask-emptiness check. Those are what the
  !! boundary-condition projectors and the marking code need, so they take a
  !! `class(bc_list_t)` and work with either concrete list.
  !!
  !! The concrete lists, `scalar_bc_list_t` and `vector_bc_list_t`, own the
  !! storage and add the corresponding `apply` interface. Keeping the storage
  !! in the concrete types is what makes the valence a compile-time property:
  !! a `vector_bc_t` cannot be appended to a scalar list at all, rather than
  !! failing at run time when the list is applied.
  type, public, abstract :: bc_list_t
   contains
     !> Get the item at the given index.
     procedure(bc_list_get_intrf), pass(this), deferred :: get
     !> Return the number of items in the list.
     procedure(bc_list_size_intrf), pass(this), deferred :: size
     !> Destructor.
     procedure(bc_list_free_intrf), pass(this), deferred :: free

     !> Get the item with the given name.
     procedure, pass(this) :: get_by_name => bc_list_get_by_name
     !> Get the item that applies to the given zone_index.
     procedure, pass(this) :: get_by_zone_index => bc_list_get_by_zone_index
     !> Check whether the list is empty
     procedure, pass(this) :: is_empty => bc_list_is_empty
     !> Return the type of a given item
     procedure, pass(this) :: bc_type => bc_list_bc_type
  end type bc_list_t

  abstract interface
     !> Get the item at the given index.
     !! @param i The index of the item to get.
     function bc_list_get_intrf(this, i) result(bc)
       import :: bc_list_t, bc_t
       class(bc_list_t), intent(in) :: this
       integer, intent(in) :: i
       class(bc_t), pointer :: bc
     end function bc_list_get_intrf
  end interface

  abstract interface
     !> Return the number of items in the list.
     pure function bc_list_size_intrf(this) result(size)
       import :: bc_list_t
       class(bc_list_t), intent(in) :: this
       integer :: size
     end function bc_list_size_intrf
  end interface

  abstract interface
     !> Destructor.
     subroutine bc_list_free_intrf(this)
       import :: bc_list_t
       class(bc_list_t), intent(inout) :: this
     end subroutine bc_list_free_intrf
  end interface

contains

  !> Get the item from a given name.
  !! @param name The name of the item to get.
  !! @return The item with the given name.
  function bc_list_get_by_name(this, name) result(bc)
    class(bc_list_t), intent(in) :: this
    class(bc_t), pointer :: bc
    character(len=*), intent(in) :: name
    integer :: i

    do i = 1, this%size()
       bc => this%get(i)
       if (bc%name .eq. trim(name)) return
    end do

    ! If the function reaches this point, no item was found
    call neko_error("Name not found in bc_list")

  end function bc_list_get_by_name

  !> Get the item from zone_index.
  !! @param zone_index where the bc applies.
  !! @return The item at the given zone_index.
  function bc_list_get_by_zone_index(this, zone_index) result(bc)
    class(bc_list_t), intent(in) :: this
    class(bc_t), pointer :: bc
    integer, intent(in) :: zone_index
    integer :: i, j

    do i = 1, this%size()
       bc => this%get(i)
       do j = 1, size(bc%zone_indices)
          if (bc%zone_indices(j) == zone_index) return
       end do
    end do

    ! If the function reaches this point, no item was found
    call neko_error("Zone index not found in bc_list")

  end function bc_list_get_by_zone_index

  !> Return the type of a given item.
  !! @param i The index of the item.
  function bc_list_bc_type(this, i) result(bc_type)
    class(bc_list_t), intent(in) :: this
    integer, intent(in) :: i
    integer :: bc_type
    class(bc_t), pointer :: bc

    bc => this%get(i)
    bc_type = bc%bc_type
  end function bc_list_bc_type

  !> Return whether the list is empty.
  function bc_list_is_empty(this) result(is_empty)
    class(bc_list_t), intent(in) :: this
    logical :: is_empty
    class(bc_t), pointer :: bc
    integer :: i

    is_empty = .true.
    do i = 1, this%size()
       bc => this%get(i)

       if (.not. allocated(bc%msk)) then
          call neko_error("bc not finalized, error in bc_list%is_empty")
       end if

       if (bc%msk(0) > 0) is_empty = .false.

    end do
  end function bc_list_is_empty

end module bc_list
