! Copyright (c) 2020-2026, The Neko Authors
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
!> Defines a boundary condition that constrains a vector field.
module vector_bc
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : glb_cmd_queue
  use bc, only : bc_t
  use field, only : field_t
  use time_state, only : time_state_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Base type for a boundary condition acting on a vector field.
  !! @details Conditions that constrain the three components of a field, such
  !! as the velocity, extend this type. It adds the vector apply interface on
  !! top of the masking and marking machinery of `bc_t`. Conditions that
  !! constrain a single field extend `scalar_bc_t` instead, and conditions that
  !! only need a boundary mask extend `bc_t` directly.
  type, public, abstract, extends(bc_t) :: vector_bc_t
   contains
     !> Apply the boundary condition to a vector field. Dispatches to the CPU
     !! or the device version.
     procedure, pass(this) :: apply_vector_generic => &
          vector_bc_apply_vector_generic
     !> Apply the boundary condition to a vector field on the CPU.
     procedure(vector_bc_apply_vector), pass(this), deferred :: apply_vector
     !> Device version of \ref apply_vector.
     procedure(vector_bc_apply_vector_dev), pass(this), deferred :: &
          apply_vector_dev
  end type vector_bc_t

  !> Pointer to a @ref `vector_bc_t`.
  type, public :: vector_bc_ptr_t
     class(vector_bc_t), pointer :: ptr => null()
  end type vector_bc_ptr_t

  abstract interface
     !> Apply the boundary condition to a vector field
     !! @param x The x comp of the field for which to apply the bc.
     !! @param y The y comp of the field for which to apply the bc.
     !! @param z The z comp of the field for which to apply the bc.
     !! @param n The size of x, y, and z.
     !! @param time Current time state.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine vector_bc_apply_vector(this, x, y, z, n, time, strong)
       import :: vector_bc_t, time_state_t
       import :: rp
       class(vector_bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout), dimension(n) :: x
       real(kind=rp), intent(inout), dimension(n) :: y
       real(kind=rp), intent(inout), dimension(n) :: z
       type(time_state_t), intent(in), optional :: time
       logical, intent(in), optional :: strong
     end subroutine vector_bc_apply_vector
  end interface

  abstract interface
     !> Apply the boundary condition to a vector field on the device.
     !! @param x_d Device pointer to the values to be applied for the x comp.
     !! @param y_d Device pointer to the values to be applied for the y comp.
     !! @param z_d Device pointer to the values to be applied for the z comp.
     !! @param time The time state.
     !! @param strong Whether we are setting a strong or a weak bc.
     !! @param strm Device stream
     subroutine vector_bc_apply_vector_dev(this, x_d, y_d, z_d, time, strong, &
          strm)
       import :: c_ptr, vector_bc_t, time_state_t
       import :: rp
       class(vector_bc_t), intent(inout), target :: this
       type(c_ptr), intent(inout) :: x_d
       type(c_ptr), intent(inout) :: y_d
       type(c_ptr), intent(inout) :: z_d
       type(time_state_t), intent(in), optional :: time
       logical, intent(in), optional :: strong
       type(c_ptr), intent(inout) :: strm
     end subroutine vector_bc_apply_vector_dev
  end interface

contains

  !> Apply the boundary condition to a vector field. Dispatches to the CPU
  !! or the device version.
  !! @param x The x comp of the field for which to apply the bc.
  !! @param y The y comp of the field for which to apply the bc.
  !! @param z The z comp of the field for which to apply the bc.
  !! @param time Current time state.
  !! @param strong Whether we are setting a strong or a weak bc.
  !! @param strm Device stream
  subroutine vector_bc_apply_vector_generic(this, x, y, z, time, strong, strm)
    class(vector_bc_t), intent(inout) :: this
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout), optional :: strm
    type(c_ptr) :: strm_
    integer :: n
    character(len=256) :: msg

    ! Get the size of the fields
    n = x%size()

    ! Ensure all fields are the same size
    if (y%size() .ne. n .or. z%size() .ne. n) then
       msg = "Fields x, y, z must have the same size in " // &
            "vector_bc_apply_vector_generic"
       call neko_error(trim(msg))
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then

       if (present(strm)) then
          strm_ = strm
       else
          strm_ = glb_cmd_queue
       end if

       call this%apply_vector_dev(x%x_d, y%x_d, z%x_d, time = time, &
            strong = strong, strm = strm_)
    else
       call this%apply_vector(x%x, y%x, z%x, n, time = time, strong = strong)
    end if

  end subroutine vector_bc_apply_vector_generic

end module vector_bc
