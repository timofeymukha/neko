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
!> Defines a boundary condition that constrains a scalar field.
module scalar_bc
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : glb_cmd_queue
  use bc, only : bc_t
  use field, only : field_t
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Base type for a boundary condition acting on a scalar field.
  !! @details Conditions that constrain a single field, such as the pressure or
  !! the temperature, extend this type. It adds the scalar apply interface on
  !! top of the masking and marking machinery of `bc_t`. Conditions that
  !! constrain a vector field extend `vector_bc_t` instead, and conditions that
  !! only need a boundary mask extend `bc_t` directly.
  type, public, abstract, extends(bc_t) :: scalar_bc_t
   contains
     !> Apply the boundary condition to a scalar field. Dispatches to the CPU
     !! or the device version.
     procedure, pass(this) :: apply_scalar_generic => &
          scalar_bc_apply_scalar_generic
     !> Apply the boundary condition to a scalar field on the CPU.
     procedure(scalar_bc_apply_scalar), pass(this), deferred :: apply_scalar
     !> Device version of \ref apply_scalar.
     procedure(scalar_bc_apply_scalar_dev), pass(this), deferred :: &
          apply_scalar_dev
  end type scalar_bc_t

  !> Pointer to a @ref `scalar_bc_t`.
  type, public :: scalar_bc_ptr_t
     class(scalar_bc_t), pointer :: ptr => null()
  end type scalar_bc_ptr_t

  abstract interface
     !> Apply the boundary condition to a scalar field
     !! @param x The field for which to apply the boundary condition.
     !! @param n The size of x.
     !! @param time Current time state.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine scalar_bc_apply_scalar(this, x, n, time, strong)
       import :: scalar_bc_t, time_state_t
       import :: rp
       class(scalar_bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout), dimension(n) :: x
       type(time_state_t), intent(in), optional :: time
       logical, intent(in), optional :: strong
     end subroutine scalar_bc_apply_scalar
  end interface

  abstract interface
     !> Apply the boundary condition to a scalar field on the device
     !! @param x_d Device pointer to the field.
     !! @param time The time state.
     !! @param strong Whether we are setting a strong or a weak bc.
     !! @param strm Device stream
     subroutine scalar_bc_apply_scalar_dev(this, x_d, time, strong, strm)
       import :: c_ptr
       import :: scalar_bc_t, time_state_t
       import :: rp
       class(scalar_bc_t), intent(inout), target :: this
       type(c_ptr), intent(inout) :: x_d
       type(time_state_t), intent(in), optional :: time
       logical, intent(in), optional :: strong
       type(c_ptr), intent(inout) :: strm
     end subroutine scalar_bc_apply_scalar_dev
  end interface

contains

  !> Apply the boundary condition to a scalar field. Dispatches to the CPU
  !! or the device version.
  !! @param x The field for which to apply the bc.
  !! @param time Current time state.
  !! @param strong Whether we are setting a strong or a weak bc.
  !! @param strm Device stream
  subroutine scalar_bc_apply_scalar_generic(this, x, time, strong, strm)
    class(scalar_bc_t), intent(inout) :: this
    type(field_t), intent(inout) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout), optional :: strm
    type(c_ptr) :: strm_
    integer :: n

    ! Get the size of the field
    n = x%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then

       if (present(strm)) then
          strm_ = strm
       else
          strm_ = glb_cmd_queue
       end if

       call this%apply_scalar_dev(x%x_d, time = time, strong = strong, &
            strm = strm_)
    else
       call this%apply_scalar(x%x, n, time = time, strong = strong)
    end if

  end subroutine scalar_bc_apply_scalar_generic

end module scalar_bc
