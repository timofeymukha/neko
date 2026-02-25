! Copyright (c) 2024, The Neko Authors
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
!> Implements `spalart_allmaras_t`.
module spalart_allmaras
  use num_types, only : rp
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use les_model, only : les_model_t
  use json_utils, only : json_get_or_default
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the spalart_allmaras LES model.
  !! @note Reference DOI: 10.1175/1520-0493(1963)091<0099:GCEWTP>2.3.CO;2
  type, public, extends(les_model_t) :: spalart_allmaras_t
     !> Model constant, defaults to 0.17.
     real(kind=rp) :: c_s
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => spalart_allmaras_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          spalart_allmaras_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => spalart_allmaras_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => spalart_allmaras_compute
  end type spalart_allmaras_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine spalart_allmaras_init(this, fluid, json)
    class(spalart_allmaras_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name
    real(kind=rp) :: c_s
    character(len=:), allocatable :: delta_type
    logical :: if_ext
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "c_s", c_s, 0.17_rp)
    call json_get_or_default(json, "extrapolation", if_ext, .false.)

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : spalart_allmaras'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c_s : ', c_s
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call spalart_allmaras_init_from_components(this, fluid, c_s, nut_name, &
         delta_type, if_ext)

  end subroutine spalart_allmaras_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c_s The model constant.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether trapolate the velocity.
  subroutine spalart_allmaras_init_from_components(this, fluid, c_s, &
       nut_name, delta_type, if_ext)
    class(spalart_allmaras_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp) :: c_s
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    call this%free()

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%c_s = c_s

  end subroutine spalart_allmaras_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine spalart_allmaras_free(this)
    class(spalart_allmaras_t), intent(inout) :: this

    call this%free_base()
  end subroutine spalart_allmaras_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine spalart_allmaras_compute(this, t, tstep)
    class(spalart_allmaras_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    type(field_t), pointer :: u, v, w, u_e, v_e, w_e

    if (this%if_ext .eqv. .true.) then
       ! Extrapolate the velocity fields
       associate(ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag, ext_bdf => this%ext_bdf)

         u => neko_registry%get_field_by_name("u")
         v => neko_registry%get_field_by_name("v")
         w => neko_registry%get_field_by_name("w")
         u_e => neko_registry%get_field_by_name("u_e")
         v_e => neko_registry%get_field_by_name("v_e")
         w_e => neko_registry%get_field_by_name("w_e")

         call this%sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
              ulag, vlag, wlag, ext_bdf%advection_coeffs, ext_bdf%nadv)

       end associate
    end if

  end subroutine spalart_allmaras_compute

end module spalart_allmaras
