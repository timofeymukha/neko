! Copyright (c) 2024-2026, The Neko Authors
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
  use operators, only : dudxyz
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use logger, only : LOG_SIZE, neko_log
  use utils, only : neko_error
  use zero_dirichlet, only : zero_dirichlet_t
  use math, only : cfill, glimax
  use gs_ops, only : GS_OP_MIN
  implicit none
  private
  character(len=*), parameter :: WALL_DISTANCE_FIELD_NAME = "wall_distance"

  !> Implements the Spalart-Allmaras RANS model.
  !! @note The scalar equation itself is solved outside this model.
  type, public, extends(les_model_t) :: spalart_allmaras_t
     !> Name of the scalar field carrying the SA working variable.
     character(len=:), allocatable :: scalar_field_name
     !> Name of the field containing the combined SA source term.
     character(len=:), allocatable :: source_field_name
     !> Whether to copy the SA scalar field into `nut`.
     logical :: sync_nut_with_scalar = .true.
     !> Minimum wall distance used in source term evaluation.
     real(kind=rp) :: min_wall_distance = 1.0e-12_rp
     !> Model constant cb1.
     real(kind=rp) :: cb1 = 0.1355_rp
     !> Model constant cb2.
     real(kind=rp) :: cb2 = 0.622_rp
     !> Model constant sigma.
     real(kind=rp) :: sigma_sa = 2.0_rp / 3.0_rp
     !> Model constant cw1.
     real(kind=rp) :: cw1 = 3.23906781678_rp
     !> Model constant cw2.
     real(kind=rp) :: cw2 = 0.3_rp
     !> Model constant cw3.
     real(kind=rp) :: cw3 = 2.0_rp
     !> Model constant cv1.
     real(kind=rp) :: cv1 = 7.1_rp
     !> Model constant ct3.
     real(kind=rp) :: ct3 = 1.2_rp
     !> Model constant ct4.
     real(kind=rp) :: ct4 = 0.5_rp
     !> Model constant kappa.
     real(kind=rp) :: kappa = 0.41_rp
     !> Zone indices for cheap wall-distance evaluation.
     integer, allocatable :: wall_distance_zone_indices(:)
     !> Pointer to the combined source field.
     type(field_t), pointer :: source_field => null()
     !> Pointer to the dynamic viscosity field.
     type(field_t), pointer :: mu => null()
     !> Pointer to the density field.
     type(field_t), pointer :: rho => null()
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
     !> Compute wall distance based on wall zone indices.
     procedure, pass(this) :: compute_cheap_dist => &
          spalart_allmaras_compute_cheap_dist_cpu
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
    character(len=:), allocatable :: scalar_field_name
    character(len=:), allocatable :: source_field_name
    integer, allocatable :: wall_distance_zone_indices(:)
    character(len=:), allocatable :: delta_type
    logical :: sync_nut_with_scalar
    logical :: if_ext
    real(kind=rp) :: cb1, cb2, sigma_sa, cw1, min_wall_distance
    real(kind=rp) :: cw2, cw3, cv1, ct3, ct4, kappa
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "scalar_field", scalar_field_name, &
         trim(nut_name))
    call json_get_or_default(json, "source_field", source_field_name, &
         "source")
    if (json%valid_path("wall_distance_zone_indices")) then
       call json_get(json, "wall_distance_zone_indices", &
            wall_distance_zone_indices)
    else
       allocate(wall_distance_zone_indices(0))
    end if
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "sync_nut_with_scalar", &
         sync_nut_with_scalar, .true.)
    call json_get_or_default(json, "cb1", cb1, 0.1355_rp)
    call json_get_or_default(json, "cb2", cb2, 0.622_rp)
    call json_get_or_default(json, "sigma", sigma_sa, 2.0_rp / 3.0_rp)
    call json_get_or_default(json, "cw1", cw1, 3.23906781678_rp)
    call json_get_or_default(json, "cw2", cw2, 0.3_rp)
    call json_get_or_default(json, "cw3", cw3, 2.0_rp)
    call json_get_or_default(json, "cv1", cv1, 7.1_rp)
    call json_get_or_default(json, "ct3", ct3, 1.2_rp)
    call json_get_or_default(json, "ct4", ct4, 0.5_rp)
    call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "min_wall_distance", min_wall_distance, &
         1.0e-12_rp)
    call json_get_or_default(json, "extrapolation", if_ext, .false.)

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : Spalart-Allmaras'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Scalar field : ', scalar_field_name
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Source field : ', source_field_name
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Wall distance field : ', &
         WALL_DISTANCE_FIELD_NAME
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call spalart_allmaras_init_from_components(this, fluid, &
         nut_name, scalar_field_name, source_field_name, &
         wall_distance_zone_indices, delta_type, if_ext, &
         sync_nut_with_scalar, cb1, cb2, sigma_sa, cw1, cw2, cw3, cv1, &
         ct3, ct4, kappa, min_wall_distance)

  end subroutine spalart_allmaras_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param scalar_field_name The name of the scalar field.
  !! @param source_field_name The name of the combined source field.
  !! @param wall_distance_zone_indices Wall zones for cheap distance.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether to extrapolate the velocity.
  !! @param sync_nut_with_scalar Whether to sync `nut` from the SA scalar.
  !! @param cb1 Model constant cb1.
  !! @param cb2 Model constant cb2.
  !! @param sigma_sa Model constant sigma.
  !! @param cw1 Model constant cw1.
  !! @param cw2 Model constant cw2.
  !! @param cw3 Model constant cw3.
  !! @param cv1 Model constant cv1.
  !! @param ct3 Model constant ct3.
  !! @param ct4 Model constant ct4.
  !! @param kappa Model constant kappa.
  !! @param min_wall_distance Minimum wall distance in source terms.
  subroutine spalart_allmaras_init_from_components(this, fluid, nut_name, &
       scalar_field_name, source_field_name, wall_distance_zone_indices, &
       delta_type, if_ext, sync_nut_with_scalar, &
       cb1, cb2, sigma_sa, cw1, cw2, cw3, cv1, ct3, ct4, kappa, &
       min_wall_distance)
    class(spalart_allmaras_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: scalar_field_name
    character(len=*), intent(in) :: source_field_name
    integer, intent(in) :: wall_distance_zone_indices(:)
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext
    logical, intent(in) :: sync_nut_with_scalar
    real(kind=rp), intent(in) :: cb1, cb2, sigma_sa, cw1, cw2, cw3
    real(kind=rp), intent(in) :: cv1, ct3, ct4, kappa, min_wall_distance
    type(field_t), pointer :: wall_distance
    real(kind=rp), pointer :: wall_distance_vec(:)

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Spalart-Allmaras model for device backend is " // &
            "not implemented yet")
    end if

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%scalar_field_name = trim(scalar_field_name)
    this%source_field_name = trim(source_field_name)
    this%sync_nut_with_scalar = sync_nut_with_scalar
    this%cb1 = cb1
    this%cb2 = cb2
    this%sigma_sa = sigma_sa
    this%cw1 = cw1
    this%cw2 = cw2
    this%cw3 = cw3
    this%cv1 = cv1
    this%ct3 = ct3
    this%ct4 = ct4
    this%kappa = kappa
    this%min_wall_distance = min_wall_distance
    allocate(this%wall_distance_zone_indices(size(wall_distance_zone_indices)))
    this%wall_distance_zone_indices = wall_distance_zone_indices
    this%mu => fluid%mu
    this%rho => fluid%rho

    call neko_registry%add_field(fluid%dm_Xh, trim(this%source_field_name), &
         ignore_existing = .true.)
    this%source_field => &
         neko_registry%get_field_by_name(trim(this%source_field_name))

    call neko_registry%add_field(fluid%dm_Xh, &
         WALL_DISTANCE_FIELD_NAME, ignore_existing = .true.)
    wall_distance => &
         neko_registry%get_field_by_name(WALL_DISTANCE_FIELD_NAME)
    if (size(this%wall_distance_zone_indices) .eq. 0) then
       call neko_error("Spalart-Allmaras requires " // &
            "wall_distance_zone_indices")
    end if
    wall_distance_vec(1:wall_distance%size()) => wall_distance%x
    call this%compute_cheap_dist(wall_distance_vec)

  end subroutine spalart_allmaras_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine spalart_allmaras_free(this)
    class(spalart_allmaras_t), intent(inout) :: this

    nullify(this%source_field)
    nullify(this%mu)
    nullify(this%rho)
    if (allocated(this%scalar_field_name)) then
       deallocate(this%scalar_field_name)
    end if
    if (allocated(this%source_field_name)) then
       deallocate(this%source_field_name)
    end if
    if (allocated(this%wall_distance_zone_indices)) then
       deallocate(this%wall_distance_zone_indices)
    end if

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
    type(field_t), pointer :: u_work, v_work, w_work
    type(field_t), pointer :: scalar, wall_distance
    type(field_t), pointer :: du_dx, du_dy, du_dz
    type(field_t), pointer :: dv_dx, dv_dy, dv_dz
    type(field_t), pointer :: dw_dx, dw_dy, dw_dz
    type(field_t), pointer :: dsa_dx, dsa_dy, dsa_dz
    integer :: idx(12)
    integer :: e, i, j, k
    real(kind=rp) :: omega_x, omega_y, omega_z, omega_mag
    real(kind=rp) :: grad_sq, d_eff, nu_tilde, nu, chi, chi3
    real(kind=rp) :: fv1, fv2, ft2, s_tilde, denom, r, g, fw
    real(kind=rp) :: cw3_6, g6

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

    u => neko_registry%get_field_by_name("u")
    v => neko_registry%get_field_by_name("v")
    w => neko_registry%get_field_by_name("w")

    if (this%if_ext) then
      u_work => neko_registry%get_field_by_name("u_e")
      v_work => neko_registry%get_field_by_name("v_e")
      w_work => neko_registry%get_field_by_name("w_e")
    else
      u_work => u
      v_work => v
      w_work => w
    end if

    scalar => neko_registry%get_field_by_name(trim(this%scalar_field_name))
    wall_distance => &
         neko_registry%get_field_by_name(WALL_DISTANCE_FIELD_NAME)
    if (.not. associated(this%mu) .or. .not. associated(this%rho)) then
       call neko_error("Spalart-Allmaras model requires associated " // &
            "fluid mu and rho fields")
    end if

    call neko_scratch_registry%request_field(du_dx, idx(1), .false.)
    call neko_scratch_registry%request_field(du_dy, idx(2), .false.)
    call neko_scratch_registry%request_field(du_dz, idx(3), .false.)
    call neko_scratch_registry%request_field(dv_dx, idx(4), .false.)
    call neko_scratch_registry%request_field(dv_dy, idx(5), .false.)
    call neko_scratch_registry%request_field(dv_dz, idx(6), .false.)
    call neko_scratch_registry%request_field(dw_dx, idx(7), .false.)
    call neko_scratch_registry%request_field(dw_dy, idx(8), .false.)
    call neko_scratch_registry%request_field(dw_dz, idx(9), .false.)
    call neko_scratch_registry%request_field(dsa_dx, idx(10), .false.)
    call neko_scratch_registry%request_field(dsa_dy, idx(11), .false.)
    call neko_scratch_registry%request_field(dsa_dz, idx(12), .false.)

    call dudxyz(du_dx%x, u_work%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(du_dy%x, u_work%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(du_dz%x, u_work%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)
    call dudxyz(dv_dx%x, v_work%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dv_dy%x, v_work%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(dv_dz%x, v_work%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)
    call dudxyz(dw_dx%x, w_work%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dw_dy%x, w_work%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(dw_dz%x, w_work%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)

    call dudxyz(dsa_dx%x, scalar%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dsa_dy%x, scalar%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(dsa_dz%x, scalar%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)

    cw3_6 = this%cw3**6

    do e = 1, this%coef%msh%nelv
       do k = 1, this%coef%Xh%lz
          do j = 1, this%coef%Xh%ly
             do i = 1, this%coef%Xh%lx
                omega_x = dw_dy%x(i, j, k, e) - dv_dz%x(i, j, k, e)
                omega_y = du_dz%x(i, j, k, e) - dw_dx%x(i, j, k, e)
                omega_z = dv_dx%x(i, j, k, e) - du_dy%x(i, j, k, e)
                omega_mag = sqrt(max(0.0_rp, &
                     omega_x * omega_x + omega_y * omega_y + omega_z * omega_z))

                grad_sq = dsa_dx%x(i, j, k, e) * dsa_dx%x(i, j, k, e) + &
                     dsa_dy%x(i, j, k, e) * dsa_dy%x(i, j, k, e) + &
                     dsa_dz%x(i, j, k, e) * dsa_dz%x(i, j, k, e)

                d_eff = max(abs(wall_distance%x(i, j, k, e)), &
                     this%min_wall_distance)
                nu_tilde = max(scalar%x(i, j, k, e), 0.0_rp)
                nu = max(this%mu%x(i, j, k, e) / &
                     max(this%rho%x(i, j, k, e), tiny(1.0_rp)), tiny(1.0_rp))

                chi = nu_tilde / nu
                chi3 = chi * chi * chi
                fv1 = chi3 / (chi3 + this%cv1**3)
                fv2 = 1.0_rp - chi / (1.0_rp + chi * fv1)
                ft2 = this%ct3 * exp(-this%ct4 * chi * chi)

                s_tilde = omega_mag + nu_tilde * fv2 / &
                     (this%kappa**2 * d_eff * d_eff)
                denom = max(s_tilde, 0.3_rp * omega_mag) * &
                     this%kappa**2 * d_eff * d_eff
                if (denom .gt. tiny(1.0_rp)) then
                   r = min(nu_tilde / denom, 10.0_rp)
                else
                   r = 10.0_rp
                end if

                g = r + this%cw2 * (r**6 - r)
                g6 = g**6
                fw = g * ((1.0_rp + cw3_6) / (g6 + cw3_6))**(1.0_rp / 6.0_rp)

                this%source_field%x(i, j, k, e) = &
                     this%cb1 * (1.0_rp - ft2) * s_tilde * nu_tilde + &
                     (this%cb2 / this%sigma_sa) * grad_sq - &
                     (this%cw1 * fw - (this%cb1 / this%kappa**2) * ft2) * &
                     (nu_tilde / d_eff)**2

                if (this%sync_nut_with_scalar) then
                   this%nut%x(i, j, k, e) = nu_tilde
                end if
             end do
          end do
       end do
    end do

    call neko_scratch_registry%relinquish_field(idx)

  end subroutine spalart_allmaras_compute


  !> Implementation of cheap_dist in Nek5000 (CPU)
  !! @param d Wall-distance field data.
  subroutine spalart_allmaras_compute_cheap_dist_cpu(this, d)
    class(spalart_allmaras_t), intent(inout) :: this
    real(kind=rp), intent(inout), target :: d(:)
    type(zero_dirichlet_t) :: bc_wall
    real(kind=rp), pointer :: d4(:, :, :, :)
    integer :: i, j, k, e, n
    integer :: ipass, nchange, max_pass
    integer :: ii, jj, kk, i0, i1, j0, j1, k0, k1
    integer :: lx, ly, lz, nel, z_idx
    integer :: m, idx
    real(kind=rp) :: dtmp, x1, y1, z1, x2, y2, z2
    integer :: change_vec(1)
    logical :: done

    lx = this%coef%dof%Xh%lx
    ly = this%coef%dof%Xh%ly
    lz = this%coef%dof%Xh%lz
    nel = this%coef%msh%nelv
    n = size(d)
    d4(1:lx, 1:ly, 1:lz, 1:nel) => d
    max_pass = 10000

    call cfill(d, huge(0.0_rp), n)

    if (size(this%wall_distance_zone_indices) .gt. 0) then
       call bc_wall%init_from_components(this%coef)
       do k = 1, size(this%wall_distance_zone_indices)
          z_idx = this%wall_distance_zone_indices(k)
          call bc_wall%mark_zone(this%coef%msh%labeled_zones(z_idx))
       end do
       call bc_wall%finalize()
       m = bc_wall%msk(0)
       do i = 1, m
          idx = bc_wall%msk(i)
          d(idx) = 0.0_rp
       end do
       call bc_wall%free()
    end if

    ipass = 1
    done = .false.
    do while (ipass <= max_pass .and. .not. done)
       nchange = 0
       do e = 1, nel
          do k = 1, lz
             do j = 1, ly
               do i = 1, lx
                   x1 = this%coef%dof%x(i, j, k, e)
                   y1 = this%coef%dof%y(i, j, k, e)
                   z1 = this%coef%dof%z(i, j, k, e)
                   i0 = max(1, i - 1)
                   i1 = min(lx, i + 1)
                   j0 = max(1, j - 1)
                   j1 = min(ly, j + 1)
                   k0 = max(1, k - 1)
                   k1 = min(lz, k + 1)
                   do kk = k0, k1
                      do jj = j0, j1
                         do ii = i0, i1
                            if (ii == i .and. jj == j .and. kk == k) then
                               cycle
                            end if
                            x2 = this%coef%dof%x(ii, jj, kk, e)
                            y2 = this%coef%dof%y(ii, jj, kk, e)
                            z2 = this%coef%dof%z(ii, jj, kk, e)
                            dtmp = d4(ii, jj, kk, e) + &
                                 sqrt((x1 - x2)**2 + &
                                 (y1 - y2)**2 + (z1 - z2)**2)
                            if (dtmp < d4(i, j, k, e)) then
                               d4(i, j, k, e) = dtmp
                               nchange = nchange + 1
                            end if
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
       call this%coef%gs_h%gs_op_vector(d, n, GS_OP_MIN)
       change_vec(1) = nchange
       if (glimax(change_vec, 1) == 0) then
          done = .true.
       end if
       ipass = ipass + 1
    end do
  end subroutine spalart_allmaras_compute_cheap_dist_cpu

end module spalart_allmaras
