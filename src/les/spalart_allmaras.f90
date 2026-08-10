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
!> Implements `spalart_allmaras_t`.
module spalart_allmaras
  use num_types, only : rp
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use les_model, only : les_model_t
  use spalart_allmaras_model, only : spalart_allmaras_cw1, &
       spalart_allmaras_pointwise, sa_default_cb1, sa_default_cb2, &
       sa_default_sigma, sa_default_cw2, sa_default_cw3, sa_default_cv1, &
       sa_default_ct3, sa_default_ct4, sa_default_kappa, sa_default_cn1
  use operators, only : dudxyz
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup_or_default
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

  !> Implements Spalart--Allmaras RANS model variants.
  !! @note The transported working-variable equation is solved as a scalar.
  type, public, extends(les_model_t) :: spalart_allmaras_t
     !> Name of the scalar field carrying the SA working variable.
     character(len=:), allocatable :: scalar_field_name
     !> Name of the field containing the non-diffusive SA source.
     character(len=:), allocatable :: source_field_name
     !> Name of the field containing nu_tilde / sigma.
     character(len=:), allocatable :: alphat_field_name
     !> Minimum wall distance used in source term evaluation.
     real(kind=rp) :: min_wall_distance = 1.0e-12_rp
     !> Standard-SA model constants.
     real(kind=rp) :: cb1 = sa_default_cb1
     real(kind=rp) :: cb2 = sa_default_cb2
     real(kind=rp) :: sigma_sa = sa_default_sigma
     real(kind=rp) :: cw1
     real(kind=rp) :: cw2 = sa_default_cw2
     real(kind=rp) :: cw3 = sa_default_cw3
     real(kind=rp) :: cv1 = sa_default_cv1
     real(kind=rp) :: ct3 = sa_default_ct3
     real(kind=rp) :: ct4 = sa_default_ct4
     real(kind=rp) :: kappa = sa_default_kappa
     real(kind=rp) :: cn1 = sa_default_cn1
     !> Whether to evaluate SA-neg terms at negative working-variable nodes.
     logical :: negative_treatment = .false.
     !> Zone indices for the temporary graph-based wall distance.
     integer, allocatable :: wall_distance_zone_indices(:)
     !> Pointwise non-diffusive source for the scalar equation.
     type(field_t), pointer :: source_field => null()
     !> Pointwise turbulent diffusivity for the scalar equation.
     type(field_t), pointer :: alphat => null()
     !> Dynamic viscosity.
     type(field_t), pointer :: mu => null()
     !> Density.
     type(field_t), pointer :: rho => null()
   contains
     procedure, pass(this) :: init => spalart_allmaras_init
     procedure, pass(this) :: init_from_components => &
          spalart_allmaras_init_from_components
     procedure, pass(this) :: free => spalart_allmaras_free
     procedure, pass(this) :: compute => spalart_allmaras_compute
     procedure, pass(this) :: compute_cheap_dist => &
          spalart_allmaras_compute_cheap_dist_cpu
  end type spalart_allmaras_t

contains

  !> Construct the model from a JSON object.
  subroutine spalart_allmaras_init(this, fluid, json)
    class(spalart_allmaras_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name, scalar_field_name
    character(len=:), allocatable :: source_field_name, alphat_field_name
    character(len=:), allocatable :: delta_type, variant
    integer, allocatable :: wall_distance_zone_indices(:)
    logical :: if_ext, negative_treatment, no_ft2, ct3_is_set
    real(kind=rp) :: cb1, cb2, sigma_sa, cw1, min_wall_distance
    real(kind=rp) :: cw2, cw3, cv1, ct3, ct4, kappa, cn1
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "scalar_field", scalar_field_name, &
         "sa_nu_tilde")
    call json_get_or_default(json, "source_field", source_field_name, &
         "sa_source")
    call json_get_or_default(json, "alphat_field", alphat_field_name, &
         "sa_alphat")
    call json_get_or_default(json, "variant", variant, "sa")
    no_ft2 = .false.
    negative_treatment = .false.
    select case (trim(variant))
    case ("sa")
    case ("sa-neg")
       negative_treatment = .true.
    case ("sa-noft2")
       no_ft2 = .true.
    case ("sa-noft2-neg")
       no_ft2 = .true.
       negative_treatment = .true.
    case default
       call neko_error("Unknown Spalart-Allmaras variant: " // &
            trim(variant))
    end select
    if (json%valid_path("wall_distance_zone_indices")) then
       call json_get(json, "wall_distance_zone_indices", &
            wall_distance_zone_indices)
    else
       allocate(wall_distance_zone_indices(0))
    end if
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "extrapolation", if_ext, .false.)
    call json_get_or_lookup_or_default(json, "cb1", cb1, sa_default_cb1)
    call json_get_or_lookup_or_default(json, "cb2", cb2, sa_default_cb2)
    call json_get_or_lookup_or_default(json, "sigma", sigma_sa, &
         sa_default_sigma)
    call json_get_or_lookup_or_default(json, "cw2", cw2, sa_default_cw2)
    call json_get_or_lookup_or_default(json, "cw3", cw3, sa_default_cw3)
    call json_get_or_lookup_or_default(json, "cv1", cv1, sa_default_cv1)
    ct3_is_set = json%valid_path("ct3")
    call json_get_or_lookup_or_default(json, "ct3", ct3, sa_default_ct3)
    call json_get_or_lookup_or_default(json, "ct4", ct4, sa_default_ct4)
    call json_get_or_lookup_or_default(json, "kappa", kappa, &
         sa_default_kappa)
    call json_get_or_lookup_or_default(json, "cn1", cn1, sa_default_cn1)
    call json_get_or_lookup_or_default(json, "min_wall_distance", &
         min_wall_distance, 1.0e-12_rp)

    if (no_ft2) then
       if (abs(ct3) .gt. tiny(1.0_rp) .and. ct3_is_set) then
          call neko_error("SA-noft2 variants require ct3 = 0")
       end if
       ct3 = 0.0_rp
    end if

    if (sigma_sa .le. 0.0_rp .or. kappa .le. 0.0_rp .or. &
         cv1 .le. 0.0_rp .or. cw3 .le. 0.0_rp .or. &
         cn1 .le. 0.0_rp .or. min_wall_distance .le. 0.0_rp) then
       call neko_error("Spalart-Allmaras constants and minimum distance " // &
            "must be positive")
    end if
    cw1 = spalart_allmaras_cw1(cb1, cb2, sigma_sa, kappa)

    call neko_log%section("Spalart-Allmaras model")
    write(log_buf, '(A,A)') "Variant : ", variant
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Working variable : ", scalar_field_name
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Eddy viscosity : ", nut_name
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Scalar diffusivity : ", alphat_field_name
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Scalar source : ", source_field_name
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') "cw1 : ", cw1
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call this%init_from_components(fluid, nut_name, scalar_field_name, &
         source_field_name, alphat_field_name, wall_distance_zone_indices, &
         delta_type, if_ext, cb1, cb2, sigma_sa, cw1, cw2, cw3, cv1, ct3, &
         ct4, kappa, cn1, min_wall_distance, negative_treatment)
  end subroutine spalart_allmaras_init

  !> Construct the model from individual components.
  subroutine spalart_allmaras_init_from_components(this, fluid, nut_name, &
       scalar_field_name, source_field_name, alphat_field_name, &
       wall_distance_zone_indices, delta_type, if_ext, cb1, cb2, sigma_sa, &
       cw1, cw2, cw3, cv1, ct3, ct4, kappa, cn1, min_wall_distance, &
       negative_treatment)
    class(spalart_allmaras_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    character(len=*), intent(in) :: nut_name, scalar_field_name
    character(len=*), intent(in) :: source_field_name, alphat_field_name
    integer, intent(in) :: wall_distance_zone_indices(:)
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: cb1, cb2, sigma_sa, cw1, cw2, cw3
    real(kind=rp), intent(in) :: cv1, ct3, ct4, kappa, cn1
    real(kind=rp), intent(in) :: min_wall_distance
    logical, intent(in) :: negative_treatment
    type(field_t), pointer :: wall_distance
    real(kind=rp), pointer :: wall_distance_vec(:)

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Spalart-Allmaras is not implemented for devices")
    end if

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%scalar_field_name = trim(scalar_field_name)
    this%source_field_name = trim(source_field_name)
    this%alphat_field_name = trim(alphat_field_name)
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
    this%cn1 = cn1
    this%min_wall_distance = min_wall_distance
    this%negative_treatment = negative_treatment
    allocate(this%wall_distance_zone_indices( &
         size(wall_distance_zone_indices)))
    this%wall_distance_zone_indices = wall_distance_zone_indices
    this%mu => fluid%mu
    this%rho => fluid%rho

    call neko_registry%add_field(fluid%dm_Xh, &
         trim(this%source_field_name), ignore_existing = .true.)
    call neko_registry%add_field(fluid%dm_Xh, &
         trim(this%alphat_field_name), ignore_existing = .true.)
    call neko_registry%add_field(fluid%dm_Xh, "wall_distance", &
         ignore_existing = .true.)

    this%source_field => neko_registry%get_field_by_name( &
         trim(this%source_field_name))
    this%alphat => neko_registry%get_field_by_name( &
         trim(this%alphat_field_name))
    wall_distance => neko_registry%get_field_by_name("wall_distance")

    if (size(this%wall_distance_zone_indices) .eq. 0) then
       call neko_error("Spalart-Allmaras requires wall zone indices")
    end if

    wall_distance_vec(1:wall_distance%size()) => wall_distance%x
    call this%compute_cheap_dist(wall_distance_vec)
  end subroutine spalart_allmaras_init_from_components

  !> Free the model.
  subroutine spalart_allmaras_free(this)
    class(spalart_allmaras_t), intent(inout) :: this

    nullify(this%source_field)
    nullify(this%alphat)
    nullify(this%mu)
    nullify(this%rho)
    if (allocated(this%scalar_field_name)) &
         deallocate(this%scalar_field_name)
    if (allocated(this%source_field_name)) &
         deallocate(this%source_field_name)
    if (allocated(this%alphat_field_name)) &
         deallocate(this%alphat_field_name)
    if (allocated(this%wall_distance_zone_indices)) &
         deallocate(this%wall_distance_zone_indices)

    call this%free_base()
  end subroutine spalart_allmaras_free

  !> Update eddy viscosity, scalar diffusivity, and scalar source fields.
  subroutine spalart_allmaras_compute(this, t, tstep)
    class(spalart_allmaras_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u, v, w, u_work, v_work, w_work
    type(field_t), pointer :: u_e, v_e, w_e, scalar, wall_distance
    type(field_t), pointer :: du_dy, du_dz, dv_dx, dv_dz, dw_dx, dw_dy
    type(field_t), pointer :: dsa_dx, dsa_dy, dsa_dz
    integer :: idx(9)
    integer :: e, i, j, k
    real(kind=rp) :: omega_x, omega_y, omega_z, omega
    real(kind=rp) :: grad_sq, d_wall, d_eff, nu_tilde, nu, sa_source

    if (this%if_ext) then
       associate(ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag, ext_bdf => this%ext_bdf)
         u => neko_registry%get_field_by_name("u")
         v => neko_registry%get_field_by_name("v")
         w => neko_registry%get_field_by_name("w")
         u_e => neko_registry%get_field_by_name("u_e")
         v_e => neko_registry%get_field_by_name("v_e")
         w_e => neko_registry%get_field_by_name("w_e")
         call this%sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
              ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)
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

    scalar => neko_registry%get_field_by_name( &
         trim(this%scalar_field_name))
    wall_distance => neko_registry%get_field_by_name("wall_distance")
    if (.not. associated(this%mu) .or. .not. associated(this%rho)) then
       call neko_error("Spalart-Allmaras requires viscosity and density")
    end if

    call neko_scratch_registry%request_field(du_dy, idx(1), .false.)
    call neko_scratch_registry%request_field(du_dz, idx(2), .false.)
    call neko_scratch_registry%request_field(dv_dx, idx(3), .false.)
    call neko_scratch_registry%request_field(dv_dz, idx(4), .false.)
    call neko_scratch_registry%request_field(dw_dx, idx(5), .false.)
    call neko_scratch_registry%request_field(dw_dy, idx(6), .false.)
    call neko_scratch_registry%request_field(dsa_dx, idx(7), .false.)
    call neko_scratch_registry%request_field(dsa_dy, idx(8), .false.)
    call neko_scratch_registry%request_field(dsa_dz, idx(9), .false.)

    call dudxyz(du_dy%x, u_work%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(du_dz%x, u_work%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)
    call dudxyz(dv_dx%x, v_work%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dv_dz%x, v_work%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)
    call dudxyz(dw_dx%x, w_work%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dw_dy%x, w_work%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(dsa_dx%x, scalar%x, this%coef%drdx, this%coef%dsdx, &
         this%coef%dtdx, this%coef)
    call dudxyz(dsa_dy%x, scalar%x, this%coef%drdy, this%coef%dsdy, &
         this%coef%dtdy, this%coef)
    call dudxyz(dsa_dz%x, scalar%x, this%coef%drdz, this%coef%dsdz, &
         this%coef%dtdz, this%coef)

    do e = 1, this%coef%msh%nelv
       do k = 1, this%coef%Xh%lz
          do j = 1, this%coef%Xh%ly
             do i = 1, this%coef%Xh%lx
                omega_x = dw_dy%x(i,j,k,e) - dv_dz%x(i,j,k,e)
                omega_y = du_dz%x(i,j,k,e) - dw_dx%x(i,j,k,e)
                omega_z = dv_dx%x(i,j,k,e) - du_dy%x(i,j,k,e)
                omega = sqrt(max(0.0_rp, omega_x**2 + omega_y**2 + &
                     omega_z**2))
                grad_sq = dsa_dx%x(i,j,k,e)**2 + &
                     dsa_dy%x(i,j,k,e)**2 + dsa_dz%x(i,j,k,e)**2

                d_wall = abs(wall_distance%x(i,j,k,e))
                d_eff = max(d_wall, this%min_wall_distance)
                nu_tilde = scalar%x(i,j,k,e)
                if (d_wall .le. this%min_wall_distance) nu_tilde = 0.0_rp
                nu = this%mu%x(i,j,k,e) / this%rho%x(i,j,k,e)

                call spalart_allmaras_pointwise(nu_tilde, nu, omega, &
                     d_eff, grad_sq, this%cb1, this%cb2, this%sigma_sa, &
                     this%cw1, this%cw2, this%cw3, this%cv1, this%ct3, &
                     this%ct4, this%kappa, this%cn1, &
                     this%negative_treatment, this%nut%x(i,j,k,e), &
                     this%alphat%x(i,j,k,e), sa_source)
                this%source_field%x(i,j,k,e) = &
                     this%rho%x(i,j,k,e) * sa_source
             end do
          end do
       end do
    end do

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine spalart_allmaras_compute

  !> Graph-based wall-distance approximation retained temporarily.
  subroutine spalart_allmaras_compute_cheap_dist_cpu(this, d)
    class(spalart_allmaras_t), intent(inout) :: this
    real(kind=rp), intent(inout), target :: d(:)
    type(zero_dirichlet_t) :: bc_wall
    real(kind=rp), pointer :: d4(:, :, :, :)
    integer :: i, j, k, e, n, ipass, nchange, max_pass
    integer :: ii, jj, kk, i0, i1, j0, j1, k0, k1
    integer :: lx, ly, lz, nel, z_idx, m, idx
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

    ipass = 1
    done = .false.
    do while (ipass .le. max_pass .and. .not. done)
       nchange = 0
       do e = 1, nel
          do k = 1, lz
             do j = 1, ly
                do i = 1, lx
                   x1 = this%coef%dof%x(i,j,k,e)
                   y1 = this%coef%dof%y(i,j,k,e)
                   z1 = this%coef%dof%z(i,j,k,e)
                   i0 = max(1, i - 1)
                   i1 = min(lx, i + 1)
                   j0 = max(1, j - 1)
                   j1 = min(ly, j + 1)
                   k0 = max(1, k - 1)
                   k1 = min(lz, k + 1)
                   do kk = k0, k1
                      do jj = j0, j1
                         do ii = i0, i1
                            if (ii .eq. i .and. jj .eq. j .and. &
                                 kk .eq. k) cycle
                            x2 = this%coef%dof%x(ii,jj,kk,e)
                            y2 = this%coef%dof%y(ii,jj,kk,e)
                            z2 = this%coef%dof%z(ii,jj,kk,e)
                            dtmp = d4(ii,jj,kk,e) + sqrt((x1 - x2)**2 &
                                 + (y1 - y2)**2 + (z1 - z2)**2)
                            if (dtmp .lt. d4(i,j,k,e)) then
                               d4(i,j,k,e) = dtmp
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
       if (glimax(change_vec, 1) .eq. 0) done = .true.
       ipass = ipass + 1
    end do
  end subroutine spalart_allmaras_compute_cheap_dist_cpu

end module spalart_allmaras
