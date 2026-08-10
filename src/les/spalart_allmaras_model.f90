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
!> Pointwise algebra for Spalart--Allmaras model variants.
module spalart_allmaras_model
  use num_types, only : rp
  implicit none
  private

  real(kind=rp), parameter, public :: sa_default_cb1 = 0.1355_rp
  real(kind=rp), parameter, public :: sa_default_cb2 = 0.622_rp
  real(kind=rp), parameter, public :: sa_default_sigma = 2.0_rp / 3.0_rp
  real(kind=rp), parameter, public :: sa_default_cw2 = 0.3_rp
  real(kind=rp), parameter, public :: sa_default_cw3 = 2.0_rp
  real(kind=rp), parameter, public :: sa_default_cv1 = 7.1_rp
  real(kind=rp), parameter, public :: sa_default_ct3 = 1.2_rp
  real(kind=rp), parameter, public :: sa_default_ct4 = 0.5_rp
  real(kind=rp), parameter, public :: sa_default_kappa = 0.41_rp
  real(kind=rp), parameter, public :: sa_default_cn1 = 16.0_rp
  real(kind=rp), parameter, public :: sa_default_c2 = 0.7_rp
  real(kind=rp), parameter, public :: sa_default_c3 = 0.9_rp

  public :: spalart_allmaras_cw1, spalart_allmaras_pointwise

contains

  !> Compute the dependent standard-SA constant cw1.
  pure function spalart_allmaras_cw1(cb1, cb2, sigma, kappa) result(cw1)
    real(kind=rp), intent(in) :: cb1, cb2, sigma, kappa
    real(kind=rp) :: cw1

    cw1 = cb1 / kappa**2 + (1.0_rp + cb2) / sigma
  end function spalart_allmaras_cw1

  !> Evaluate all pointwise SA terms.
  !!
  !! The returned source excludes the variable-diffusion divergence term,
  !! which is handled by the scalar solver using `alphat`. It includes the
  !! non-conservative cb2 gradient term. When `negative_treatment` is true,
  !! the SA-neg terms are evaluated locally wherever `nu_tilde` is negative.
  pure subroutine spalart_allmaras_pointwise(nu_tilde, nu, omega, d, &
       grad_sq, cb1, cb2, sigma, cw1, cw2, cw3, cv1, ct3, ct4, kappa, &
       cn1, negative_treatment, nut, alphat, source, fv1_out, fv2_out, &
       ft2_out, s_tilde_out, r_out, fw_out, production_out, &
       cross_diffusion_out, destruction_out)
    real(kind=rp), intent(in) :: nu_tilde, nu, omega, d, grad_sq
    real(kind=rp), intent(in) :: cb1, cb2, sigma, cw1, cw2, cw3
    real(kind=rp), intent(in) :: cv1, ct3, ct4, kappa, cn1
    logical, intent(in) :: negative_treatment
    real(kind=rp), intent(out) :: nut, alphat, source
    real(kind=rp), intent(out), optional :: fv1_out, fv2_out, ft2_out
    real(kind=rp), intent(out), optional :: s_tilde_out, r_out, fw_out
    real(kind=rp), intent(out), optional :: production_out
    real(kind=rp), intent(out), optional :: cross_diffusion_out
    real(kind=rp), intent(out), optional :: destruction_out
    real(kind=rp) :: chi, chi3, fn, fv1, fv2, ft2
    real(kind=rp) :: s_bar, s_tilde, denom, r, g, fw
    real(kind=rp) :: production, cross_diffusion, destruction

    chi = nu_tilde / nu
    chi3 = chi**3
    cross_diffusion = (cb2 / sigma) * grad_sq

    if (negative_treatment .and. nu_tilde .lt. 0.0_rp) then
       ! SA-neg is a local continuation of SA through negative transient
       ! states. The signed destruction quantity is subtracted below, so it
       ! is negative on this branch and contributes positively to the source.
       fn = (cn1 + chi3) / (cn1 - chi3)
       fv1 = 0.0_rp
       fv2 = 0.0_rp
       ft2 = ct3
       s_tilde = omega
       r = 0.0_rp
       fw = 0.0_rp
       nut = 0.0_rp
       alphat = nu_tilde * fn / sigma
       production = cb1 * (1.0_rp - ct3) * omega * nu_tilde
       destruction = -cw1 * (nu_tilde / d)**2
    else
       fv1 = chi3 / (chi3 + cv1**3)
       fv2 = 1.0_rp - chi / (1.0_rp + chi * fv1)
       ft2 = ct3 * exp(-ct4 * chi**2)

       ! NASA TMR Note 1(c). This piecewise continuation replaces S_tilde,
       ! rather than merely clipping it in the denominator of r. Therefore
       ! the same modified value is used by both production and destruction.
       s_bar = nu_tilde * fv2 / (kappa**2 * d**2)
       if (s_bar .ge. -sa_default_c2 * omega) then
          s_tilde = omega + s_bar
       else
          denom = (sa_default_c3 - 2.0_rp * sa_default_c2) * omega &
               - s_bar
          s_tilde = omega + omega * &
               (sa_default_c2**2 * omega + sa_default_c3 * s_bar) / denom
       end if

       denom = s_tilde * kappa**2 * d**2
       if (denom .gt. tiny(1.0_rp)) then
          r = min(nu_tilde / denom, 10.0_rp)
       else
          r = 10.0_rp
       end if

       g = r + cw2 * (r**6 - r)
       fw = g * ((1.0_rp + cw3**6) / &
            (g**6 + cw3**6))**(1.0_rp / 6.0_rp)

       nut = nu_tilde * fv1
       alphat = nu_tilde / sigma
       production = cb1 * (1.0_rp - ft2) * s_tilde * nu_tilde
       destruction = (cw1 * fw - (cb1 / kappa**2) * ft2) * &
            (nu_tilde / d)**2
    end if
    source = production + cross_diffusion - destruction

    if (present(fv1_out)) fv1_out = fv1
    if (present(fv2_out)) fv2_out = fv2
    if (present(ft2_out)) ft2_out = ft2
    if (present(s_tilde_out)) s_tilde_out = s_tilde
    if (present(r_out)) r_out = r
    if (present(fw_out)) fw_out = fw
    if (present(production_out)) production_out = production
    if (present(cross_diffusion_out)) &
         cross_diffusion_out = cross_diffusion
    if (present(destruction_out)) destruction_out = destruction
  end subroutine spalart_allmaras_pointwise

end module spalart_allmaras_model
