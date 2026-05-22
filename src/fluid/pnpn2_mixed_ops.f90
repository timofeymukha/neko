! Copyright (c) 2025, The Neko Authors
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
!> Implements type `pnpn2_mixed_ops_t`.
module pnpn2_mixed_ops
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use math, only : add2, rzero
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use space, only : GL, GLL, space_t
  use speclib, only : DGLLGL, IGLLM
  use tensor, only : tnsr3d
  use utils, only : neko_error
  implicit none
  private

  !> Exact mixed \f$ X_h \leftrightarrow Y_h \f$ operators for Nek-style PnPn-2.
  type, public :: pnpn2_mixed_ops_t
     !> Velocity function space on the GLL mesh.
     type(space_t), pointer :: Xh => null()
     !> Pressure function space on the GL mesh.
     type(space_t), pointer :: Yh => null()
     !> Velocity coefficients on \f$ X_h \f$.
     type(coef_t), pointer :: coef_Xh => null()
     !> Pressure coefficients on \f$ Y_h \f$.
     type(coef_t), pointer :: coef_Yh => null()
     !> Nek `I12`: GLL -> GL interpolation.
     real(kind=rp), allocatable :: i12(:,:)
     !> Literal transpose of `I12`.
     real(kind=rp), allocatable :: i12t(:,:)
     !> Nek `D12`: GLL -> GL derivative.
     real(kind=rp), allocatable :: d12(:,:)
     !> Literal transpose of `D12`.
     real(kind=rp), allocatable :: d12t(:,:)
     !> Work buffer on \f$ Y_h \f$.
     type(field_t) :: work_Yh
     !> Work buffer on \f$ X_h \f$.
     type(field_t) :: work_Xh
   contains
     !> Constructor.
     procedure, pass(this) :: init => pnpn2_mixed_ops_init
     !> Destructor.
     procedure, pass(this) :: free => pnpn2_mixed_ops_free
     !> Apply the mixed weak gradient \f$ G^T : Y_h \rightarrow X_h \f$.
     procedure, pass(this) :: opgradt => pnpn2_mixed_ops_opgradt
     !> Apply the mixed weak divergence \f$ D : X_h \rightarrow Y_h \f$.
     procedure, pass(this) :: opdiv => pnpn2_mixed_ops_opdiv
  end type pnpn2_mixed_ops_t

contains

  !> Initialize the exact mixed `12`-family operators for the unequal-order path.
  subroutine pnpn2_mixed_ops_init(this, Xh, Yh, dm_Xh, dm_Yh, c_Xh, c_Yh)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    type(space_t), target, intent(inout) :: Xh
    type(space_t), target, intent(inout) :: Yh
    type(dofmap_t), intent(inout) :: dm_Xh
    type(dofmap_t), intent(inout) :: dm_Yh
    type(coef_t), target, intent(inout) :: c_Xh
    type(coef_t), target, intent(inout) :: c_Yh

    call this%free()

    if (Xh%t .ne. GLL) then
      call neko_error('PnPn-2 mixed operators require a GLL velocity space.')
    end if
    if (Yh%t .ne. GL) then
      call neko_error('PnPn-2 mixed operators require a GL pressure space.')
    end if

    allocate(this%i12(Yh%lx, Xh%lx))
    allocate(this%i12t(Xh%lx, Yh%lx))
    allocate(this%d12(Yh%lx, Xh%lx))
    allocate(this%d12t(Xh%lx, Yh%lx))

    call IGLLM(this%i12, this%i12t, Xh%zg(1,1), Yh%zg(1,1), Xh%lx, Yh%lx, &
         Xh%lx, Yh%lx)
    call DGLLGL(this%d12, this%d12t, Xh%zg(1,1), Yh%zg(1,1), this%i12, &
         Xh%lx, Yh%lx, Xh%lx, Yh%lx)

    call this%work_Yh%init(dm_Yh, 'pnpn2_work_Yh')
    call this%work_Xh%init(dm_Xh, 'pnpn2_work_Xh')

    this%Xh => Xh
    this%Yh => Yh
    this%coef_Xh => c_Xh
    this%coef_Yh => c_Yh
  end subroutine pnpn2_mixed_ops_init

  !> Free mixed-space operator storage.
  subroutine pnpn2_mixed_ops_free(this)
    class(pnpn2_mixed_ops_t), intent(inout) :: this

    if (allocated(this%i12)) then
       deallocate(this%i12)
    end if
    if (allocated(this%i12t)) then
       deallocate(this%i12t)
    end if
    if (allocated(this%d12)) then
       deallocate(this%d12)
    end if
    if (allocated(this%d12t)) then
       deallocate(this%d12t)
    end if

    call this%work_Yh%free()
    call this%work_Xh%free()

    nullify(this%Xh)
    nullify(this%Yh)
    nullify(this%coef_Xh)
    nullify(this%coef_Yh)
  end subroutine pnpn2_mixed_ops_free

  !> Apply the Nek-style mixed weak gradient from \f$ Y_h \f$ to \f$ X_h \f$.
  subroutine pnpn2_mixed_ops_opgradt(this, gx, gy, gz, p)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: gx(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(inout) :: gy(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(inout) :: gz(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: p(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2_mixed_ops_t is currently CPU-only.')
    end if

    call mixed_grad_component(this, gx, p, this%coef_Yh%drdx, this%coef_Yh%dsdx, &
         this%coef_Yh%dtdx)
    call mixed_grad_component(this, gy, p, this%coef_Yh%drdy, this%coef_Yh%dsdy, &
         this%coef_Yh%dtdy)
    if (this%coef_Xh%msh%gdim .eq. 3) then
      call mixed_grad_component(this, gz, p, this%coef_Yh%drdz, &
           this%coef_Yh%dsdz, this%coef_Yh%dtdz)
    else
      call rzero(gz, this%coef_Xh%dof%size())
    end if
  end subroutine pnpn2_mixed_ops_opgradt

  !> Apply the Nek-style mixed weak divergence from \f$ X_h \f$ to \f$ Y_h \f$.
  subroutine pnpn2_mixed_ops_opdiv(this, div, u, v, w)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: div(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: u(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: v(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: w(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_error('pnpn2_mixed_ops_t is currently CPU-only.')
    end if

    call rzero(div, this%coef_Yh%dof%size())

    call mixed_div_component(this, div, u, this%coef_Yh%drdx, this%coef_Yh%dsdx, &
         this%coef_Yh%dtdx)
    call mixed_div_component(this, div, v, this%coef_Yh%drdy, this%coef_Yh%dsdy, &
         this%coef_Yh%dtdy)
    if (this%coef_Xh%msh%gdim .eq. 3) then
      call mixed_div_component(this, div, w, this%coef_Yh%drdz, &
           this%coef_Yh%dsdz, this%coef_Yh%dtdz)
    end if
  end subroutine pnpn2_mixed_ops_opdiv

  !> Apply one Cartesian component of Nek's `opgradt`.
  subroutine mixed_grad_component(this, out, p, dr, ds, dt)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: out(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: p(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: dr(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: ds(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: dt(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    integer :: n_x, n_y

    n_x = this%coef_Xh%dof%size()
    n_y = this%coef_Yh%dof%size()

    call rzero(out, n_x)

    call build_weighted_pressure(this, p, dr)
    call tnsr3d(out, this%Xh%lx, this%work_Yh%x, this%Yh%lx, this%d12t, &
         this%i12, this%i12, this%coef_Xh%msh%nelv)

    call build_weighted_pressure(this, p, ds)
    call tnsr3d(this%work_Xh%x, this%Xh%lx, this%work_Yh%x, this%Yh%lx, &
         this%i12t, this%d12, this%i12, this%coef_Xh%msh%nelv)
    call add2(out, this%work_Xh%x, n_x)

    if (this%coef_Xh%msh%gdim .eq. 3) then
      call build_weighted_pressure(this, p, dt)
      call tnsr3d(this%work_Xh%x, this%Xh%lx, this%work_Yh%x, this%Yh%lx, &
           this%i12t, this%i12, this%d12, this%coef_Xh%msh%nelv)
      call add2(out, this%work_Xh%x, n_x)
    end if
  end subroutine mixed_grad_component

  !> Apply one Cartesian component of Nek's `opdiv`.
  subroutine mixed_div_component(this, div, vel, dr, ds, dt)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: div(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: vel(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: dr(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: ds(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: dt(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    integer :: n_y

    n_y = this%coef_Yh%dof%size()

    call tnsr3d(this%work_Yh%x, this%Yh%lx, vel, this%Xh%lx, this%d12, &
         this%i12t, this%i12t, this%coef_Xh%msh%nelv)
    call accumulate_weighted_divergence(this, div, this%work_Yh%x, dr, n_y)

    call tnsr3d(this%work_Yh%x, this%Yh%lx, vel, this%Xh%lx, this%i12, &
         this%d12t, this%i12t, this%coef_Xh%msh%nelv)
    call accumulate_weighted_divergence(this, div, this%work_Yh%x, ds, n_y)

    if (this%coef_Xh%msh%gdim .eq. 3) then
      call tnsr3d(this%work_Yh%x, this%Yh%lx, vel, this%Xh%lx, this%i12, &
           this%i12t, this%d12t, this%coef_Xh%msh%nelv)
      call accumulate_weighted_divergence(this, div, this%work_Yh%x, dt, n_y)
    end if
  end subroutine mixed_div_component

  !> Form the weighted pressure-grid scalar used by Nek `cdtp`.
  subroutine build_weighted_pressure(this, p, metric)
    class(pnpn2_mixed_ops_t), intent(inout) :: this
    real(kind=rp), intent(in) :: p(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    real(kind=rp), intent(in) :: metric(this%Yh%lx, this%Yh%ly, this%Yh%lz, &
         this%coef_Xh%msh%nelv)
    integer :: i, n_y

    n_y = this%coef_Yh%dof%size()
    do concurrent (i = 1:n_y)
      this%work_Yh%x(i,1,1,1) = this%coef_Yh%B(i,1,1,1) * &
           this%coef_Yh%jacinv(i,1,1,1) * p(i,1,1,1) * metric(i,1,1,1)
    end do
  end subroutine build_weighted_pressure

  !> Add one weighted pressure-grid contribution to Nek `opdiv`.
  subroutine accumulate_weighted_divergence(this, div, contribution, metric, n_y)
    class(pnpn2_mixed_ops_t), intent(in) :: this
    real(kind=rp), intent(inout) :: div(:,:,:,:)
    real(kind=rp), intent(in) :: contribution(:,:,:,:)
    real(kind=rp), intent(in) :: metric(:,:,:,:)
    integer, intent(in) :: n_y
    integer :: i

    do concurrent (i = 1:n_y)
      div(i,1,1,1) = div(i,1,1,1) + this%coef_Yh%B(i,1,1,1) * &
           this%coef_Yh%jacinv(i,1,1,1) * contribution(i,1,1,1) * &
           metric(i,1,1,1)
    end do
  end subroutine accumulate_weighted_divergence

end module pnpn2_mixed_ops
