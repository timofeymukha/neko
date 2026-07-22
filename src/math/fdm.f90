! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
!> Type for the Fast Diagonalization connected with the schwarz overlapping solves.
module fdm
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM
  use num_types, only : rp, sp, dp, qp
  use mesh, only : mesh_t
  use space, only : space_t
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use fdm_cpu, only : fdm_do_fast_cpu
  use fdm_device, only : fdm_do_fast_device
  use fdm_sx, only : fdm_do_fast_sx
  use fdm_xsmm, only : fdm_do_fast_xsmm
  use utils, only : neko_error, neko_warning
  use comm, only : pe_rank
  use math, only : vlmax
  use device, only : glb_cmd_queue, DEVICE_TO_HOST, HOST_TO_DEVICE, &
       device_memcpy, device_map, device_free
  use fast3d, only : semhat
  use tensor, only : trsp
  use math, only : rzero, row_zero
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: fdm_t
     real(kind=rp), allocatable :: s(:,:,:,:)
     real(kind=rp), allocatable :: d(:,:)
     type(c_ptr) :: s_d = C_NULL_PTR
     type(c_ptr) :: d_d = C_NULL_PTR
     real(kind=rp), allocatable :: len_lr(:), len_ls(:), len_lt(:)
     real(kind=rp), allocatable :: len_mr(:), len_ms(:), len_mt(:)
     real(kind=rp), allocatable :: len_rr(:), len_rs(:), len_rt(:)
     real(kind=rp), allocatable :: swplen(:,:,:,:)
     type(c_ptr) :: swplen_d = C_NULL_PTR
     type(space_t), pointer :: Xh => null()
     type(dofmap_t), pointer :: dof => null()
     type(gs_t), pointer :: gs_h => null()
     type(mesh_t), pointer :: msh => null()
    integer :: nl = 0
   contains
     procedure, pass(this) :: init => fdm_init
     procedure, pass(this) :: init_sem => fdm_init_sem
     procedure, pass(this) :: free => fdm_free
     procedure, pass(this) :: compute => fdm_compute
  end type fdm_t

  interface sygv
     module procedure sp_sygv, dp_sygv, qp_sygv
  end interface sygv

contains

  subroutine fdm_init(this, Xh, dof, gs_h, length_fdm)
    class(fdm_t), intent(inout) :: this
    type(space_t), target, intent(inout) :: Xh
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs_h
    type(fdm_t), intent(in), optional :: length_fdm
    !We only really use ah, bh
    real(kind=rp), dimension((Xh%lx)**2) :: ah, bh, ch, dh, zh
    real(kind=rp), dimension((Xh%lx)**2) :: dph, jph, bgl, zglhat, dgl, jgl, wh
    integer :: nl, n, nelv

    n = Xh%lx - 1 !Polynomnial degree
    nl = Xh%lx + 2 !Schwarz!
    nelv = dof%msh%nelv
    call fdm_free(this)
    this%nl = nl
    allocate(this%s(nl*nl, 2, dof%msh%gdim, dof%msh%nelv))
    allocate(this%d(nl**3, dof%msh%nelv))
    allocate(this%swplen(Xh%lx, Xh%lx, Xh%lx, dof%msh%nelv))
    allocate(this%len_lr(nelv), this%len_ls(nelv), this%len_lt(nelv))
    allocate(this%len_mr(nelv), this%len_ms(nelv), this%len_mt(nelv))
    allocate(this%len_rr(nelv), this%len_rs(nelv), this%len_rt(nelv))

    ! Zeroing here enables easier debugging since then
    ! MPI messages in GS are deterministic
    call rzero(this%swplen, Xh%lxyz * dof%msh%nelv)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%s, this%s_d, nl*nl*2*dof%msh%gdim*dof%msh%nelv)
       call device_map(this%d, this%d_d, nl**dof%msh%gdim*dof%msh%nelv)
       call device_map(this%swplen, this%swplen_d, Xh%lxyz*dof%msh%nelv)
    end if

    call semhat(ah, bh, ch, dh, zh, dph, jph, bgl, zglhat, dgl, jgl, n, wh)
    this%Xh => Xh
    this%dof => dof
    this%gs_h => gs_h
    this%msh => dof%msh

    if (present(length_fdm)) then
       this%len_lr = length_fdm%len_lr
       this%len_ls = length_fdm%len_ls
       this%len_lt = length_fdm%len_lt
       this%len_mr = length_fdm%len_mr
       this%len_ms = length_fdm%len_ms
       this%len_mt = length_fdm%len_mt
       this%len_rr = length_fdm%len_rr
       this%len_rs = length_fdm%len_rs
       this%len_rt = length_fdm%len_rt
    else
       call swap_lengths(this, dof%x, dof%y, dof%z, dof%msh%nelv, &
            dof%msh%gdim)
    end if

    call fdm_setup_fast(this, ah, bh, nl, n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%s, this%s_d, &
            nl*nl*2*dof%msh%gdim*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%d, this%d_d, &
            nl**dof%msh%gdim*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%swplen, this%swplen_d, &
            Xh%lxyz*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
    end if
  end subroutine fdm_init

  subroutine fdm_init_sem(this, Xh, dof, gs_h)
   class(fdm_t), intent(inout) :: this
   type(space_t), target, intent(inout) :: Xh
   type(dofmap_t), target, intent(in) :: dof
   type(gs_t), target, intent(inout) :: gs_h
   integer :: nl, n, nelv
   real(kind=rp) :: ah(0:Xh%lx-1, 0:Xh%lx-1), bh(0:Xh%lx-1)
   real(kind=rp) :: ch(0:Xh%lx-1, 0:Xh%lx-1), dh(0:Xh%lx-1, 0:Xh%lx-1)
   real(kind=rp) :: zh(0:Xh%lx-1), wh(0:2*Xh%lx-1)
   real(kind=rp) :: dph(0:Xh%lx-1, 1:Xh%lx-2), jph(0:Xh%lx-1, 1:Xh%lx-2)
   real(kind=rp) :: bgl(1:Xh%lx-2), zglhat(1:Xh%lx-2)
   real(kind=rp) :: dgl(1:Xh%lx-2, 0:Xh%lx-1), jgl(1:Xh%lx-2, 0:Xh%lx-1)

   n = Xh%lx - 1
   nl = Xh%lx
   nelv = dof%msh%nelv
   call fdm_free(this)
   this%nl = nl
   allocate(this%s(nl*nl, 2, dof%msh%gdim, dof%msh%nelv))
   allocate(this%d(nl**3, dof%msh%nelv))
   allocate(this%swplen(Xh%lx, Xh%ly, Xh%lz, dof%msh%nelv))
   allocate(this%len_lr(nelv), this%len_ls(nelv), this%len_lt(nelv))
   allocate(this%len_mr(nelv), this%len_ms(nelv), this%len_mt(nelv))
   allocate(this%len_rr(nelv), this%len_rs(nelv), this%len_rt(nelv))

   call rzero(this%swplen, Xh%lxyz * dof%msh%nelv)

   if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_map(this%s, this%s_d, nl*nl*2*dof%msh%gdim*dof%msh%nelv)
      call device_map(this%d, this%d_d, nl**dof%msh%gdim*dof%msh%nelv)
      call device_map(this%swplen, this%swplen_d, Xh%lxyz*dof%msh%nelv)
   end if

   call semhat(ah, bh, ch, dh, zh, dph, jph, bgl, zglhat, dgl, jgl, n, wh)
   call fdm_do_semhat_weight(jgl, dgl, bgl, n)
   this%Xh => Xh
   this%dof => dof
   this%gs_h => gs_h
   this%msh => dof%msh

   call swap_lengths(this, dof%x, dof%y, dof%z, dof%msh%nelv, dof%msh%gdim)
   call fdm_setup_fast_sem(this, bh, dgl, jgl, nl, n)

   if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(this%s, this%s_d, &
           nl*nl*2*dof%msh%gdim*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
      call device_memcpy(this%d, this%d_d, &
           nl**dof%msh%gdim*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
      call device_memcpy(this%swplen, this%swplen_d, &
           Xh%lxyz*dof%msh%nelv, HOST_TO_DEVICE, sync = .false.)
   end if
  end subroutine fdm_init_sem

  subroutine swap_lengths(this, x, y, z, nelv, gdim)
    type(fdm_t), intent(inout) :: this
    integer, intent(in) :: gdim, nelv
    real(kind=rp), dimension(this%Xh%lxyz, nelv), intent(in) :: x, y, z
    integer :: j, k, e, n2, nz0, nzn, nx, lx1, n

    associate(l => this%swplen, Xh => this%Xh, &
         llr => this%len_lr, lls => this%len_ls, llt => this%len_lt, &
         lmr => this%len_mr, lms => this%len_ms, lmt => this%len_mt, &
         lrr => this%len_rr, lrs => this%len_rs, lrt => this%len_rt)
      lx1 = this%Xh%lx
      n2 = lx1 - 1
      nz0 = 1
      nzn = 1
      nx = lx1 - 2
      if (gdim .eq. 3) then
         nz0 = 0
         nzn = n2
      end if
      call plane_space(lmr, lms, lmt, 0, n2, Xh%wx, x, y, z, &
           nx, n2, nz0, nzn, nelv, gdim)
      n = n2 + 1
      if (gdim .eq. 3) then
         do e = 1, nelv
            do j = 2, n2
               do k = 2, n2
                  l(1, k, j, e) = lmr(e)
                  l(n, k, j, e) = lmr(e)
                  l(k, 1, j, e) = lms(e)
                  l(k, n, j, e) = lms(e)
                  l(k, j, 1, e) = lmt(e)
                  l(k, j, n, e) = lmt(e)
               end do
            end do
         end do
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_memcpy(l, this%swplen_d, this%dof%size(), &
                 HOST_TO_DEVICE, sync = .false.)
            call this%gs_h%op(l, this%dof%size(), GS_OP_ADD)
            call device_memcpy(l, this%swplen_d, this%dof%size(), &
                 DEVICE_TO_HOST, sync = .true.)
         else
            call this%gs_h%op(l, this%dof%size(), GS_OP_ADD)
         end if

         do e = 1, nelv
            llr(e) = l(1, 2, 2, e) - lmr(e)
            lrr(e) = l(n, 2, 2, e) - lmr(e)
            lls(e) = l(2, 1, 2, e) - lms(e)
            lrs(e) = l(2, n, 2, e) - lms(e)
            llt(e) = l(2, 2, 1, e) - lmt(e)
            lrt(e) = l(2, 2, n, e) - lmt(e)
         end do
      else
         do e = 1, nelv
            do j = 2, n2
               l(1, j, 1, e) = lmr(e)
               l(n, j, 1, e) = lmr(e)
               l(j, 1, 1, e) = lms(e)
               l(j, n, 1, e) = lms(e)
            end do
         end do

         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_memcpy(l, this%swplen_d, this%dof%size(), &
                 HOST_TO_DEVICE, sync = .false.)
            call this%gs_h%op(l, this%dof%size(), GS_OP_ADD)
            call device_memcpy(l, this%swplen_d, this%dof%size(), &
                 DEVICE_TO_HOST, sync = .true.)
         else
            call this%gs_h%op(l, this%dof%size(), GS_OP_ADD)
         end if

         do e = 1, nelv
            llr(e) = l(1, 2, 1, e) - lmr(e)
            lrr(e) = l(n, 2, 1, e) - lmr(e)
            lls(e) = l(2, 1, 1, e) - lms(e)
            lrs(e) = l(2, n, 1, e) - lms(e)
         end do
      end if
    end associate
  end subroutine swap_lengths

  !> Here, spacing is based on harmonic mean.  pff 2/10/07
  !! We no longer base this on the finest grid, but rather
  !! the dofmap we are working with, Karp 210112
  subroutine plane_space(lr, ls, lt, i1, i2, w, x, y, z, &
       nx, nxn, nz0, nzn, nelv, gdim)
    integer, intent(in) :: nxn, nzn, i1, i2, nelv, gdim, nx, nz0
    real(kind=rp), intent(inout) :: lr(nelv), ls(nelv), lt(nelv)
    real(kind=rp), intent(inout) :: w(nx)
    real(kind=rp), intent(in) :: x(0:nxn, 0:nxn, nz0:nzn, nelv)
    real(kind=rp), intent(in) :: y(0:nxn, 0:nxn, nz0:nzn, nelv)
    real(kind=rp), intent(in) :: z(0:nxn, 0:nxn, nz0:nzn, nelv)
    real(kind=rp) :: lr2, ls2, lt2, weight, wsum
    integer :: ny, nz, j1, k1, j2, k2, i, j, k, ie
    ny = nx
    nz = nx
    j1 = i1
    k1 = i1
    j2 = i2
    k2 = i2
    !   Now, for each element, compute lr,ls,lt between specified planes
    do ie = 1, nelv
       if (gdim .eq. 3) then
          lr2 = 0d0
          wsum = 0d0
          do k = 1, nz
             do j = 1, ny
                weight = w(j)*w(k)
                lr2 = lr2 + weight /&
                     ( (x(i2, j, k, ie) - x(i1, j, k, ie))**2&
                     + (y(i2, j, k, ie) - y(i1, j, k, ie))**2&
                     + (z(i2, j, k, ie) - z(i1, j, k, ie))**2 )
                wsum = wsum + weight
             end do
          end do
          lr2 = lr2/wsum
          lr(ie) = 1d0/sqrt(lr2)
          ls2 = 0d0
          wsum = 0d0
          do k = 1, nz
             do i = 1, nx
                weight = w(i)*w(k)
                ls2 = ls2 + weight / &
                     ( (x(i, j2, k, ie) - x(i, j1, k, ie))**2 &
                     + (y(i, j2, k, ie) - y(i, j1, k, ie))**2 &
                     + (z(i, j2, k, ie) - z(i, j1, k, ie))**2 )
                wsum = wsum + weight
             end do
          end do
          ls2 = ls2/wsum
          ls(ie) = 1d0/sqrt(ls2)
          lt2 = 0d0
          wsum = 0d0
          do j = 1, ny
             do i = 1, nx
                weight = w(i)*w(j)
                lt2 = lt2 + weight / &
                     ( (x(i, j, k2, ie) - x(i, j, k1, ie))**2 &
                     + (y(i, j, k2, ie) - y(i, j, k1, ie))**2 &
                     + (z(i, j, k2, ie) - z(i, j, k1, ie))**2 )
                wsum = wsum + weight
             end do
          end do
          lt2 = lt2/wsum
          lt(ie) = 1d0/sqrt(lt2)
       else ! 2D
          lr2 = 0d0
          wsum = 0d0
          do j = 1, ny
             weight = w(j)
             lr2 = lr2 + weight / &
                  ( (x(i2, j, 1, ie) - x(i1, j, 1, ie))**2 &
                  + (y(i2, j, 1, ie) - y(i1, j, 1, ie))**2 )
             wsum = wsum + weight
          end do
          lr2 = lr2/wsum
          lr(ie) = 1d0/sqrt(lr2)
          ls2 = 0d0
          wsum = 0d0
          do i = 1, nx
             weight = w(i)
             ls2 = ls2 + weight / &
                  ( (x(i, j2, 1, ie) - x(i, j1, 1, ie))**2 &
                  + (y(i, j2, 1, ie) - y(i, j1, 1, ie))**2 )
             wsum = wsum + weight
          end do
          ls2 = ls2/wsum
          ls(ie) = 1d0/sqrt(ls2)
       end if
    end do
    ie = 1014
  end subroutine plane_space

  !> Setup the arrays s, d needed for the fast evaluation of the system
  subroutine fdm_setup_fast(this, ah, bh, nl, n)
    integer, intent(in) :: nl, n
    type(fdm_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: ah(n+1, n+1), bh(n+1)
    real(kind=rp), dimension(2*this%Xh%lx + 4) :: lr, ls, lt
    integer :: i, j, k
    integer :: ie, il, nr, ns, nt
    integer :: lbr, rbr, lbs, rbs, lbt, rbt
    real(kind=rp) :: eps, diag

    associate(s => this%s, d => this%d, &
         llr => this%len_lr, lls => this%len_ls, llt => this%len_lt, &
         lmr => this%len_mr, lms => this%len_ms, lmt => this%len_mt, &
         lrr => this%len_rr, lrs => this%len_rs, lrt => this%len_rt)
      do ie = 1, this%dof%msh%nelv
         lbr = this%dof%msh%facet_type(1, ie)
         rbr = this%dof%msh%facet_type(2, ie)
         lbs = this%dof%msh%facet_type(3, ie)
         rbs = this%dof%msh%facet_type(4, ie)
         lbt = this%dof%msh%facet_type(5, ie)
         rbt = this%dof%msh%facet_type(6, ie)

         nr = nl
         ns = nl
         nt = nl
         call fdm_setup_fast1d(s(1, 1, 1, ie), lr, nr, lbr, rbr, &
              llr(ie), lmr(ie), lrr(ie), ah, bh, n)
         call fdm_setup_fast1d(s(1, 1, 2, ie), ls, ns, lbs, rbs, &
              lls(ie), lms(ie), lrs(ie), ah, bh, n)
         if (this%dof%msh%gdim .eq. 3) then
            call fdm_setup_fast1d(s(1, 1, 3, ie), lt, nt, lbt, rbt, &
                 llt(ie), lmt(ie), lrt(ie), ah, bh, n)
         end if

         il = 1
         if (.not. this%dof%msh%gdim .eq. 3) then
            eps = 1d-5 * (vlmax(lr(2), nr - 2) + vlmax(ls(2), ns - 2))
            do j = 1, ns
               do i = 1, nr
                  diag = lr(i) + ls(j)
                  if (diag .gt. eps) then
                     d(il, ie) = 1.0_rp / diag
                  else
                     d(il, ie) = 0.0_rp
                  end if
                  il = il + 1
               end do
            end do
         else
            eps = 1d-5 * (vlmax(lr(2), nr - 2) + &
                 vlmax(ls(2), ns - 2) + vlmax(lt(2), nt - 2))
            do k = 1, nt
               do j = 1, ns
                  do i = 1, nr
                     diag = lr(i) + ls(j) + lt(k)
                     if (diag .gt. eps) then
                        d(il, ie) = 1.0_rp / diag
                     else
                        d(il, ie) = 0.0_rp
                     end if
                     il = il + 1
                  end do
               end do
            end do
         end if
      end do
    end associate

  end subroutine fdm_setup_fast

  subroutine fdm_do_semhat_weight(jgl, dgl, bgl, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: jgl(1:n-1, 0:n), dgl(1:n-1, 0:n)
    real(kind=rp), intent(in) :: bgl(1:n-1)
    integer :: i, j

    do j = 0, n
       do i = 1, n - 1
          jgl(i, j) = bgl(i) * jgl(i, j)
          dgl(i, j) = bgl(i) * dgl(i, j)
       end do
    end do
  end subroutine fdm_do_semhat_weight

  subroutine fdm_setup_fast_sem(this, bh, dgl, jgl, nl, n)
    integer, intent(in) :: nl, n
    type(fdm_t), intent(inout) :: this
    real(kind=rp), intent(in) :: bh(0:n)
    real(kind=rp), intent(in) :: dgl(1:n-1, 0:n), jgl(1:n-1, 0:n)
    real(kind=rp), dimension(2*this%Xh%lx + 4) :: lr, ls, lt
    integer :: i, j, k
    integer :: ie, il, nr, ns, nt
    integer :: lbr, rbr, lbs, rbs, lbt, rbt
    real(kind=rp) :: eps, diag

    associate(s => this%s, d => this%d, &
         llr => this%len_lr, lls => this%len_ls, llt => this%len_lt, &
         lmr => this%len_mr, lms => this%len_ms, lmt => this%len_mt, &
         lrr => this%len_rr, lrs => this%len_rs, lrt => this%len_rt)
      do ie = 1, this%dof%msh%nelv
         lbr = this%dof%msh%facet_type(1, ie)
         rbr = this%dof%msh%facet_type(2, ie)
         lbs = this%dof%msh%facet_type(3, ie)
         rbs = this%dof%msh%facet_type(4, ie)
         lbt = this%dof%msh%facet_type(5, ie)
         rbt = this%dof%msh%facet_type(6, ie)

         nr = nl
         ns = nl
         nt = nl
         call fdm_setup_fast1d_sem(s(1, 1, 1, ie), lr, nr, lbr, rbr, &
              llr(ie), lmr(ie), lrr(ie), bh, dgl, jgl, n)
         call fdm_setup_fast1d_sem(s(1, 1, 2, ie), ls, ns, lbs, rbs, &
              lls(ie), lms(ie), lrs(ie), bh, dgl, jgl, n)
         if (this%dof%msh%gdim .eq. 3) then
            call fdm_setup_fast1d_sem(s(1, 1, 3, ie), lt, nt, lbt, rbt, &
                 llt(ie), lmt(ie), lrt(ie), bh, dgl, jgl, n)
         end if

         il = 1
         if (.not. this%dof%msh%gdim .eq. 3) then
            eps = 1d-5 * (vlmax(lr(2), nr - 2) + vlmax(ls(2), ns - 2))
            do j = 1, ns
               do i = 1, nr
                  diag = lr(i) + ls(j)
                  if (diag .gt. eps) then
                     d(il, ie) = 1.0_rp / diag
                  else
                     d(il, ie) = 0.0_rp
                  end if
                  il = il + 1
               end do
            end do
         else
            eps = 1d-5 * (vlmax(lr(2), nr - 2) + &
                 vlmax(ls(2), ns - 2) + vlmax(lt(2), nt - 2))
            do k = 1, nt
               do j = 1, ns
                  do i = 1, nr
                     diag = lr(i) + ls(j) + lt(k)
                     if (diag .gt. eps) then
                        d(il, ie) = 1.0_rp / diag
                     else
                        d(il, ie) = 0.0_rp
                     end if
                     il = il + 1
                  end do
               end do
            end do
         end if
      end do
    end associate

  end subroutine fdm_setup_fast_sem

  subroutine fdm_setup_fast1d_sem(s, lam, nl, lbc, rbc, ll, lm, lr, bh, dgl, jgl, n)
    integer, intent(in) :: nl, lbc, rbc, n
    real(kind=rp), intent(inout) :: s(0:n, 0:n, 2), lam(nl), ll, lm, lr
    real(kind=rp), intent(in) :: bh(0:n)
    real(kind=rp), intent(in) :: dgl(1:n-1, 0:n), jgl(1:n-1, 0:n)
    integer :: bb0, bb1, eb0, eb1
    logical :: l, r
    real(kind=rp) :: b(0:n, 0:n)

    if (lbc .eq. 2 .or. lbc .eq. 3) then
       eb0 = 1
    else
       eb0 = 0
    end if
    if (rbc .eq. 2 .or. rbc .eq. 3) then
       eb1 = n - 1
    else
       eb1 = n
    end if
    if (lbc .eq. 2) then
       bb0 = 1
    else
       bb0 = 0
    end if
    if (rbc .eq. 2) then
       bb1 = n - 1
    else
       bb1 = n
    end if

    l = (lbc .eq. 0)
    r = (rbc .eq. 0)

    call fdm_setup_fast1d_sem_op(s(0, 0, 1), eb0, eb1, l, r, ll, lm, lr, bh, &
         dgl, 0, n)
    call fdm_setup_fast1d_sem_op(b, bb0, bb1, l, r, ll, lm, lr, bh, jgl, 1, n)

    call generalev(s(0, 0, 1), b, lam, nl, nl)
    if (.not. l) call row_zero(s(0, 0, 1), nl, nl, 1)
    if (.not. r) call row_zero(s(0, 0, 1), nl, nl, nl)
    call trsp(s(0, 0, 2), nl, s(0, 0, 1), nl)
  end subroutine fdm_setup_fast1d_sem

  subroutine fdm_setup_fast1d_sem_op(g, b0, b1, l, r, ll, lm, lr, bh, jgl, jscl, n)
    integer, intent(in) :: b0, b1, jscl, n
    logical, intent(in) :: l, r
    real(kind=rp), intent(inout) :: g(0:n, 0:n)
    real(kind=rp), intent(in) :: bh(0:n), jgl(1:n-1, 0:n)
    real(kind=rp), intent(in) :: ll, lm, lr
    real(kind=rp) :: bl(0:n), bm(0:n), br(0:n)
    real(kind=rp) :: gl, gm, gr, gll, glm, gmm, gmr, grr
    integer :: i, j, k

    if (jscl .eq. 0) then
       gl = 1.0_rp
       gm = 1.0_rp
       gr = 1.0_rp
    else
       gl = 0.5_rp * ll
       gm = 0.5_rp * lm
       gr = 0.5_rp * lr
    end if
    gll = gl * gl
    glm = gl * gm
    gmm = gm * gm
    gmr = gm * gr
    grr = gr * gr

    do i = 1, n - 1
       bm(i) = 2.0_rp / (lm * bh(i))
    end do
    if (b0 .eq. 0) then
       bm(0) = 0.5_rp * lm * bh(0)
       if (l) bm(0) = bm(0) + 0.5_rp * ll * bh(n)
       bm(0) = 1.0_rp / bm(0)
    end if
    if (b1 .eq. n) then
       bm(n) = 0.5_rp * lm * bh(n)
       if (r) bm(n) = bm(n) + 0.5_rp * lr * bh(0)
       bm(n) = 1.0_rp / bm(n)
    end if

    if (l) then
       do i = 0, n - 1
          bl(i) = 2.0_rp / (ll * bh(i))
       end do
       bl(n) = bm(0)
    end if
    if (r) then
       do i = 1, n
          br(i) = 2.0_rp / (lr * bh(i))
       end do
       br(0) = bm(n)
    end if

    call rzero(g, (n + 1) * (n + 1))
    do j = 1, n - 1
       do i = 1, n - 1
          do k = b0, b1
             g(i, j) = g(i, j) + gmm * jgl(i, k) * bm(k) * jgl(j, k)
          end do
       end do
    end do

    if (l) then
       do i = 1, n - 1
          g(i, 0) = glm * jgl(i, 0) * bm(0) * jgl(n - 1, n)
          g(0, i) = g(i, 0)
       end do
       do i = 0, n
          g(0, 0) = g(0, 0) + gll * jgl(n - 1, i) * bl(i) * jgl(n - 1, i)
       end do
    else
       g(0, 0) = 1.0_rp
    end if

    if (r) then
       do i = 1, n - 1
          g(i, n) = gmr * jgl(i, n) * bm(n) * jgl(1, 0)
          g(n, i) = g(i, n)
       end do
       do i = 0, n
          g(n, n) = g(n, n) + grr * jgl(1, i) * br(i) * jgl(1, i)
       end do
    else
       g(n, n) = 1.0_rp
    end if
  end subroutine fdm_setup_fast1d_sem_op

  subroutine fdm_setup_fast1d(s, lam, nl, lbc, rbc, ll, lm, lr, ah, bh, n)
    integer, intent(in) :: nl, lbc, rbc, n
    real(kind=rp), intent(inout) :: s(nl, nl, 2), lam(nl), ll, lm, lr
    real(kind=rp), intent(inout) :: ah(0:n, 0:n), bh(0:n)
    integer :: lx1, lxm
    real(kind=rp) :: b(2*(n+3)**2)

    lx1 = n + 1
    lxm = lx1 + 2

    call fdm_setup_fast1d_a(s, lbc, rbc, ll, lm, lr, ah, n)
    call fdm_setup_fast1d_b(b, lbc, rbc, ll, lm, lr, bh, n)
    call generalev(s, b, lam, nl, lx1)
    if (lbc .gt. 0) call row_zero(s, nl, nl, 1)
    if (lbc .eq. 1) call row_zero(s, nl, nl, 2)
    if (rbc .gt. 0) call row_zero(s, nl, nl, nl)
    if (rbc .eq. 1) call row_zero(s, nl, nl, nl-1)

    call trsp(s(1, 1, 2), nl, s, nl)

  end subroutine fdm_setup_fast1d

  !> Solve the generalized eigenvalue problem /$ A x = lam B x/$
  !! A -- symm.
  !! B -- symm., pos. definite
  subroutine generalev(a, b, lam, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: a(n, n), b(n, n), lam(n)
    integer :: lbw, lw
    real(kind=rp) :: bw(4*(lx+2)**3)

    lbw = 4*(lx+2)**3
    lw = n*n
    call sygv(a, b, lam, n, lx, bw, lbw)

  end subroutine generalev

  subroutine sp_sygv(a, b, lam, n, lx, bw, lbw)
    integer, intent(in) :: n, lx, lbw
    real(kind=sp), intent(inout) :: a(n, n), b(n, n), lam(n)
    real(kind=sp) :: bw(4*(lx+2)**3)
    integer :: info = 0
    call ssygv(1, 'V', 'U', n, a, n, b, n, lam, bw, lbw, info)
  end subroutine sp_sygv

  subroutine dp_sygv(a, b, lam, n, lx, bw, lbw)
    integer, intent(in) :: n, lx, lbw
    real(kind=dp), intent(inout) :: a(n, n), b(n, n), lam(n)
    real(kind=dp) :: bw(4*(lx+2)**3)
    integer :: info = 0
    call dsygv(1, 'V', 'U', n, a, n, b, n, lam, bw, lbw, info)
  end subroutine dp_sygv

  subroutine qp_sygv(a, b, lam, n, lx, bw, lbw)
    integer, intent(in) :: n, lx, lbw
    real(kind=qp), intent(inout) :: a(n, n), b(n, n), lam(n)
    real(kind=dp) :: a2(n, n), b2(n, n), lam2(n)
    real(kind=qp) :: bw(4*(lx+2)**3)
    real(kind=dp) :: bw2(4*(lx+2)**3)
    integer :: info = 0

    a2 = real(a, dp)
    b2 = real(b, dp)
    lam2 = real(lam, dp)
    call dsygv(1, 'V', 'U', n, a2, n, b2, n, lam2, bw2, lbw, info)
    a = real(a2, qp)
    b = real(b2, qp)
    lam = real(lam2, qp)
    if (pe_rank .eq. 0) then
       call neko_warning('Real precision choice not supported for fdm, ' // &
            'treating it as double')
    end if

  end subroutine qp_sygv

  subroutine fdm_setup_fast1d_a(a, lbc, rbc, ll, lm, lr, ah, n)
    integer, intent(in) ::lbc, rbc, n
    real(kind=rp), intent(inout) :: a(0:n+2, 0:n+2), ll, lm, lr
    real(kind=rp), intent(inout) :: ah(0:n, 0:n)
    real(kind=rp) :: fac
    integer :: i, j, i0, i1

    i0 = 0
    if (lbc .eq. 1) i0 = 1
    i1 = n
    if (rbc .eq. 1) i1 = n - 1

    call rzero(a, (n+3) * (n+3))

    fac = 2.0_rp / lm
    a(1, 1) = 1.0_rp
    a(n+1, n+1) = 1.0_rp

    do j = i0, i1
       do i = i0, i1
          a(i+1, j+1) = fac * ah(i, j)
       end do
    end do

    if (lbc .eq. 0) then
       fac = 2.0_rp / ll
       a(0, 0) = fac * ah(n-1, n-1)
       a(1, 0) = fac * ah(n, n-1)
       a(0, 1) = fac * ah(n-1, n)
       a(1, 1) = a(1, 1) + fac * ah(n, n)
    else
       a(0, 0) = 1.0_rp
    end if

    if (rbc .eq. 0) then
       fac = 2.0_rp / lr
       a(n+1, n+1) = a(n+1, n+1) + fac * ah(0, 0)
       a(n+2, n+1) = fac * ah(1, 0)
       a(n+1, n+2) = fac * ah(0, 1)
       a(n+2, n+2) = fac * ah(1, 1)
    else
       a(n+2, n+2) = 1.0_rp
    end if

  end subroutine fdm_setup_fast1d_a

  subroutine fdm_setup_fast1d_b(b, lbc, rbc, ll, lm, lr, bh, n)
    integer, intent(in) :: lbc, rbc, n
    real(kind=rp), intent(inout) :: b(0:n+2, 0:n+2), ll, lm, lr
    real(kind=rp), intent(inout) :: bh(0:n)
    real(kind=rp) :: fac
    integer :: i, i0, i1

    i0 = 0
    if (lbc .eq. 1) i0 = 1
    i1 = n
    if (rbc .eq. 1) i1 = n - 1

    call rzero(b, (n + 3) * (n + 3))

    fac = 0.5_rp * lm
    b(1, 1) = 1.0_rp
    b(n+1, n+1) = 1.0_rp

    do i = i0, i1
       b(i+1, i+1) = fac * bh(i)
    end do

    if (lbc .eq. 0) then
       fac = 0.5_rp * ll
       b(0, 0) = fac * bh(n-1)
       b(1, 1) = b(1, 1) + fac * bh(n)
    else
       b(0, 0) = 1.0_rp
    end if

    if (rbc .eq. 0) then
       fac = 0.5_rp * lr
       b(n+1, n+1) = b(n+1, n+1) + fac * bh(0)
       b(n+2, n+2) = fac * bh(1)
    else
       b(n+2, n+2) = 1.0_rp
    end if

  end subroutine fdm_setup_fast1d_b

  subroutine fdm_free(this)
    class(fdm_t), intent(inout) :: this

    if (allocated(this%s)) then
       deallocate(this%s)
    end if

    if (allocated(this%d)) then
       deallocate(this%d)
    end if

    if (allocated(this%len_lr)) then
       deallocate(this%len_lr)
    end if

    if (allocated(this%len_ls)) then
       deallocate(this%len_ls)
    end if

    if (allocated(this%len_lt)) then
       deallocate(this%len_lt)
    end if

    if (allocated(this%len_mr)) then
       deallocate(this%len_mr)
    end if

    if (allocated(this%len_ms)) then
       deallocate(this%len_ms)
    end if

    if (allocated(this%len_mt)) then
       deallocate(this%len_mt)
    end if

    if (allocated(this%len_rr)) then
       deallocate(this%len_rr)
    end if

    if (allocated(this%len_rs)) then
       deallocate(this%len_rs)
    end if

    if (allocated(this%len_rt)) then
       deallocate(this%len_rt)
    end if

    if (allocated(this%swplen)) then
       deallocate(this%swplen)
    end if

    this%nl = 0
    nullify(this%Xh)
    nullify(this%dof)
    nullify(this%gs_h)
    nullify(this%msh)

    if (c_associated(this%s_d)) then
       call device_free(this%s_d)
    end if
    if (c_associated(this%d_d)) then
       call device_free(this%d_d)
    end if
    if (c_associated(this%swplen_d)) then
       call device_free(this%swplen_d)
    end if

  end subroutine fdm_free

  subroutine fdm_compute(this, e, r, stream)
    class(fdm_t), intent(inout) :: this
    real(kind=rp), dimension(this%nl**this%msh%gdim, this%msh%nelv), &
         intent(inout) :: e, r
    type(c_ptr), optional :: stream
    type(c_ptr) :: strm

    if (present(stream)) then
       strm = stream
    else
       strm = glb_cmd_queue
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       call fdm_do_fast_sx(e, r, this%s, this%d, &
            this%nl, this%msh%gdim, this%msh%nelv)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call fdm_do_fast_xsmm(e, r, this%s, this%d, &
            this%nl, this%msh%gdim, this%msh%nelv)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call fdm_do_fast_device(e, r, this%s, this%d, &
            this%nl, this%msh%gdim, this%msh%nelv, strm)
    else
       call fdm_do_fast_cpu(e, r, this%s, this%d, &
            this%nl, this%msh%gdim, this%msh%nelv)
    end if

  end subroutine fdm_compute


end module fdm
