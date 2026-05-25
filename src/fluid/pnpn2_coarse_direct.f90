! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   * Neither the name of the authors nor the names of its contributors may be
!     used to endorse or promote products derived from this software without
!     specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!> Explicit unique-dof direct coarse solve for the Pn/Pn-2 HSMG preconditioner.
module pnpn2_coarse_direct
  use ax_product, only : ax_t, ax_helm_factory
  use bc_list, only : bc_list_t
  use coefs, only : coef_t
  use comm, only : NEKO_COMM
  use dofmap, only : dofmap_t
  use field, only : field_t
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX, MPI_SUM
  use num_types, only : rp
  use space, only : space_t
  use utils, only : neko_error
  implicit none
  private

  type, public :: pnpn2_coarse_direct_t
     integer :: n_local = 0
     integer :: n_global = 0
     integer :: max_gid = 0
     integer :: lx_fine = 0
     integer :: lz_fine = 0
     integer :: nloc_crs = 0
     integer, allocatable :: gid_to_idx(:)
     integer, allocatable :: local_to_idx(:)
     integer, allocatable :: ipiv(:)
     real(kind=rp), allocatable :: lu(:,:)
     real(kind=rp), allocatable :: rhs(:)
     real(kind=rp), allocatable :: sol(:)
     real(kind=rp), allocatable :: basis(:,:,:,:)
     logical :: initialized = .false.
     type(field_t) :: u
     type(field_t) :: w
     class(ax_t), allocatable :: ax
   contains
     procedure, pass(this) :: init => pnpn2_coarse_direct_init
     procedure, pass(this) :: free => pnpn2_coarse_direct_free
     procedure, pass(this) :: solve => pnpn2_coarse_direct_solve
  end type pnpn2_coarse_direct_t

contains

  subroutine pnpn2_coarse_direct_init(this, coef_fine, dm_fine, dm_crs, bclst)
    class(pnpn2_coarse_direct_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef_fine
    type(dofmap_t), intent(inout), target :: dm_fine
    type(dofmap_t), intent(inout), target :: dm_crs
    type(bc_list_t), intent(inout) :: bclst
    integer :: ierr, i, e, gid, idx, iloc, jloc
    integer :: local_max_gid
    integer, allocatable :: present_local(:), present_global(:)
    real(kind=rp), allocatable :: amat(:,:), amat_global(:,:), local_amat(:,:,:)
    real(kind=rp), allocatable :: h1_save(:), h2_save(:)
    logical :: ifh2_save
    integer :: ncoef

    if (this%initialized) call this%free()

    this%n_local = dm_crs%size()
    this%lx_fine = coef_fine%Xh%lx
    this%lz_fine = coef_fine%Xh%lz
    this%nloc_crs = dm_crs%Xh%lx * dm_crs%Xh%ly * dm_crs%Xh%lz
    local_max_gid = 0
    do i = 1, this%n_local
      local_max_gid = max(local_max_gid, int(dm_crs%dof(i,1,1,1)))
    end do
    call MPI_Allreduce(local_max_gid, this%max_gid, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    if (this%max_gid .le. 0) then
      call neko_error('PnPn-2 direct coarse solver found no coarse dofs.')
    end if

    allocate(present_local(this%max_gid))
    allocate(present_global(this%max_gid))
    allocate(this%gid_to_idx(this%max_gid))
    allocate(this%local_to_idx(this%n_local))
    present_local = 0
    this%gid_to_idx = 0
    do i = 1, this%n_local
      gid = int(dm_crs%dof(i,1,1,1))
      if (gid .gt. 0) present_local(gid) = 1
    end do
    call MPI_Allreduce(present_local, present_global, this%max_gid, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)

    this%n_global = 0
    do gid = 1, this%max_gid
      if (present_global(gid) .ne. 0) then
        this%n_global = this%n_global + 1
        this%gid_to_idx(gid) = this%n_global
      end if
    end do

    do i = 1, this%n_local
      gid = int(dm_crs%dof(i,1,1,1))
      this%local_to_idx(i) = this%gid_to_idx(gid)
    end do

    allocate(amat(this%n_global, this%n_global))
    allocate(amat_global(this%n_global, this%n_global))
    allocate(local_amat(this%nloc_crs, this%nloc_crs, coef_fine%msh%nelv))
    allocate(this%lu(this%n_global, this%n_global))
    allocate(this%rhs(this%n_global))
    allocate(this%sol(this%n_global))
    allocate(this%ipiv(this%n_global))
    allocate(this%basis(this%lx_fine, this%lx_fine, this%lz_fine, this%nloc_crs))
    call this%u%init(dm_fine, 'pnpn2_coarse_direct_u')
    call this%w%init(dm_fine, 'pnpn2_coarse_direct_w')
    call ax_helm_factory(this%ax, full_formulation = .false.)
    call generate_crs_basis(this%basis, coef_fine%Xh)
    ncoef = coef_fine%dof%size()
    allocate(h1_save(ncoef), h2_save(ncoef))
    h1_save = coef_fine%h1(1:ncoef,1,1,1)
    h2_save = coef_fine%h2(1:ncoef,1,1,1)
    ifh2_save = coef_fine%ifh2
    coef_fine%h1(1:ncoef,1,1,1) = 1.0_rp
    coef_fine%h2(1:ncoef,1,1,1) = 0.0_rp
    coef_fine%ifh2 = .false.

    amat = 0.0_rp
    local_amat = 0.0_rp
    do jloc = 1, this%nloc_crs
      this%u%x = 0.0_rp
      this%w%x = 0.0_rp
      do e = 1, coef_fine%msh%nelv
        this%u%x(:,:,:,e) = this%basis(:,:,:,jloc)
      end do
      call this%ax%compute(this%w%x, this%u%x, coef_fine, this%u%msh, this%u%Xh)

      do e = 1, coef_fine%msh%nelv
        do iloc = 1, this%nloc_crs
          local_amat(iloc,jloc,e) = sum(this%basis(:,:,:,iloc) * this%w%x(:,:,:,e))
        end do
      end do


    end do
    do e = 1, coef_fine%msh%nelv
      do jloc = 1, this%nloc_crs
        do iloc = 1, this%nloc_crs
          gid = this%local_to_idx(local_crs_id(dm_crs%Xh, iloc, e))
          idx = this%local_to_idx(local_crs_id(dm_crs%Xh, jloc, e))
          amat(gid, idx) = amat(gid, idx) + local_amat(iloc, jloc, e)
        end do
      end do
    end do
    coef_fine%h1(1:ncoef,1,1,1) = h1_save
    coef_fine%h2(1:ncoef,1,1,1) = h2_save
    coef_fine%ifh2 = ifh2_save

    call MPI_Allreduce(amat, amat_global, this%n_global * this%n_global, &
         MPI_REAL_RP(), MPI_SUM, NEKO_COMM, ierr)

    ! The pressure coarse operator has a constant nullspace for fully periodic
    ! cases. Pinning one coarse dof gives one representative; the outer pressure
    ! Krylov solve orthogonalizes the preconditioned vector afterwards.
    amat_global(1,:) = 0.0_rp
    amat_global(:,1) = 0.0_rp
    amat_global(1,1) = 1.0_rp

    this%lu = amat_global
    call dense_lu_factor(this%lu, this%ipiv, this%n_global)
    this%initialized = .true.

    deallocate(amat)
    deallocate(amat_global)
    deallocate(local_amat)
    deallocate(h1_save)
    deallocate(h2_save)
    deallocate(present_local)
    deallocate(present_global)
  end subroutine pnpn2_coarse_direct_init

  subroutine pnpn2_coarse_direct_free(this)
    class(pnpn2_coarse_direct_t), intent(inout) :: this

    if (.not. this%initialized) return

    if (allocated(this%gid_to_idx)) deallocate(this%gid_to_idx)
    if (allocated(this%local_to_idx)) deallocate(this%local_to_idx)
    if (allocated(this%ipiv)) deallocate(this%ipiv)
    if (allocated(this%lu)) deallocate(this%lu)
    if (allocated(this%rhs)) deallocate(this%rhs)
    if (allocated(this%sol)) deallocate(this%sol)
    if (allocated(this%basis)) deallocate(this%basis)
    call this%u%free()
    call this%w%free()
    if (allocated(this%ax)) deallocate(this%ax)
    this%n_local = 0
    this%n_global = 0
    this%max_gid = 0
    this%lx_fine = 0
    this%lz_fine = 0
    this%nloc_crs = 0
    this%initialized = .false.
  end subroutine pnpn2_coarse_direct_free

  subroutine pnpn2_coarse_direct_solve(this, e, r, coef)
   class(pnpn2_coarse_direct_t), intent(inout) :: this
   real(kind=rp), intent(inout) :: e(this%n_local)
   real(kind=rp), intent(in) :: r(this%n_local)
   type(coef_t), intent(in) :: coef
   integer :: ierr, i, idx

    if (.not. allocated(this%lu)) then
      call neko_error('PnPn-2 direct coarse solver used before initialization.')
    end if

    this%rhs = 0.0_rp
    do i = 1, this%n_local
      idx = this%local_to_idx(i)
      this%rhs(idx) = this%rhs(idx) + r(i)
    end do

    call MPI_Allreduce(this%rhs, this%sol, this%n_global, MPI_REAL_RP(), &
         MPI_SUM, NEKO_COMM, ierr)

    this%sol(1) = 0.0_rp
    call dense_lu_solve(this%lu, this%ipiv, this%sol, this%n_global)

    do i = 1, this%n_local
      e(i) = this%sol(this%local_to_idx(i))
    end do
  end subroutine pnpn2_coarse_direct_solve

  subroutine dense_lu_factor(a, ipiv, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(n,n)
    integer, intent(out) :: ipiv(n)
    integer :: k, p, i
    real(kind=rp) :: pivot, tmp

    do k = 1, n
      p = k
      pivot = abs(a(k,k))
      do i = k + 1, n
        if (abs(a(i,k)) .gt. pivot) then
          pivot = abs(a(i,k))
          p = i
        end if
      end do
      if (pivot .le. tiny(1.0_rp)) then
        call neko_error('Singular PnPn-2 direct coarse matrix.')
      end if
      ipiv(k) = p
      if (p .ne. k) then
        do i = 1, n
          tmp = a(k,i)
          a(k,i) = a(p,i)
          a(p,i) = tmp
        end do
      end if
      do i = k + 1, n
        a(i,k) = a(i,k) / a(k,k)
        a(i,k + 1:n) = a(i,k + 1:n) - a(i,k) * a(k,k + 1:n)
      end do
    end do
  end subroutine dense_lu_factor

  subroutine dense_lu_solve(a, ipiv, b, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: a(n,n)
    integer, intent(in) :: ipiv(n)
    real(kind=rp), intent(inout) :: b(n)
    integer :: i
    real(kind=rp) :: tmp

    do i = 1, n
      if (ipiv(i) .ne. i) then
        tmp = b(i)
        b(i) = b(ipiv(i))
        b(ipiv(i)) = tmp
      end if
    end do

    do i = 2, n
      b(i) = b(i) - dot_product(a(i,1:i-1), b(1:i-1))
    end do

    do i = n, 1, -1
      if (i .lt. n) then
        b(i) = b(i) - dot_product(a(i,i+1:n), b(i+1:n))
      end if
      b(i) = b(i) / a(i,i)
    end do
  end subroutine dense_lu_solve

  function MPI_REAL_RP() result(dtype)
    use comm, only : MPI_REAL_PRECISION
    use mpi_f08, only : MPI_Datatype
    type(MPI_Datatype) :: dtype

    dtype = MPI_REAL_PRECISION
  end function MPI_REAL_RP

  integer function local_crs_id(Xh_crs, iloc, e) result(id)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc, e
    integer :: i, j, k

    i = local_crs_i(Xh_crs, iloc)
    j = local_crs_j(Xh_crs, iloc)
    k = local_crs_k(Xh_crs, iloc)
    id = i + (j - 1) * Xh_crs%lx + (k - 1) * Xh_crs%lx * Xh_crs%ly + &
         (e - 1) * Xh_crs%lx * Xh_crs%ly * Xh_crs%lz
  end function local_crs_id

  integer function local_crs_i(Xh_crs, iloc) result(i)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    i = mod(iloc - 1, Xh_crs%lx) + 1
  end function local_crs_i

  integer function local_crs_j(Xh_crs, iloc) result(j)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    j = mod((iloc - 1) / Xh_crs%lx, Xh_crs%ly) + 1
  end function local_crs_j

  integer function local_crs_k(Xh_crs, iloc) result(k)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    k = ((iloc - 1) / (Xh_crs%lx * Xh_crs%ly)) + 1
  end function local_crs_k

  subroutine generate_crs_basis(basis, Xh_space)
    type(space_t), intent(in) :: Xh_space
    real(kind=rp), intent(inout) :: basis(:,:,:,:)
    real(kind=rp), allocatable :: z0(:), z1(:), zr(:), zs(:), zt(:)
    integer :: i, j, p, q, kk

    allocate(z0(size(basis,1)), z1(size(basis,1)))
    allocate(zr(size(basis,1)), zs(size(basis,1)), zt(size(basis,1)))

    do i = 1, size(basis,1)
      z0(i) = 0.5_rp * (1.0_rp - Xh_space%zg(i,1))
      z1(i) = 0.5_rp * (1.0_rp + Xh_space%zg(i,1))
    end do

    basis = 0.0_rp
    do j = 1, size(basis, 4)
      zr = z0
      zs = z0
      zt = z0

      if (mod(j, 2) .eq. 0) zr = z1
      if (j .eq. 3 .or. j .eq. 4 .or. j .eq. 7 .or. j .eq. 8) zs = z1
      if (j .gt. 4) zt = z1

      if (size(basis,3) .gt. 1) then
        do kk = 1, size(basis,3)
          do q = 1, size(basis,2)
            do p = 1, size(basis,1)
              basis(p,q,kk,j) = zr(p) * zs(q) * zt(kk)
            end do
          end do
        end do
      else
        do q = 1, size(basis,2)
          do p = 1, size(basis,1)
            basis(p,q,1,j) = zr(p) * zs(q)
          end do
        end do
      end if
    end do

    deallocate(z0, z1, zr, zs, zt)
  end subroutine generate_crs_basis

end module pnpn2_coarse_direct
