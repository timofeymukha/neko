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
  use gather_scatter, only : gs_t, GS_OP_ADD
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX, MPI_SUM
  use num_types, only : rp
  use utils, only : neko_error
  implicit none
  private

  type, public :: pnpn2_coarse_direct_t
     integer :: n_local = 0
     integer :: n_global = 0
     integer :: max_gid = 0
     integer, allocatable :: gid_to_idx(:)
     integer, allocatable :: local_to_idx(:)
     integer, allocatable :: ipiv(:)
     real(kind=rp), allocatable :: lu(:,:)
     real(kind=rp), allocatable :: rhs(:)
     real(kind=rp), allocatable :: sol(:)
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

  subroutine pnpn2_coarse_direct_init(this, coef, dm, gs, bclst)
    class(pnpn2_coarse_direct_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dm
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bclst
    integer :: ierr, i, j, gid, idx
    integer :: local_max_gid
    integer, allocatable :: present_local(:), present_global(:)
    real(kind=rp), allocatable :: amat(:,:), amat_global(:,:)

    if (this%initialized) call this%free()

    this%n_local = dm%size()
    local_max_gid = 0
    do i = 1, this%n_local
      local_max_gid = max(local_max_gid, int(dm%dof(i,1,1,1)))
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
      gid = int(dm%dof(i,1,1,1))
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
      gid = int(dm%dof(i,1,1,1))
      this%local_to_idx(i) = this%gid_to_idx(gid)
    end do

    allocate(amat(this%n_global, this%n_global))
    allocate(amat_global(this%n_global, this%n_global))
    allocate(this%lu(this%n_global, this%n_global))
    allocate(this%rhs(this%n_global))
    allocate(this%sol(this%n_global))
    allocate(this%ipiv(this%n_global))

    call this%u%init(dm, 'pnpn2_coarse_direct_u')
    call this%w%init(dm, 'pnpn2_coarse_direct_w')
    call ax_helm_factory(this%ax, full_formulation = .false.)

    amat = 0.0_rp
    do j = 1, this%n_global
      this%u%x = 0.0_rp
      this%w%x = 0.0_rp
      do i = 1, this%n_local
        if (this%local_to_idx(i) .eq. j) then
          this%u%x(i,1,1,1) = 1.0_rp
        end if
      end do

      call this%ax%compute(this%w%x, this%u%x, coef, this%u%msh, this%u%Xh)
      call gs%op(this%w%x, this%n_local, GS_OP_ADD)
      call bclst%apply_scalar(this%w%x, this%n_local)

      do i = 1, this%n_local
        idx = this%local_to_idx(i)
        amat(idx,j) = amat(idx,j) + this%w%x(i,1,1,1) * coef%mult(i,1,1,1)
      end do
    end do

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
    call this%u%free()
    call this%w%free()
    if (allocated(this%ax)) deallocate(this%ax)
    this%n_local = 0
    this%n_global = 0
    this%max_gid = 0
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
      this%rhs(idx) = this%rhs(idx) + r(i) * coef%mult(i,1,1,1)
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

end module pnpn2_coarse_direct
