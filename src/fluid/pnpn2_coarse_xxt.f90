! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
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
!> Direct Nek5000 XXT coarse solver for the Pn/Pn-2 HSMG path.
!!
!! This is a Fortran port of the local algebra and separator handling in
!! `crs_xxt.c`.  The setup consumes the rank-local coarse-grid user dofs and the
!! rank-local Galerkin triplets, discovers the same separator topology that Nek
!! uses, factors the local `A_ll` block, builds the separator basis `X`, and
!! then applies the same `crs_solve` branches in the coarse solve.
module pnpn2_coarse_xxt
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only : MPI_Allgather, MPI_Allreduce, MPI_Send, MPI_Recv, &
       MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE, MPI_Datatype, MPI_MAX, &
       MPI_SUM
  use num_types, only : i8, rp
  use utils, only : neko_error
  implicit none
  private

  type, public :: xxt_csr_mat_t
    integer :: n = 0
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: col_ind(:)
    real(kind=rp), allocatable :: val(:)
   contains
    procedure, pass(this) :: free => xxt_csr_mat_free
  end type xxt_csr_mat_t

  type, public :: xxt_sparse_cholesky_t
    integer :: n = 0
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: col_ind(:)
    real(kind=rp), allocatable :: l_val(:)
    real(kind=rp), allocatable :: d(:)
   contains
    procedure, pass(this) :: free => xxt_sparse_cholesky_free
    procedure, pass(this) :: factor => xxt_sparse_cholesky_factor
    procedure, pass(this) :: solve => xxt_sparse_cholesky_solve
  end type xxt_sparse_cholesky_t

  type, public :: pnpn2_coarse_xxt_t
    integer :: pcoord = 0
    integer :: nsep = 0
    integer :: plevels = 0
    integer, allocatable :: pother(:)
    logical :: null_space = .false.
    integer :: un = 0
    integer :: cn = 0
    integer :: ln = 0
    integer :: sn = 0
    integer :: xn = 0
    integer :: ln_fac = 0
    integer, allocatable :: sep_size(:)
    integer, allocatable :: perm_u2c(:)
    real(kind=rp), allocatable :: share_weight(:)
    type(xxt_sparse_cholesky_t) :: fac_a_ll
    type(xxt_csr_mat_t) :: a_sl
    integer, allocatable :: xp(:)
    real(kind=rp), allocatable :: x(:)
    real(kind=rp), allocatable :: vl(:), vc(:), vx(:), combuf(:)
    logical :: tree_initialized = .false.
   contains
    procedure, pass(this) :: init_tree => pnpn2_coarse_xxt_init_tree
    procedure, pass(this) :: setup => pnpn2_coarse_xxt_setup
    procedure, pass(this) :: solve => pnpn2_coarse_xxt_solve
    procedure, pass(this) :: free => pnpn2_coarse_xxt_free
  end type pnpn2_coarse_xxt_t

  type :: xxt_dof_t
    integer(kind=i8) :: id = 0_i8
    integer :: level = 0
    integer :: count = 0
  end type xxt_dof_t

  type :: xxt_yale_mat_t
    integer :: i = 0
    integer :: j = 0
    real(kind=rp) :: v = 0.0_rp
  end type xxt_yale_mat_t

contains

  subroutine xxt_csr_mat_free(this)
    class(xxt_csr_mat_t), intent(inout) :: this
    if (allocated(this%row_ptr)) deallocate(this%row_ptr)
    if (allocated(this%col_ind)) deallocate(this%col_ind)
    if (allocated(this%val)) deallocate(this%val)
    this%n = 0
  end subroutine xxt_csr_mat_free

  subroutine xxt_sparse_cholesky_free(this)
    class(xxt_sparse_cholesky_t), intent(inout) :: this
    if (allocated(this%row_ptr)) deallocate(this%row_ptr)
    if (allocated(this%col_ind)) deallocate(this%col_ind)
    if (allocated(this%l_val)) deallocate(this%l_val)
    if (allocated(this%d)) deallocate(this%d)
    this%n = 0
  end subroutine xxt_sparse_cholesky_free

  !> Sparse LDL^T solve: x = A^{-1} b.
  !!
  !! Executes the three-sweep solve matching Nek5000's
  !! \c sparse_cholesky_solve:
  !!   1. Forward sweep  : (I - L) x = b
  !!   2. Diagonal scale : x = D^{-1} x
  !!   3. Backward sweep : (I - L)^T x = x
  subroutine xxt_sparse_cholesky_solve(this, x, b)
    class(xxt_sparse_cholesky_t), intent(in) :: this
    real(kind=rp), intent(inout) :: x(:)
    real(kind=rp), intent(in) :: b(:)
    integer :: i, p, n
    real(kind=rp) :: xi

    n = this%n
    if (n .eq. 0) return

    ! Forward sweep: x(i) = b(i) + sum_{j in row i of L} L(p)*x(j)
    do i = 1, n
      xi = b(i)
      do p = this%row_ptr(i), this%row_ptr(i + 1) - 1
        xi = xi + this%l_val(p) * x(this%col_ind(p))
      end do
      x(i) = xi
    end do

    ! Diagonal scaling: x(i) *= D(i)  [D stores D^{-1} from factorization]
    do concurrent (i = 1:n)
      x(i) = x(i) * this%d(i)
    end do

    ! Backward sweep: for i from n to 1, x(j) += L(p)*x(i) for j in row i of L
    do i = n, 1, -1
      xi = x(i)
      do p = this%row_ptr(i), this%row_ptr(i + 1) - 1
        x(this%col_ind(p)) = x(this%col_ind(p)) + this%l_val(p) * xi
      end do
    end do
  end subroutine xxt_sparse_cholesky_solve

  !> Sparse LDL^T factorization: A = (I-L) D^{-1} (I-L)^T.
  !!
  !! Direct Fortran port of Nek5000's \c factor_symbolic + \c factor_numeric
  !! from \c crs_xxt.c.  The input CSR may contain both upper and lower
  !! triangles; only strictly lower-triangle entries (j < i) are used.
  !!
  !! On output:
  !!   - \c this%row_ptr / \c this%col_ind  hold the sparsity pattern of L
  !!   - \c this%l_val  holds the non-unit entries of L
  !!   - \c this%d(i)   holds 1/a_{ii}^{(i)} (the inverse pivot used in solve)
  subroutine xxt_sparse_cholesky_factor(this, row_ptr, col_ind, val)
    class(xxt_sparse_cholesky_t), intent(inout) :: this
    integer, intent(in) :: row_ptr(:), col_ind(:)
    real(kind=rp), intent(in) :: val(:)
    integer, allocatable :: visit(:), parent(:)
    real(kind=rp), allocatable :: y_work(:)
    integer :: n, i, j, k, p, q, count, nz
    real(kind=rp) :: a_diag, yj, lij

    call this%free()

    n = size(row_ptr) - 1
    this%n = n
    if (n .eq. 0) return

    ! ---------------------------------------------------------------
    ! Symbolic phase: build the sparsity pattern of L
    ! (Nek5000's factor_symbolic)
    ! ---------------------------------------------------------------
    allocate(visit(n), parent(n))
    visit = 0

    ! --- First pass: count total fill-in entries (nz) ---
    nz = 0
    do i = 1, n
      visit(i) = i
      parent(i) = n + 1        ! n+1 = "no parent" sentinel
      do p = row_ptr(i), row_ptr(i + 1) - 1
        j = col_ind(p)
        if (j .ge. i) exit     ! only strictly lower triangle
        ! Walk elimination tree from j upward
        do while (visit(j) .ne. i)
          nz = nz + 1
          visit(j) = i
          if (parent(j) .eq. n + 1) then
            parent(j) = i
            exit
          end if
          j = parent(j)
        end do
      end do
    end do

    allocate(this%row_ptr(n + 1))
    allocate(this%col_ind(nz))

    ! --- Second pass: fill col_ind and row_ptr ---
    ! (visit is now in the same state as after the first pass; we
    !  reset visit(i)=i per row and re-run the identical traversal)
    this%row_ptr(1) = 1
    do i = 1, n
      visit(i) = i
      count = 0
      do p = row_ptr(i), row_ptr(i + 1) - 1
        j = col_ind(p)
        if (j .ge. i) exit     ! only strictly lower triangle
        do while (visit(j) .ne. i)
          count = count + 1
          this%col_ind(this%row_ptr(i) + count - 1) = j
          visit(j) = i
          j = parent(j)
        end do
      end do
      ! Sort column indices for this row
      if (count .gt. 1) then
        call sort_int(this%col_ind, this%row_ptr(i), this%row_ptr(i) + count - 1)
      end if
      this%row_ptr(i + 1) = this%row_ptr(i) + count
    end do

    deallocate(visit, parent)

    ! ---------------------------------------------------------------
    ! Numeric phase: compute L entries and diagonal D^{-1}
    ! (Nek5000's factor_numeric)
    ! ---------------------------------------------------------------
    allocate(this%d(n))
    allocate(this%l_val(nz))
    allocate(y_work(n))
    allocate(visit(n))

    this%l_val = 0.0_rp
    y_work = 0.0_rp
    visit = 0                  ! 0 = "not active for the current row"

    do i = 1, n
      visit(i) = 0             ! sentinel: diagonal, not in active column set
      a_diag = 0.0_rp

      ! Zero y_work for L's row i and mark those columns active
      do p = this%row_ptr(i), this%row_ptr(i + 1) - 1
        j = this%col_ind(p)
        y_work(j) = 0.0_rp
        visit(j) = i
      end do

      ! Load strictly lower triangle of A into y_work; grab diagonal
      do p = row_ptr(i), row_ptr(i + 1) - 1
        j = col_ind(p)
        if (j .ge. i) then
          if (j .eq. i) a_diag = val(p)
          exit
        end if
        y_work(j) = -val(p)
      end do

      ! Compute L entries for row i using previous rows
      do p = this%row_ptr(i), this%row_ptr(i + 1) - 1
        j = this%col_ind(p)
        yj = y_work(j)
        ! Scatter contributions from L's row j
        do q = this%row_ptr(j), this%row_ptr(j + 1) - 1
          k = this%col_ind(q)
          if (visit(k) .eq. i) yj = yj + this%l_val(q) * y_work(k)
        end do
        y_work(j) = yj
        lij = this%d(j) * yj
        this%l_val(p) = lij
        a_diag = a_diag - yj * lij
      end do

      this%d(i) = 1.0_rp / a_diag
    end do

    deallocate(y_work, visit)
  end subroutine xxt_sparse_cholesky_factor

  subroutine pnpn2_coarse_xxt_init_tree(this, rank, nproc, null_space)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    integer, intent(in) :: rank, nproc
    logical, intent(in) :: null_space
    integer :: n, c, odd, base, level, targ

    if (allocated(this%pother)) deallocate(this%pother)
    this%pcoord = 0
    this%nsep = 0
    this%plevels = 0
    this%tree_initialized = .false.

    this%null_space = null_space
    n = nproc
    c = 1
    odd = 0
    base = 0
    level = 0

    do while (n .gt. 1)
      level = level + 1
      odd = ishft(odd, 1) + iand(n, 1)
      c = ishft(c, 1)
      n = ishft(n, -1)
      if (rank .ge. base + n) then
        c = ior(c, 1)
        base = base + n
        n = n + iand(odd, 1)
      end if
    end do

    this%pcoord = c
    this%nsep = level + 1
    this%plevels = this%nsep - 1
    if (this%plevels .gt. 0) allocate(this%pother(this%plevels))

    n = nproc
    c = this%pcoord
    odd = 0
    do while (n .gt. 1)
      odd = ishft(odd, 1) + iand(n, 1)
      n = ishft(n, -1)
    end do

    n = 1
    odd = 0
    base = 0
    do level = 1, this%plevels
      if (iand(c, 1) .eq. 1) then
        targ = rank - (n - iand(odd, 1))
        this%pother(level) = -(targ + 1)
        this%plevels = level
        exit
      else
        this%pother(level) = rank + n
        c = ishft(c, -1)
        n = ishft(n, 1) + iand(odd, 1)
        odd = ishft(odd, -1)
      end if
    end do

    this%tree_initialized = .true.
  end subroutine pnpn2_coarse_xxt_init_tree

  subroutine pnpn2_coarse_xxt_free(this)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this

    if (allocated(this%pother)) deallocate(this%pother)
    if (allocated(this%sep_size)) deallocate(this%sep_size)
    if (allocated(this%perm_u2c)) deallocate(this%perm_u2c)
    if (allocated(this%share_weight)) deallocate(this%share_weight)
    call this%fac_a_ll%free()
    call this%a_sl%free()
    if (allocated(this%xp)) deallocate(this%xp)
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%vl)) deallocate(this%vl)
    if (allocated(this%vc)) deallocate(this%vc)
    if (allocated(this%vx)) deallocate(this%vx)
    if (allocated(this%combuf)) deallocate(this%combuf)

    this%pcoord = 0
    this%nsep = 0
    this%plevels = 0
    this%un = 0
    this%cn = 0
    this%ln = 0
    this%sn = 0
    this%xn = 0
    this%ln_fac = 0
    this%null_space = .false.
    this%tree_initialized = .false.
  end subroutine pnpn2_coarse_xxt_free

  subroutine pnpn2_coarse_xxt_setup(this, id, nz, ai, aj, av)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    integer(kind=i8), intent(in) :: id(:)
    integer, intent(in) :: nz
    integer, intent(in) :: ai(nz), aj(nz)
    real(kind=rp), intent(in) :: av(nz)
    type(xxt_dof_t), allocatable :: dof(:)
    type(xxt_csr_mat_t) :: a_ll, a_ss
    integer, allocatable :: perm_x2c(:)
    integer :: xcol, max_sep

    if (.not. this%tree_initialized) then
      call this%init_tree(pe_rank, pe_size, this%null_space)
    end if

    if (allocated(this%sep_size)) deallocate(this%sep_size)
    if (allocated(this%perm_u2c)) deallocate(this%perm_u2c)
    if (allocated(this%share_weight)) deallocate(this%share_weight)
    call this%fac_a_ll%free()
    call this%a_sl%free()
    if (allocated(this%xp)) deallocate(this%xp)
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%vl)) deallocate(this%vl)
    if (allocated(this%vc)) deallocate(this%vc)
    if (allocated(this%vx)) deallocate(this%vx)
    if (allocated(this%combuf)) deallocate(this%combuf)

    call xxt_discover_dofs(this, id, dof)
    call xxt_discover_sep_sizes(this, dof)

    if (this%null_space) call xxt_build_share_weights(this, dof)
    call xxt_discover_sep_ids(this, dof, perm_x2c)

    this%ln_fac = this%ln
    if (this%null_space .and. this%xn .eq. 0) then
      this%ln_fac = this%ln - 1
      call xxt_separate_matrix(nz, ai, aj, av, this%perm_u2c, this%ln_fac, 1, &
           a_ll, this%a_sl, a_ss)
    else
      call xxt_separate_matrix(nz, ai, aj, av, this%perm_u2c, this%ln, this%sn, &
           a_ll, this%a_sl, a_ss)
    end if

    if (this%ln_fac .gt. 0) then
      call this%fac_a_ll%factor(a_ll%row_ptr, a_ll%col_ind, a_ll%val)
    else
      call this%fac_a_ll%free()
    end if
    call a_ll%free()

    allocate(this%vl(max(0, this%ln)))
    allocate(this%vc(max(0, this%cn)))
    allocate(this%vx(max(0, this%xn)))
    max_sep = max(1, this%xn)
    allocate(this%combuf(max_sep))

    xcol = this%xn
    if (this%null_space .and. xcol .gt. 0) xcol = xcol - 1
    if (xcol .gt. 0) then
      call xxt_orthogonalize(this, a_ss, perm_x2c)
    else
      allocate(this%xp(1), this%x(0))
      this%xp(1) = 1
    end if

    call a_ss%free()
    if (allocated(perm_x2c)) deallocate(perm_x2c)
    if (allocated(dof)) deallocate(dof)
  end subroutine pnpn2_coarse_xxt_setup

  subroutine xxt_diag_stats(a, diag_count, diag_min)
    type(xxt_csr_mat_t), intent(in) :: a
    integer, intent(out) :: diag_count
    real(kind=rp), intent(out) :: diag_min
    integer :: i, p

    diag_count = 0
    diag_min = huge(1.0_rp)
    do i = 1, a%n
      do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
        if (a%col_ind(p) .eq. i) then
          diag_count = diag_count + 1
          diag_min = min(diag_min, a%val(p))
          exit
        end if
      end do
    end do
    if (diag_count .eq. 0) diag_min = 0.0_rp
  end subroutine xxt_diag_stats

  subroutine pnpn2_coarse_xxt_solve(this, x, b)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    real(kind=rp), intent(out) :: x(:)
    real(kind=rp), intent(in) :: b(:)
    integer :: i, p, xcol
    real(kind=rp) :: s

    if (size(x) .ne. this%un .or. size(b) .ne. this%un) then
      call neko_error('XXT solve received vectors of invalid size.')
    end if

    this%vc = 0.0_rp
    do i = 1, this%un
      p = this%perm_u2c(i)
      if (p .ge. 1) this%vc(p) = this%vc(p) + b(i)
    end do

    xcol = this%xn
    if (xcol .gt. 0 .and. (.not. this%null_space .or. xcol .gt. 1)) then
      if (this%null_space) xcol = xcol - 1
      if (this%ln .gt. 0) then
        this%vl = 0.0_rp
        call this%fac_a_ll%solve(this%vc(1:this%ln), this%vc(1:this%ln))
      end if
      call xxt_apply_m_asl(this%vc(this%ln + 1:this%cn), this%sn, this, &
           this%vc(1:this%ln))
      call xxt_apply_xt(this%vx(1:xcol), xcol, this, this%vc(this%ln + 1:this%cn))
      call xxt_apply_qqt(this, this%vx(1:xcol), xcol, 0)
      call xxt_apply_x(this%vc(this%ln + 1:this%cn), this%sn, this, this%vx(1:xcol), xcol)
      this%vl = 0.0_rp
      call xxt_apply_p_als(this%vl, this, this%vc(this%ln + 1:this%cn), this%sn)
      if (this%ln .gt. 0) call this%fac_a_ll%solve(this%vl, this%vl)
      this%vc(1:this%ln) = this%vc(1:this%ln) - this%vl
    else
      if (this%ln_fac .gt. 0) then
        call this%fac_a_ll%solve(this%vc(1:this%ln_fac), this%vc(1:this%ln_fac))
      end if
      if (this%null_space) then
        if (this%xn .eq. 0) then
          this%vc(this%ln) = 0.0_rp
        else if (this%sn .eq. 1) then
          this%vc(this%ln + 1) = 0.0_rp
        end if
      end if
    end if

    if (this%null_space) then
      s = dot_product(this%share_weight, this%vc)
      s = xxt_sum_scalar(this, s, this%xn, 0)
      this%vc = this%vc - s
    end if

    do i = 1, this%un
      p = this%perm_u2c(i)
      if (p .ge. 1) then
        x(i) = this%vc(p)
      else
        x(i) = 0.0_rp
      end if
    end do
  end subroutine pnpn2_coarse_xxt_solve

  subroutine xxt_discover_dofs(this, id, dof)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    integer(kind=i8), intent(in) :: id(:)
    type(xxt_dof_t), allocatable, intent(out) :: dof(:)
    integer, allocatable :: order(:), old_to_new(:), owner_local(:), owner_global(:)
    integer, allocatable :: count_local(:), count_global(:), dof_idx(:)
    integer, allocatable :: owner_all(:,:)
    integer(kind=i8), allocatable :: uniq_id(:)
    integer :: i, cn, max_gid, ierr, gid, p

    this%un = size(id)
    allocate(this%perm_u2c(this%un))
    this%perm_u2c = -1

    if (this%un .eq. 0) then
      this%cn = 0
      allocate(dof(0))
      return
    end if

    allocate(order(this%un))
    do i = 1, this%un
      order(i) = i
    end do
    call sort_idx_i8(order, id, 1, this%un)

    allocate(uniq_id(this%un), dof_idx(this%un))
    cn = 0
    do i = 1, this%un
      if (id(order(i)) .eq. 0_i8) then
        this%perm_u2c(order(i)) = -1
      else
        if (cn .eq. 0 .or. id(order(i)) .ne. uniq_id(cn)) then
          cn = cn + 1
          uniq_id(cn) = id(order(i))
        end if
        dof_idx(order(i)) = cn
      end if
    end do

    this%cn = cn
    allocate(dof(cn))
    do i = 1, cn
      dof(i)%id = uniq_id(i)
    end do

    max_gid = 0
    do i = 1, cn
      max_gid = max(max_gid, int(dof(i)%id))
    end do
    call MPI_Allreduce(max_gid, max_gid, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    allocate(owner_local(max_gid), owner_global(max_gid))
    allocate(owner_all(max_gid, pe_size))
    allocate(count_local(max_gid), count_global(max_gid))
    owner_local = 0
    owner_global = 0
    count_local = 0
    do i = 1, cn
      gid = int(dof(i)%id)
      owner_local(gid) = this%pcoord
      count_local(gid) = 1
    end do

    call MPI_Allgather(owner_local, max_gid, MPI_INTEGER, owner_all, max_gid, &
         MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allreduce(count_local, count_global, max_gid, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)

    do gid = 1, max_gid
      do p = 1, pe_size
        owner_global(gid) = xxt_bpr(owner_global(gid), owner_all(gid, p))
      end do
    end do

    do i = 1, cn
      gid = int(dof(i)%id)
      dof(i)%count = count_global(gid)
      dof(i)%level = this%nsep - 1 - ilog2_pos(owner_global(gid))
    end do

    deallocate(owner_local, owner_global, owner_all, count_local, count_global)

    deallocate(order)
    allocate(order(cn), old_to_new(cn))
    do i = 1, cn
      order(i) = i
    end do
    call sort_idx_dof(order, dof, 1, cn)
    do i = 1, cn
      old_to_new(order(i)) = i
    end do
    call permute_dofs(dof, order)
    do i = 1, this%un
      if (id(i) .eq. 0_i8) then
        this%perm_u2c(i) = -1
      else
        this%perm_u2c(i) = old_to_new(dof_idx(i))
      end if
    end do

    deallocate(order, old_to_new, uniq_id, dof_idx)
  end subroutine xxt_discover_dofs

  subroutine xxt_discover_sep_sizes(this, dof)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    type(xxt_dof_t), intent(in) :: dof(:)
    real(kind=rp), allocatable :: v(:), recv(:)
    integer :: i, lvl, s, other, ierr

    allocate(v(this%nsep), recv(this%nsep))
    v = 0.0_rp
    recv = 0.0_rp
    do i = 1, size(dof)
      v(dof(i)%level + 1) = v(dof(i)%level + 1) + &
           1.0_rp / real(dof(i)%count, rp)
    end do

    do lvl = 0, this%plevels - 1
      other = this%pother(lvl + 1)
      s = this%nsep - (lvl + 1)
      if (other .lt. 0) then
        call MPI_Send(v(lvl + 2), s, MPI_REAL_RP(), -other - 1, s, NEKO_COMM, ierr)
      else
        call MPI_Recv(recv(lvl + 2), s, MPI_REAL_RP(), other, s, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
        v(lvl + 2:this%nsep) = v(lvl + 2:this%nsep) + recv(lvl + 2:this%nsep)
      end if
    end do

    do lvl = this%plevels - 1, 0, -1
      other = this%pother(lvl + 1)
      s = this%nsep - (lvl + 1)
      if (other .lt. 0) then
        call MPI_Recv(v(lvl + 2), s, MPI_REAL_RP(), -other - 1, s, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
      else
        call MPI_Send(v(lvl + 2), s, MPI_REAL_RP(), other, s, NEKO_COMM, ierr)
      end if
    end do

    allocate(this%sep_size(this%nsep))
    this%xn = 0
    do i = 1, this%nsep
      this%sep_size(i) = int(v(i) + 0.1_rp)
      this%xn = this%xn + this%sep_size(i)
    end do
    this%ln = this%sep_size(1)
    this%sn = this%cn - this%ln
    this%xn = this%xn - this%ln

    deallocate(v, recv)
  end subroutine xxt_discover_sep_sizes

  subroutine xxt_build_share_weights(this, dof)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    type(xxt_dof_t), intent(in) :: dof(:)
    real(kind=rp) :: count
    integer :: i

    count = 0.0_rp
    do i = 1, this%cn
      count = count + 1.0_rp / real(dof(i)%count, rp)
    end do
    count = 1.0_rp / xxt_sum_scalar(this, count, this%xn, 0)
    allocate(this%share_weight(this%cn))
    do i = 1, this%cn
      this%share_weight(i) = count / real(dof(i)%count, rp)
    end do
  end subroutine xxt_build_share_weights

  subroutine xxt_discover_sep_ids(this, dof, perm_x2c)
    class(pnpn2_coarse_xxt_t), intent(in) :: this
    type(xxt_dof_t), intent(in) :: dof(:)
    integer, allocatable, intent(out) :: perm_x2c(:)
    integer(kind=i8), allocatable :: xid(:), recv(:), work(:)
    integer :: lvl, s, size_max, other, ierr, offset

    allocate(perm_x2c(max(0, this%xn)))
    if (this%xn .eq. 0) return

    size_max = 0
    do lvl = 2, this%nsep
      size_max = max(size_max, this%sep_size(lvl))
    end do
    allocate(xid(this%xn), recv(this%xn), work(2 * max(1, size_max)))
    call xxt_init_sep_ids(this, dof, xid)

    if (this%plevels .gt. 0) then
      s = this%xn
      do lvl = 0, this%plevels - 1
        other = this%pother(lvl + 1)
        if (lvl .eq. 0) then
          offset = 0
        else
          offset = sum(this%sep_size(2:lvl+1))
        end if
        if (other .lt. 0) then
          call MPI_Send(xid(offset + 1), s, MPI_INTEGER8, &
               -other - 1, s, NEKO_COMM, ierr)
        else
          call MPI_Recv(recv(1), s, MPI_INTEGER8, other, s, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
          call xxt_merge_sep_ids(this, xid(offset + 1:), &
               recv(1:s), work, lvl + 2)
        end if
        if (lvl .eq. this%plevels - 1) exit
        if (this%sep_size(lvl + 2) .ge. s) exit
        s = s - this%sep_size(lvl + 2)
      end do

      do lvl = lvl, 0, -1
        other = this%pother(lvl + 1)
        if (lvl .eq. 0) then
          offset = 0
        else
          offset = sum(this%sep_size(2:lvl+1))
        end if
        s = this%xn - offset
        if (other .lt. 0) then
          call MPI_Recv(xid(offset + 1), s, MPI_INTEGER8, &
               -other - 1, s, NEKO_COMM, MPI_STATUS_IGNORE, ierr)
        else
          call MPI_Send(xid(offset + 1), s, MPI_INTEGER8, &
               other, s, NEKO_COMM, ierr)
        end if
      end do
    end if

    call xxt_find_perm_x2c(this, dof, xid, perm_x2c)
    deallocate(xid, recv, work)
  end subroutine xxt_discover_sep_ids

  subroutine xxt_init_sep_ids(this, dof, xid)
    class(pnpn2_coarse_xxt_t), intent(in) :: this
    type(xxt_dof_t), intent(in) :: dof(:)
    integer(kind=i8), intent(out) :: xid(:)
    integer :: i, p, lvl, idx

    xid = 0_i8
    if (this%nsep .eq. 1) return
    p = 1
    lvl = 1
    idx = this%ln + 1
    do while (lvl .le. this%nsep - 1)
      do while (idx .le. this%cn .and. dof(idx)%level .eq. lvl)
        xid(p) = dof(idx)%id
        p = p + 1
        idx = idx + 1
      end do
      p = sum(this%sep_size(2:lvl+1)) + 1
      lvl = lvl + 1
    end do
  end subroutine xxt_init_sep_ids

  subroutine xxt_merge_sep_ids(this, sep_id, other, work, s0)
    class(pnpn2_coarse_xxt_t), intent(in) :: this
    integer(kind=i8), intent(inout) :: sep_id(:)
    integer(kind=i8), intent(in) :: other(:)
    integer(kind=i8), intent(inout) :: work(:)
    integer, intent(in) :: s0
    integer :: s, size, p, q, nwork

    p = 1
    q = 1
    do s = s0, this%nsep
      size = this%sep_size(s)
      if (size .le. 0) cycle
      work(1:size) = sep_id(p:p + size - 1)
      work(size + 1:2 * size) = other(q:q + size - 1)
      call sort_i8(work, 1, 2 * size)
      nwork = xxt_unique_nonzero(work, 2 * size)
      if (nwork .gt. size) then
        call neko_error('XXT separator-id merge produced inconsistent separator sizes.')
      end if
      sep_id(p:p + nwork - 1) = work(1:nwork)
      if (nwork .lt. size) sep_id(p + nwork:p + size - 1) = 0_i8
      p = p + size
      q = q + size
    end do
  end subroutine xxt_merge_sep_ids

  integer function xxt_unique_nonzero(v, n) result(m)
    integer(kind=i8), intent(inout) :: v(:)
    integer, intent(in) :: n
    integer :: i
    integer(kind=i8) :: last

    m = 0
    last = -1_i8
    do i = 1, n
      if (v(i) .eq. 0_i8) cycle
      if (m .eq. 0 .or. v(i) .ne. last) then
        m = m + 1
        v(m) = v(i)
        last = v(i)
      end if
    end do
  end function xxt_unique_nonzero

  subroutine xxt_find_perm_x2c(this, dof, xid, perm_x2c)
    class(pnpn2_coarse_xxt_t), intent(in) :: this
    type(xxt_dof_t), intent(in) :: dof(:)
    integer(kind=i8), intent(in) :: xid(:)
    integer, intent(out) :: perm_x2c(:)
    integer :: i, j, h

    i = this%ln + 1
    h = 1
    do j = 1, this%xn
      if (i .le. this%cn .and. xid(j) .eq. dof(i)%id) then
        perm_x2c(j) = i - this%ln
        i = i + 1
      else
        perm_x2c(j) = -1
      end if
    end do
  end subroutine xxt_find_perm_x2c

  subroutine xxt_allocate_x(this, perm_x2c)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    integer, intent(in) :: perm_x2c(:)
    integer :: xcol, i, h

    xcol = this%xn
    if (this%null_space .and. xcol .gt. 0) xcol = xcol - 1
    allocate(this%xp(xcol + 1))
    this%xp(1) = 1
    h = 0
    do i = 1, xcol
      if (perm_x2c(i) .ne. -1) h = h + 1
      this%xp(i + 1) = this%xp(i) + h
    end do
    allocate(this%x(max(0, this%xp(xcol + 1) - 1)))
  end subroutine xxt_allocate_x

  subroutine xxt_orthogonalize(this, a_ss, perm_x2c)
    class(pnpn2_coarse_xxt_t), intent(inout) :: this
    type(xxt_csr_mat_t), intent(in) :: a_ss
    integer, intent(in) :: perm_x2c(:)
    real(kind=rp), allocatable :: vs(:), svs(:), vx(:), vl(:)
    integer :: i, j, ns, ui, xcol, start_idx
    real(kind=rp) :: ytsy

    call xxt_allocate_x(this, perm_x2c)

    xcol = this%xn
    if (this%null_space .and. xcol .gt. 0) xcol = xcol - 1
    allocate(vs(max(1, this%sn)), svs(max(1, this%sn)), vx(max(1, xcol)), &
         vl(max(1, this%ln)))

    do i = 1, xcol
      ns = this%xp(i + 1) - this%xp(i)
      ui = perm_x2c(i)

      if (ui .eq. -1) then
        if (i .gt. 1) vx(1:i - 1) = 0.0_rp
      else
        call xxt_apply_s_col(vs, this, a_ss, ui, vl)
        if (i .gt. 1) call xxt_apply_xt(vx(1:i - 1), i - 1, this, vs)
      end if

      if (i .gt. 1) then
        call xxt_apply_qqt(this, vx(1:i - 1), i - 1, xcol - i)
        call xxt_apply_x(vs, ns, this, vx(1:i - 1), i - 1)
      else
        vs(1:ns) = 0.0_rp
      end if

      if (ui .ne. -1) vs(ui) = 1.0_rp
      call xxt_apply_s(svs, ns, this, a_ss, vs, vl)
      ytsy = dot_product(vs(1:ns), svs(1:ns))
      ytsy = xxt_sum_scalar(this, ytsy, i, xcol - i)
      if (ytsy .lt. epsilon(1.0_rp) / 128.0_rp) then
        ytsy = 0.0_rp
      else
        ytsy = 1.0_rp / sqrt(ytsy)
      end if
      start_idx = this%xp(i)
      do j = 1, ns
        this%x(start_idx + j - 1) = ytsy * vs(j)
      end do
    end do

    deallocate(vs, svs, vx, vl)
  end subroutine xxt_orthogonalize

  subroutine xxt_apply_p_als(vl, data, vs, ns)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    integer, intent(in) :: ns
    real(kind=rp), intent(inout) :: vl(:)
    real(kind=rp), intent(in) :: vs(:)
    integer :: i, p

    do i = 1, ns
      do p = data%a_sl%row_ptr(i), data%a_sl%row_ptr(i + 1) - 1
        vl(data%a_sl%col_ind(p)) = vl(data%a_sl%col_ind(p)) + &
             data%a_sl%val(p) * vs(i)
      end do
    end do
  end subroutine xxt_apply_p_als

  subroutine xxt_apply_m_asl(vs, ns, data, vl)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    integer, intent(in) :: ns
    real(kind=rp), intent(inout) :: vs(:)
    real(kind=rp), intent(in) :: vl(:)
    integer :: i, p

    do i = 1, ns
      do p = data%a_sl%row_ptr(i), data%a_sl%row_ptr(i + 1) - 1
        vs(i) = vs(i) - data%a_sl%val(p) * vl(data%a_sl%col_ind(p))
      end do
    end do
  end subroutine xxt_apply_m_asl

  subroutine xxt_apply_s_col(vs, data, a_ss, ei, vl)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    type(xxt_csr_mat_t), intent(in) :: a_ss
    integer, intent(in) :: ei
    real(kind=rp), intent(inout) :: vs(:)
    real(kind=rp), intent(inout) :: vl(:)
    integer :: p, pe, j

    if (ei .gt. 1) vs(1:ei - 1) = 0.0_rp
    do p = a_ss%row_ptr(ei), a_ss%row_ptr(ei + 1) - 1
      j = a_ss%col_ind(p)
      if (j .ge. ei) exit
      vs(j) = -a_ss%val(p)
    end do
    if (data%ln .gt. 0) vl = 0.0_rp
    do p = data%a_sl%row_ptr(ei), data%a_sl%row_ptr(ei + 1) - 1
      vl(data%a_sl%col_ind(p)) = -data%a_sl%val(p)
    end do
    if (data%ln .gt. 0) call data%fac_a_ll%solve(vl, vl)
    call xxt_apply_m_asl(vs, ei - 1, data, vl)
  end subroutine xxt_apply_s_col

  subroutine xxt_apply_s(svs, ns, data, a_ss, vs, vl)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    type(xxt_csr_mat_t), intent(in) :: a_ss
    integer, intent(in) :: ns
    real(kind=rp), intent(out) :: svs(:)
    real(kind=rp), intent(in) :: vs(:)
    real(kind=rp), intent(inout) :: vl(:)
    integer :: i, p

    do i = 1, ns
      svs(i) = 0.0_rp
      do p = a_ss%row_ptr(i), a_ss%row_ptr(i + 1) - 1
        if (a_ss%col_ind(p) .gt. ns) exit
        svs(i) = svs(i) + a_ss%val(p) * vs(a_ss%col_ind(p))
      end do
    end do
    if (data%ln .gt. 0) then
      vl = 0.0_rp
      call xxt_apply_p_als(vl, data, vs, ns)
      call data%fac_a_ll%solve(vl, vl)
      call xxt_apply_m_asl(svs, ns, data, vl)
    end if
  end subroutine xxt_apply_s

  subroutine xxt_apply_xt(vx, nx, data, vs)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    integer, intent(in) :: nx
    real(kind=rp), intent(out) :: vx(:)
    real(kind=rp), intent(in) :: vs(:)
    integer :: i

    do i = 1, nx
      vx(i) = dot_product(vs(1:data%xp(i + 1) - data%xp(i)), &
           data%x(data%xp(i):data%xp(i + 1) - 1))
    end do
  end subroutine xxt_apply_xt

  subroutine xxt_apply_x(vs, ns, data, vx, nx)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    integer, intent(in) :: ns, nx
    real(kind=rp), intent(out) :: vs(:)
    real(kind=rp), intent(in) :: vx(:)
    integer :: i, j, n

    vs(1:ns) = 0.0_rp
    do i = 1, nx
      n = data%xp(i + 1) - data%xp(i)
      do j = 1, n
        vs(j) = vs(j) + data%x(data%xp(i) + j - 1) * vx(i)
      end do
    end do
  end subroutine xxt_apply_x

  subroutine xxt_apply_qqt(data, v, n, tag)
    class(pnpn2_coarse_xxt_t), intent(inout) :: data
    integer, intent(in) :: n, tag
    real(kind=rp), intent(inout) :: v(:)
    integer :: lvl, offset, size, ss, other, ierr
    real(kind=rp), allocatable :: recv(:)

    if (n .eq. 0 .or. data%plevels .eq. 0) return

    allocate(recv(n))
    offset = 0
    size = n
    do lvl = 0, data%plevels - 1
      other = data%pother(lvl + 1)
      if (other .lt. 0) then
        call MPI_Send(v(offset + 1), size, MPI_REAL_RP(), -other - 1, &
             2 * tag, NEKO_COMM, ierr)
      else
        call MPI_Recv(recv(1), size, MPI_REAL_RP(), other, 2 * tag, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
        v(offset + 1:offset + size) = v(offset + 1:offset + size) + recv(1:size)
      end if
      ss = data%sep_size(lvl + 2)
      if (ss .ge. size .or. lvl .eq. data%plevels - 1) exit
      offset = offset + ss
      size = size - ss
    end do

    do lvl = lvl, 0, -1
      other = data%pother(lvl + 1)
      if (other .lt. 0) then
        call MPI_Recv(v(offset + 1), size, MPI_REAL_RP(), -other - 1, &
             2 * tag, NEKO_COMM, MPI_STATUS_IGNORE, ierr)
      else
        call MPI_Send(v(offset + 1), size, MPI_REAL_RP(), other, 2 * tag, &
             NEKO_COMM, ierr)
      end if
      if (lvl .eq. 0) exit
      ss = data%sep_size(lvl + 1)
      offset = offset - ss
      size = size + ss
    end do

    deallocate(recv)
  end subroutine xxt_apply_qqt

  real(kind=rp) function xxt_sum_scalar(data, v, n, tag) result(total)
    class(pnpn2_coarse_xxt_t), intent(in) :: data
    real(kind=rp), intent(in) :: v
    integer, intent(in) :: n, tag
    real(kind=rp) :: r
    integer :: lvl, size, ss, other, ierr

    total = v
    if (n .eq. 0 .or. data%plevels .eq. 0) return

    size = n
    do lvl = 0, data%plevels - 1
      other = data%pother(lvl + 1)
      if (other .lt. 0) then
        call MPI_Send(total, 1, MPI_REAL_RP(), -other - 1, 2 * tag + 1, NEKO_COMM, ierr)
      else
        call MPI_Recv(r, 1, MPI_REAL_RP(), other, 2 * tag + 1, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
        total = total + r
      end if
      ss = data%sep_size(lvl + 2)
      if (ss .ge. size .or. lvl .eq. data%plevels - 1) exit
      size = size - ss
    end do

    do lvl = lvl, 0, -1
      other = data%pother(lvl + 1)
      if (other .lt. 0) then
        call MPI_Recv(total, 1, MPI_REAL_RP(), -other - 1, 2 * tag + 1, NEKO_COMM, &
             MPI_STATUS_IGNORE, ierr)
      else
        call MPI_Send(total, 1, MPI_REAL_RP(), other, 2 * tag + 1, NEKO_COMM, ierr)
      end if
      if (lvl .eq. 0) exit
      ss = data%sep_size(lvl + 1)
      size = size + ss
    end do
  end function xxt_sum_scalar

  subroutine xxt_separate_matrix(nz, ai, aj, av, perm, ln, sn, a_ll, a_sl, a_ss)
    integer, intent(in) :: nz, ai(nz), aj(nz), perm(:), ln, sn
    real(kind=rp), intent(in) :: av(nz)
    type(xxt_csr_mat_t), intent(inout) :: a_ll, a_sl, a_ss
    type(xxt_yale_mat_t), allocatable :: mat_ll(:), mat_sl(:), mat_ss(:)
    integer :: k, i, j, n_ll, n_sl, n_ss

    call a_ll%free()
    call a_sl%free()
    call a_ss%free()

    allocate(mat_ll(2 * max(1, nz)))
    allocate(mat_sl(2 * max(1, nz)))
    allocate(mat_ss(2 * max(1, nz)))
    n_ll = 0
    n_sl = 0
    n_ss = 0

    do k = 1, nz
      i = perm(ai(k))
      j = perm(aj(k))
      if (i .lt. 1 .or. j .lt. 1) cycle
      if (abs(av(k)) .le. tiny(1.0_rp)) cycle

      if (i .le. ln) then
        if (j .le. ln) then
          n_ll = n_ll + 1
          mat_ll(n_ll)%i = i
          mat_ll(n_ll)%j = j
          mat_ll(n_ll)%v = av(k)
        end if
      else
        if (j .le. ln) then
          n_sl = n_sl + 1
          mat_sl(n_sl)%i = i - ln
          mat_sl(n_sl)%j = j
          mat_sl(n_sl)%v = av(k)
        else
          n_ss = n_ss + 1
          mat_ss(n_ss)%i = i - ln
          mat_ss(n_ss)%j = j - ln
          mat_ss(n_ss)%v = av(k)
        end if
      end if
    end do

    call xxt_condense_matrix(mat_ll, n_ll, ln, a_ll)
    call xxt_condense_matrix(mat_sl, n_sl, sn, a_sl)
    call xxt_condense_matrix(mat_ss, n_ss, sn, a_ss)

    deallocate(mat_ll, mat_sl, mat_ss)
  end subroutine xxt_separate_matrix

  subroutine xxt_condense_matrix(mat, n_entries, nr, out)
    type(xxt_yale_mat_t), intent(inout) :: mat(:)
    integer, intent(inout) :: n_entries
    integer, intent(in) :: nr
    type(xxt_csr_mat_t), intent(inout) :: out
    integer :: k, nz, row, pos

    call out%free()

    if (nr .eq. 0) then
      out%n = 0
      allocate(out%row_ptr(1), out%col_ind(0), out%val(0))
      out%row_ptr(1) = 1
      return
    end if

    if (n_entries .gt. 1) call sort_yale(mat, 1, n_entries)

    if (n_entries .gt. 1) then
      nz = 1
      do k = 2, n_entries
        if (mat(k)%i .eq. mat(nz)%i .and. mat(k)%j .eq. mat(nz)%j) then
          mat(nz)%v = mat(nz)%v + mat(k)%v
        else
          nz = nz + 1
          if (nz .ne. k) mat(nz) = mat(k)
        end if
      end do
      n_entries = nz
    end if

    out%n = nr
    allocate(out%row_ptr(nr + 1), out%col_ind(max(1, n_entries)), &
         out%val(max(1, n_entries)))
    out%row_ptr = 0

    if (n_entries .gt. 0) then
      do k = 1, n_entries
        out%row_ptr(mat(k)%i + 1) = out%row_ptr(mat(k)%i + 1) + 1
        out%col_ind(k) = mat(k)%j
        out%val(k) = mat(k)%v
      end do
    end if

    out%row_ptr(1) = 1
    do row = 1, nr
      out%row_ptr(row + 1) = out%row_ptr(row + 1) + out%row_ptr(row)
    end do

    if (n_entries .eq. 0) then
      deallocate(out%col_ind, out%val)
      allocate(out%col_ind(0), out%val(0))
    else
      do row = 1, nr
        do pos = out%row_ptr(row), out%row_ptr(row + 1) - 2
          if (out%col_ind(pos) .gt. out%col_ind(pos + 1)) then
            call neko_error('XXT condense_matrix produced unsorted CSR columns.')
          end if
        end do
      end do
    end if
  end subroutine xxt_condense_matrix

  recursive subroutine sort_int(v, lo, hi)
    integer, intent(inout) :: v(:)
    integer, intent(in) :: lo, hi
    integer :: i, j, pivot, tmp

    if (lo .ge. hi) return
    i = lo
    j = hi
    pivot = v((lo + hi) / 2)
    do
      do while (v(i) .lt. pivot)
        i = i + 1
      end do
      do while (v(j) .gt. pivot)
        j = j - 1
      end do
      if (i .le. j) then
        tmp = v(i)
        v(i) = v(j)
        v(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do
    if (lo .lt. j) call sort_int(v, lo, j)
    if (i .lt. hi) call sort_int(v, i, hi)
  end subroutine sort_int

  recursive subroutine sort_i8(v, lo, hi)
    integer(kind=i8), intent(inout) :: v(:)
    integer, intent(in) :: lo, hi
    integer :: i, j
    integer(kind=i8) :: pivot, tmp

    if (lo .ge. hi) return
    i = lo
    j = hi
    pivot = v((lo + hi) / 2)
    do
      do while (v(i) .lt. pivot)
        i = i + 1
      end do
      do while (v(j) .gt. pivot)
        j = j - 1
      end do
      if (i .le. j) then
        tmp = v(i)
        v(i) = v(j)
        v(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do
    if (lo .lt. j) call sort_i8(v, lo, j)
    if (i .lt. hi) call sort_i8(v, i, hi)
  end subroutine sort_i8

  recursive subroutine sort_idx_i8(idx, key, lo, hi)
    integer, intent(inout) :: idx(:)
    integer(kind=i8), intent(in) :: key(:)
    integer, intent(in) :: lo, hi
    integer :: i, j, tmp
    integer(kind=i8) :: pivot

    if (lo .ge. hi) return
    i = lo
    j = hi
    pivot = key(idx((lo + hi) / 2))
    do
      do while (key(idx(i)) .lt. pivot)
        i = i + 1
      end do
      do while (key(idx(j)) .gt. pivot)
        j = j - 1
      end do
      if (i .le. j) then
        tmp = idx(i)
        idx(i) = idx(j)
        idx(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do
    if (lo .lt. j) call sort_idx_i8(idx, key, lo, j)
    if (i .lt. hi) call sort_idx_i8(idx, key, i, hi)
  end subroutine sort_idx_i8

  recursive subroutine sort_idx_dof(idx, dof, lo, hi)
    integer, intent(inout) :: idx(:)
    type(xxt_dof_t), intent(in) :: dof(:)
    integer, intent(in) :: lo, hi
    integer :: i, j, tmp
    type(xxt_dof_t) :: pivot

    if (lo .ge. hi) return
    i = lo
    j = hi
    pivot = dof(idx((lo + hi) / 2))
    do
      do while (dof_less(dof(idx(i)), pivot))
        i = i + 1
      end do
      do while (dof_less(pivot, dof(idx(j))))
        j = j - 1
      end do
      if (i .le. j) then
        tmp = idx(i)
        idx(i) = idx(j)
        idx(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do
    if (lo .lt. j) call sort_idx_dof(idx, dof, lo, j)
    if (i .lt. hi) call sort_idx_dof(idx, dof, i, hi)
  end subroutine sort_idx_dof

  subroutine permute_dofs(dof, order)
    type(xxt_dof_t), intent(inout) :: dof(:)
    integer, intent(in) :: order(:)
    type(xxt_dof_t), allocatable :: tmp(:)
    integer :: i

    allocate(tmp(size(dof)))
    do i = 1, size(dof)
      tmp(i) = dof(order(i))
    end do
    dof = tmp
    deallocate(tmp)
  end subroutine permute_dofs

  logical function dof_less(a, b) result(is_less)
    type(xxt_dof_t), intent(in) :: a, b
    is_less = (a%level .lt. b%level) .or. &
         (a%level .eq. b%level .and. a%count .lt. b%count) .or. &
         (a%level .eq. b%level .and. a%count .eq. b%count .and. a%id .lt. b%id)
  end function dof_less

  integer function ilog2_pos(v) result(lg)
    integer, intent(in) :: v
    integer :: t
    if (v .le. 0) call neko_error('XXT encountered a non-positive processor code.')
    lg = -1
    t = v
    do while (t .gt. 0)
      t = ishft(t, -1)
      lg = lg + 1
    end do
  end function ilog2_pos

  integer function xxt_bpr(a, b) result(r)
    integer, intent(in) :: a, b
    integer :: aa, bb

    r = a
    if (b .eq. 0) return
    if (r .eq. 0) then
      r = b
      return
    end if

    aa = r
    bb = b
    do
      if (aa .eq. bb) exit
      if (aa .lt. bb) then
        bb = ishft(bb, -1)
      else
        aa = ishft(aa, -1)
      end if
    end do
    r = aa
  end function xxt_bpr

  recursive subroutine sort_yale(v, lo, hi)
    type(xxt_yale_mat_t), intent(inout) :: v(:)
    integer, intent(in) :: lo, hi
    integer :: i, j
    type(xxt_yale_mat_t) :: pivot, tmp

    if (lo .ge. hi) return
    i = lo
    j = hi
    pivot = v((lo + hi) / 2)
    do
      do while (yale_less(v(i), pivot))
        i = i + 1
      end do
      do while (yale_less(pivot, v(j)))
        j = j - 1
      end do
      if (i .le. j) then
        tmp = v(i)
        v(i) = v(j)
        v(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do

    if (lo .lt. j) call sort_yale(v, lo, j)
    if (i .lt. hi) call sort_yale(v, i, hi)
  end subroutine sort_yale

  logical function yale_less(a, b) result(is_less)
    type(xxt_yale_mat_t), intent(in) :: a, b
    is_less = (a%i .lt. b%i) .or. (a%i .eq. b%i .and. a%j .lt. b%j)
  end function yale_less

  function MPI_REAL_RP() result(dtype)
    type(MPI_Datatype) :: dtype
    dtype = MPI_REAL_PRECISION
  end function MPI_REAL_RP

end module pnpn2_coarse_xxt
