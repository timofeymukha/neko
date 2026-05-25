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
!!
!! This module owns the lowest level of the milestone Pn/Pn-2 HSMG pressure
!! preconditioner.
!!
!! Conceptually, the module performs three tasks:
!! 1. construct the coarse-space basis functions on the fine H1 grid,
!! 2. assemble the condensed coarse Galerkin operator explicitly, and
!! 3. solve the condensed coarse system with a local dense LU factorisation.
!!
!! The implementation follows the same mathematical idea as Nek5000's coarse
!! H1 setup:
!! - the coarse basis is trilinear (or bilinear in 2D),
!! - the coarse operator is assembled by local Galerkin projection,
!! - shared coarse degrees of freedom are condensed to a unique global index,
!! - periodic/null-space cases are handled by removing the weighted coarse mean
!!   after solve.
!!
!! Unlike Nek5000's XXT-based coarse solve, this milestone implementation keeps
!! an explicit dense matrix on every rank.  That makes the data flow much easier
!! to inspect and document, which is useful while the Pn/Pn-2 HSMG path is
!! still being validated.
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

  !> Explicit condensed coarse solver used by the Pn/Pn-2 HSMG path.
  !!
  !! The type stores both the data needed to assemble the operator and the data
  !! needed to apply the factorised coarse solve repeatedly during GMRES:
  !! - mappings from local coarse dofs to condensed unique indices,
  !! - the dense LU factors and pivot array,
  !! - the trilinear coarse basis on the fine grid, and
  !! - work vectors for right-hand sides and solutions.
  type, public :: pnpn2_coarse_direct_t
    !> Number of local coarse-grid entries, including duplicates on shared dofs.
    integer :: n_local = 0
    !> Number of condensed unique coarse-grid unknowns.
    integer :: n_global = 0
    !> Largest raw coarse-grid global id seen before compaction.
    integer :: max_gid = 0
    !> Fine-grid point count in one horizontal direction.
    integer :: lx_fine = 0
    !> Fine-grid point count in the third direction.
    integer :: lz_fine = 0
    !> Number of local coarse nodes per element, i.e. `lx*ly*lz` on `Xh_crs`.
    integer :: nloc_crs = 0
    !> Map from raw coarse-grid ids to compact dense-system indices.
    integer, allocatable :: gid_to_idx(:)
    !> Map from each local coarse entry to its compact dense-system index.
    integer, allocatable :: local_to_idx(:)
    !> Pivot array produced by the dense LU factorisation.
    integer, allocatable :: ipiv(:)
    !> LU factors of the condensed coarse matrix.
    real(kind=rp), allocatable :: lu(:,:)
    !> Condensed coarse right-hand side work vector.
    real(kind=rp), allocatable :: rhs(:)
    !> Condensed coarse solution work vector.
    real(kind=rp), allocatable :: sol(:)
    !> Normalised reciprocal sharing counts used for mean removal.
    real(kind=rp), allocatable :: share_weight(:)
    !> Trilinear coarse basis functions sampled on the fine H1 grid.
    real(kind=rp), allocatable :: basis(:,:,:,:)
    !> True if the condensed matrix has a constant null-space.
    logical :: null_space = .false.
    !> True once the factorisation has been built successfully.
    logical :: initialized = .false.
    !> Fine-grid trial field used during Galerkin assembly.
    type(field_t) :: u
    !> Fine-grid operator image used during Galerkin assembly.
    type(field_t) :: w
    !> Fine-grid Helmholtz/Laplacian operator used during assembly.
    class(ax_t), allocatable :: ax
   contains
     !> Assemble and factorise the condensed coarse operator.
     procedure, pass(this) :: init => pnpn2_coarse_direct_init
     !> Release the factored coarse operator and all work storage.
     procedure, pass(this) :: free => pnpn2_coarse_direct_free
     !> Apply the already-factorised condensed coarse solve.
     procedure, pass(this) :: solve => pnpn2_coarse_direct_solve
  end type pnpn2_coarse_direct_t

contains

  !> Assemble and factorise the condensed coarse operator.
  !!
  !! @param coef_fine Fine-grid coefficient object used to apply the H1 operator
  !! during Galerkin assembly.
  !! @param dm_fine Fine-grid dofmap associated with `coef_fine`.
  !! @param dm_crs Coarse-grid dofmap on the trilinear/trilinear coarse space.
  !! @param null_space True if the condensed coarse operator has a constant
  !! null-space and therefore needs weighted mean removal after solve.
  !!
  !! The routine performs the full coarse setup:
  !! 1. compact raw coarse-grid ids into dense unique indices;
  !! 2. compute coarse sharing counts for null-space mean removal;
  !! 3. generate trilinear basis functions on the fine grid;
  !! 4. assemble local Galerkin element matrices;
  !! 5. condense them into one dense coarse matrix; and
  !! 6. factorise that matrix once with dense LU.
  subroutine pnpn2_coarse_direct_init(this, coef_fine, dm_fine, dm_crs, null_space)
    class(pnpn2_coarse_direct_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef_fine
    type(dofmap_t), intent(inout), target :: dm_fine
    type(dofmap_t), intent(inout), target :: dm_crs
    integer :: ierr, i, e, gid, idx, iloc, jloc
    integer :: local_max_gid
    logical, intent(in) :: null_space
    integer, allocatable :: present_local(:), present_global(:)
    integer, allocatable :: share_count_local(:), share_count_global(:)
    real(kind=rp), allocatable :: amat(:,:), amat_global(:,:), local_amat(:,:,:)
    real(kind=rp), allocatable :: h1_save(:), h2_save(:)
    logical :: ifh2_save
    integer :: ncoef
    real(kind=rp) :: weight_sum

    if (this%initialized) call this%free()

    ! Stage 1:
    ! Record the level sizes and determine the raw coarse-grid id range.  The
    ! coarse dofmap stores raw global ids, but the dense coarse solve needs a
    ! compact `1:n_global` numbering.
    this%n_local = dm_crs%size()
    this%lx_fine = coef_fine%Xh%lx
    this%lz_fine = coef_fine%Xh%lz
    this%nloc_crs = dm_crs%Xh%lx * dm_crs%Xh%ly * dm_crs%Xh%lz
    this%null_space = null_space
    local_max_gid = 0
    do i = 1, this%n_local
      local_max_gid = max(local_max_gid, int(dm_crs%dof(i,1,1,1)))
    end do
    call MPI_Allreduce(local_max_gid, this%max_gid, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    if (this%max_gid .le. 0) then
      call neko_error('PnPn-2 direct coarse solver found no coarse dofs.')
    end if

    ! Stage 2:
    ! Discover which raw coarse ids are present anywhere in the run and count
    ! how many local coarse entries contribute to each of them.  The latter is
    ! used only in null-space cases, where we remove the weighted coarse mean
    ! after solve.
    allocate(present_local(this%max_gid))
    allocate(present_global(this%max_gid))
    allocate(this%gid_to_idx(this%max_gid))
    allocate(this%local_to_idx(this%n_local))
    allocate(share_count_local(this%max_gid))
    allocate(share_count_global(this%max_gid))
    present_local = 0
    share_count_local = 0
    this%gid_to_idx = 0
    do i = 1, this%n_local
      gid = int(dm_crs%dof(i,1,1,1))
      if (gid .gt. 0) then
        present_local(gid) = 1
        share_count_local(gid) = share_count_local(gid) + 1
      end if
    end do
    call MPI_Allreduce(present_local, present_global, this%max_gid, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(share_count_local, share_count_global, this%max_gid, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

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

    ! Stage 3:
    ! Allocate the matrix/factorisation storage and create the fine-grid work
    ! fields used to apply the operator to coarse basis functions.
    allocate(amat(this%n_global, this%n_global))
    allocate(amat_global(this%n_global, this%n_global))
    allocate(local_amat(this%nloc_crs, this%nloc_crs, coef_fine%msh%nelv))
    allocate(this%lu(this%n_global, this%n_global))
    allocate(this%rhs(this%n_global))
    allocate(this%sol(this%n_global))
    allocate(this%share_weight(this%n_global))
    allocate(this%ipiv(this%n_global))
    allocate(this%basis(this%lx_fine, this%lx_fine, this%lz_fine, this%nloc_crs))
    call this%u%init(dm_fine, 'pnpn2_coarse_direct_u')
    call this%w%init(dm_fine, 'pnpn2_coarse_direct_w')
    call ax_helm_factory(this%ax, full_formulation = .false.)
    call generate_crs_basis(this%basis, coef_fine%Xh)

    ! The coarse operator is an H1 operator.  Save the current fine-grid
    ! coefficients, force `h1 = 1` and `h2 = 0`, then restore the original
    ! values after assembly.
    this%share_weight = 0.0_rp
    if (this%null_space) then
      do gid = 1, this%max_gid
        idx = this%gid_to_idx(gid)
        if (idx .gt. 0) then
          this%share_weight(idx) = 1.0_rp / real(share_count_global(gid), rp)
        end if
      end do
      weight_sum = sum(this%share_weight)
      if (weight_sum .le. tiny(1.0_rp)) then
        call neko_error('PnPn-2 direct coarse solver found zero null-space weight.')
      end if
      this%share_weight = this%share_weight / weight_sum
    end if
    ncoef = coef_fine%dof%size()
    allocate(h1_save(ncoef), h2_save(ncoef))
    h1_save = coef_fine%h1(1:ncoef,1,1,1)
    h2_save = coef_fine%h2(1:ncoef,1,1,1)
    ifh2_save = coef_fine%ifh2
    coef_fine%h1(1:ncoef,1,1,1) = 1.0_rp
    coef_fine%h2(1:ncoef,1,1,1) = 0.0_rp
    coef_fine%ifh2 = .false.

    ! Stage 4:
    ! Assemble local Galerkin matrices.
    !
    ! For each coarse basis function `b_j`, apply the fine-grid H1 operator and
    ! then form the local inner products `b_i^T A b_j`.
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

    ! Stage 5:
    ! Condense the element-local matrices into one dense global coarse matrix
    ! using the compact unique-dof numbering.
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

    ! Stage 6:
    ! Sum the dense matrix contributions across MPI ranks and prepare the final
    ! factorisation matrix.  In null-space cases we keep one pinned
    ! representative during the dense solve and then remove the weighted mean
    ! from the resulting solution.
    call MPI_Allreduce(amat, amat_global, this%n_global * this%n_global, &
         MPI_REAL_RP(), MPI_SUM, NEKO_COMM, ierr)

    if (this%null_space) then
      this%lu = amat_global
      this%lu(1,:) = 0.0_rp
      this%lu(:,1) = 0.0_rp
      this%lu(1,1) = 1.0_rp
    else
      this%lu = amat_global
    end if

    call dense_lu_factor(this%lu, this%ipiv, this%n_global)
    this%initialized = .true.

    deallocate(amat)
    deallocate(amat_global)
    deallocate(local_amat)
    deallocate(h1_save)
    deallocate(h2_save)
    deallocate(present_local)
    deallocate(present_global)
    deallocate(share_count_local)
    deallocate(share_count_global)
  end subroutine pnpn2_coarse_direct_init

  !> Release the factorised coarse operator and all work storage.
  subroutine pnpn2_coarse_direct_free(this)
    class(pnpn2_coarse_direct_t), intent(inout) :: this

    if (.not. this%initialized) return

    if (allocated(this%gid_to_idx)) deallocate(this%gid_to_idx)
    if (allocated(this%local_to_idx)) deallocate(this%local_to_idx)
    if (allocated(this%ipiv)) deallocate(this%ipiv)
    if (allocated(this%lu)) deallocate(this%lu)
    if (allocated(this%rhs)) deallocate(this%rhs)
    if (allocated(this%sol)) deallocate(this%sol)
    if (allocated(this%share_weight)) deallocate(this%share_weight)
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
    this%null_space = .false.
    this%initialized = .false.
  end subroutine pnpn2_coarse_direct_free

  !> Apply the already-factorised condensed coarse solve.
  !!
  !! @param e Output coarse correction on the local coarse grid.
  !! @param r Input coarse residual on the local coarse grid.
  !! @param coef Coefficient object of the coarse level.  It is not currently
  !! needed by the solve itself, but keeping the argument preserves the usual
  !! coarse-solver call shape in the surrounding code.
  !!
  !! The solve itself is straightforward:
  !! 1. sum duplicated local coarse entries into the condensed right-hand side;
  !! 2. allreduce that condensed vector so each rank sees the same global system;
  !! 3. apply the dense LU solve; and
  !! 4. expand the condensed solution back to each local coarse entry.
  !!
  !! In null-space cases a weighted coarse mean is removed after solve so that
  !! the returned correction lives in a representative zero-mean subspace.
  subroutine pnpn2_coarse_direct_solve(this, e, r, coef)
   class(pnpn2_coarse_direct_t), intent(inout) :: this
   real(kind=rp), intent(inout) :: e(this%n_local)
   real(kind=rp), intent(in) :: r(this%n_local)
   type(coef_t), intent(in) :: coef
   integer :: ierr, i, idx
   real(kind=rp) :: mean_sol

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

    if (this%null_space) this%sol(1) = 0.0_rp
    call dense_lu_solve(this%lu, this%ipiv, this%sol, this%n_global)

    if (this%null_space) then
      mean_sol = dot_product(this%share_weight, this%sol)
      this%sol = this%sol - mean_sol
    end if

    do i = 1, this%n_local
      e(i) = this%sol(this%local_to_idx(i))
    end do
  end subroutine pnpn2_coarse_direct_solve

  !> Factor a dense matrix in-place with partial pivoting.
  !!
  !! @param a Dense square matrix overwritten by its LU factors.
  !! @param ipiv Pivot array recording row swaps.
  !! @param n Matrix dimension.
  !!
  !! This is a small, self-contained LU factorisation specialised for the coarse
  !! matrices used here.  It keeps the implementation transparent and avoids
  !! introducing another dependency into the milestone branch.
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

  !> Solve a dense LU-factored linear system in-place.
  !!
  !! @param a Dense matrix containing LU factors.
  !! @param ipiv Pivot array returned by [dense_lu_factor](@ref
  !! pnpn2_coarse_direct::dense_lu_factor).
  !! @param b Right-hand side on entry and solution on exit.
  !! @param n System dimension.
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

  !> Return the MPI datatype matching Neko's runtime real precision.
  function MPI_REAL_RP() result(dtype)
    use comm, only : MPI_REAL_PRECISION
    use mpi_f08, only : MPI_Datatype
    type(MPI_Datatype) :: dtype

    dtype = MPI_REAL_PRECISION
  end function MPI_REAL_RP

  !> Convert a coarse local index and element number to the linear dof index.
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

  !> Recover the first tensor-product index from a linear coarse local index.
  integer function local_crs_i(Xh_crs, iloc) result(i)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    i = mod(iloc - 1, Xh_crs%lx) + 1
  end function local_crs_i

  !> Recover the second tensor-product index from a linear coarse local index.
  integer function local_crs_j(Xh_crs, iloc) result(j)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    j = mod((iloc - 1) / Xh_crs%lx, Xh_crs%ly) + 1
  end function local_crs_j

  !> Recover the third tensor-product index from a linear coarse local index.
  integer function local_crs_k(Xh_crs, iloc) result(k)
    type(space_t), intent(in) :: Xh_crs
    integer, intent(in) :: iloc
    k = ((iloc - 1) / (Xh_crs%lx * Xh_crs%ly)) + 1
  end function local_crs_k

  !> Generate the trilinear coarse basis sampled on the fine H1 grid.
  !!
  !! @param basis Basis array written in-place.
  !! @param Xh_space Fine H1 space on which the basis is sampled.
  !!
  !! The basis numbering follows Nek's corner ordering:
  !! - basis functions 1 and 2 differ in the first coordinate,
  !! - 3/4 and 7/8 additionally flip the second coordinate,
  !! - 5--8 additionally flip the third coordinate.
  !!
  !! Each basis function is therefore simply the tensor product of one-dimensional
  !! endpoint polynomials `z0 = (1-z)/2` and `z1 = (1+z)/2`.
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
