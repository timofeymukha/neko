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
!> Dedicated coarse H1 solver for the Pn/Pn-2 pressure HSMG path.
!!
!! Historically this module hosted the explicit dense coarse solve that Neko used
!! while the Pn/Pn-2 HSMG path was being brought up.  To reduce the remaining
!! distance to Nek5000's `crs_solve(xxth(...), e, r)` path, the implementation
!! now keeps the same public interface but builds and applies a sparse
!! coarse-grid operator instead:
!!
!! 1. build the trilinear/bilinear coarse basis on the fine H1 grid;
!! 2. assemble rank-local Galerkin triplets on the local coarse user dofs; and
!! 3. hand those local dofs and triplets to the Nek-style XXT backend.
!!
!! This keeps the public Neko interface unchanged while making the coarse solve
!! follow the same local-user-space setup structure as Nek5000 `crs_xxt.c`.
module pnpn2_coarse_direct
  use ax_product, only : ax_t, ax_helm_factory
  use coefs, only : coef_t
  use comm, only : pe_rank, pe_size
  use dofmap, only : dofmap_t
  use field, only : field_t
  use mpi_f08, only : MPI_Datatype
  use num_types, only : i8, rp
  use pnpn2_coarse_xxt, only : pnpn2_coarse_xxt_t
  use space, only : space_t
  use utils, only : neko_error
  implicit none
  private

  !> Sparse coarse solver used by the Pn/Pn-2 HSMG path.
  !!
  !! The public type name is kept for compatibility with the surrounding
  !! pressure-HSMG code, but the owned data is now CRS/CSR-oriented rather than
  !! dense-LU-oriented.
  type, public :: pnpn2_coarse_direct_t
    !> Number of local coarse-grid entries, including duplicated shared dofs.
    integer :: n_local = 0
    !> Number of unique coarse-grid unknowns after global compaction.
    integer :: n_global = 0
    !> Largest raw coarse-grid global id before compaction.
    integer :: max_gid = 0
    !> Fine-grid point count in one horizontal direction.
    integer :: lx_fine = 0
    !> Fine-grid point count in the third direction.
    integer :: lz_fine = 0
    !> Number of local coarse nodes per element, i.e. `lx*ly*lz` on `Xh_crs`.
    integer :: nloc_crs = 0
    !> Number of nonzeros in the replicated CSR matrix.
    integer :: nnz = 0
    !> Legacy CG iteration cap kept only while the old helper routines remain in-tree.
    integer :: cg_max_iter = 0
    !> Number of condensed coarse unknowns owned by separator 0.
    integer :: n_local_unique = 0
    !> Map from raw coarse-grid ids to compact unique indices.
    integer, allocatable :: gid_to_idx(:)
    !> Map from each local coarse entry to its compact unique index.
    integer, allocatable :: local_to_idx(:)
    !> CSR row pointers.
    integer, allocatable :: row_ptr(:)
    !> CSR column indices.
    integer, allocatable :: col_ind(:)
    !> Legacy Jacobi diagonal kept only while the old helper routines remain in-tree.
    real(kind=rp), allocatable :: diag_inv(:)
    !> CSR nonzero values.
    real(kind=rp), allocatable :: vals(:)
    !> Normalized weights used to remove the constant null mode.
    real(kind=rp), allocatable :: share_weight(:)
    !> Condensed coarse right-hand side.
    real(kind=rp), allocatable :: rhs(:)
    !> Condensed coarse solution.
    real(kind=rp), allocatable :: sol(:)
    !> Legacy work vectors kept only while the old helper routines remain in-tree.
    real(kind=rp), allocatable :: res(:)
    real(kind=rp), allocatable :: p(:)
    real(kind=rp), allocatable :: ap(:)
    real(kind=rp), allocatable :: z(:)
    !> Trilinear coarse basis sampled on the fine H1 grid.
    real(kind=rp), allocatable :: basis(:,:,:,:)
    !> Directly ported Nek5000 XXT state under construction.
    type(pnpn2_coarse_xxt_t) :: xxt
    !> True if the coarse operator has a constant null-space.
    logical :: null_space = .false.
    !> True once the sparse operator has been constructed successfully.
    logical :: initialized = .false.
    !> Fine-grid trial field used during Galerkin assembly.
    type(field_t) :: u
    !> Fine-grid operator image used during Galerkin assembly.
    type(field_t) :: w
    !> Fine-grid H1 operator used during Galerkin assembly.
    class(ax_t), allocatable :: ax
   contains
    procedure, pass(this) :: init => pnpn2_coarse_direct_init
    procedure, pass(this) :: free => pnpn2_coarse_direct_free
    procedure, pass(this) :: solve => pnpn2_coarse_direct_solve
  end type pnpn2_coarse_direct_t

contains

  !> Assemble the sparse coarse operator and the XXT-style condensed solve.
  subroutine pnpn2_coarse_direct_init(this, coef_fine, dm_fine, dm_crs, null_space)
    class(pnpn2_coarse_direct_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef_fine
    type(dofmap_t), intent(inout), target :: dm_fine
    type(dofmap_t), intent(inout), target :: dm_crs
    logical, intent(in) :: null_space
    integer :: i, e, iloc, jloc, ncoef
    integer :: nelv, local_nnz
    integer :: row, col, pos
    integer, allocatable :: local_rows(:), local_cols(:)
    integer(kind=i8), allocatable :: local_ids(:), local_keys(:)
    real(kind=rp), allocatable :: local_vals(:)
    real(kind=rp), allocatable :: h1_save(:), h2_save(:)
    real(kind=rp) :: val
    logical :: ifh2_save

    call this%free()

    nelv = coef_fine%msh%nelv
    this%n_local = dm_crs%size()
    this%lx_fine = coef_fine%Xh%lx
    this%lz_fine = coef_fine%Xh%lz
    this%nloc_crs = dm_crs%Xh%lx * dm_crs%Xh%ly * dm_crs%Xh%lz
    this%null_space = null_space
    call this%xxt%init_tree(pe_rank, pe_size, null_space)

    allocate(local_ids(this%n_local))
    do i = 1, this%n_local
      local_ids(i) = dm_crs%dof(i,1,1,1)
    end do

    allocate(this%basis(this%lx_fine, this%lx_fine, this%lz_fine, this%nloc_crs))
    call generate_crs_basis(this%basis, coef_fine%Xh)
    call this%u%init(dm_fine, 'pnpn2_coarse_direct_u')
    call this%w%init(dm_fine, 'pnpn2_coarse_direct_w')
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ncoef = coef_fine%dof%size()
    allocate(h1_save(ncoef), h2_save(ncoef))
    h1_save = coef_fine%h1(1:ncoef,1,1,1)
    h2_save = coef_fine%h2(1:ncoef,1,1,1)
    ifh2_save = coef_fine%ifh2
    coef_fine%h1(1:ncoef,1,1,1) = 1.0_rp
    coef_fine%h2(1:ncoef,1,1,1) = 0.0_rp
    coef_fine%ifh2 = .false.

    allocate(local_keys(this%nloc_crs * this%nloc_crs * nelv))
    allocate(local_vals(this%nloc_crs * this%nloc_crs * nelv))
    local_nnz = 0

    do jloc = 1, this%nloc_crs
      this%u%x = 0.0_rp
      this%w%x = 0.0_rp
      do e = 1, nelv
        this%u%x(:,:,:,e) = this%basis(:,:,:,jloc)
      end do
      call this%ax%compute(this%w%x, this%u%x, coef_fine, this%u%msh, this%u%Xh)

      do e = 1, nelv
        col = local_crs_id(dm_crs%Xh, jloc, e)
        do iloc = 1, this%nloc_crs
          row = local_crs_id(dm_crs%Xh, iloc, e)
          val = sum(this%basis(:,:,:,iloc) * this%w%x(:,:,:,e))
          if (abs(val) .gt. 32.0_rp * epsilon(1.0_rp)) then
            local_nnz = local_nnz + 1
            local_keys(local_nnz) = encode_key(row, col, this%n_local)
            local_vals(local_nnz) = val
          end if
        end do
      end do
    end do

    coef_fine%h1(1:ncoef,1,1,1) = h1_save
    coef_fine%h2(1:ncoef,1,1,1) = h2_save
    coef_fine%ifh2 = ifh2_save
    deallocate(h1_save, h2_save)

    if (local_nnz .gt. 1) then
      call sort_key_value(local_keys, local_vals, 1, local_nnz)
    end if
    call compress_sorted_entries(local_keys, local_vals, local_nnz)
    this%nnz = local_nnz
    allocate(local_rows(this%nnz), local_cols(this%nnz))
    do i = 1, this%nnz
      local_rows(i) = decode_row(local_keys(i), this%n_local)
      local_cols(i) = decode_col(local_keys(i), this%n_local)
    end do
    call this%xxt%setup(local_ids, this%nnz, local_rows, local_cols, local_vals)
    this%initialized = .true.

    deallocate(local_ids, local_keys, local_vals, local_rows, local_cols)
  end subroutine pnpn2_coarse_direct_init

  !> Release all coarse-grid storage.
  subroutine pnpn2_coarse_direct_free(this)
    class(pnpn2_coarse_direct_t), intent(inout) :: this

    if (allocated(this%gid_to_idx)) deallocate(this%gid_to_idx)
    if (allocated(this%local_to_idx)) deallocate(this%local_to_idx)
    if (allocated(this%row_ptr)) deallocate(this%row_ptr)
    if (allocated(this%col_ind)) deallocate(this%col_ind)
    if (allocated(this%diag_inv)) deallocate(this%diag_inv)
    if (allocated(this%vals)) deallocate(this%vals)
    if (allocated(this%share_weight)) deallocate(this%share_weight)
    if (allocated(this%rhs)) deallocate(this%rhs)
    if (allocated(this%sol)) deallocate(this%sol)
    if (allocated(this%res)) deallocate(this%res)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%ap)) deallocate(this%ap)
    if (allocated(this%z)) deallocate(this%z)
    if (allocated(this%basis)) deallocate(this%basis)
    call this%xxt%free()
    call this%u%free()
    call this%w%free()
    if (allocated(this%ax)) deallocate(this%ax)

    this%n_local = 0
    this%n_global = 0
    this%max_gid = 0
    this%lx_fine = 0
    this%lz_fine = 0
    this%nloc_crs = 0
    this%nnz = 0
    this%n_local_unique = 0
    this%null_space = .false.
    this%initialized = .false.
  end subroutine pnpn2_coarse_direct_free

  !> Solve the condensed coarse system and expand it back to local dofs.
  subroutine pnpn2_coarse_direct_solve(this, e, r, coef)
    class(pnpn2_coarse_direct_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: e(this%n_local)
    real(kind=rp), intent(in) :: r(this%n_local)
    type(coef_t), intent(in) :: coef
    integer :: ierr, i, idx

    if (.not. this%initialized) then
      call neko_error('PnPn-2 coarse solver used before initialization.')
    end if

    call this%xxt%solve(e, r)
  end subroutine pnpn2_coarse_direct_solve

  !> Sort coupled key/value arrays in ascending key order.
  recursive subroutine sort_key_value(keys, vals, lo, hi)
    integer, intent(in) :: lo, hi
    integer(kind=i8), intent(inout) :: keys(:)
    real(kind=rp), intent(inout) :: vals(:)
    integer :: i, j
    integer(kind=i8) :: pivot, key_tmp
    real(kind=rp) :: val_tmp

    if (lo .ge. hi) return

    i = lo
    j = hi
    pivot = keys((lo + hi) / 2)
    do
      do while (keys(i) .lt. pivot)
        i = i + 1
      end do
      do while (keys(j) .gt. pivot)
        j = j - 1
      end do
      if (i .le. j) then
        key_tmp = keys(i)
        keys(i) = keys(j)
        keys(j) = key_tmp
        val_tmp = vals(i)
        vals(i) = vals(j)
        vals(j) = val_tmp
        i = i + 1
        j = j - 1
      end if
      if (i .gt. j) exit
    end do

    if (lo .lt. j) call sort_key_value(keys, vals, lo, j)
    if (i .lt. hi) call sort_key_value(keys, vals, i, hi)
  end subroutine sort_key_value

  !> Compress a sorted key/value list by summing duplicate keys in place.
  subroutine compress_sorted_entries(keys, vals, n)
    integer, intent(inout) :: n
    integer(kind=i8), intent(inout) :: keys(:)
    real(kind=rp), intent(inout) :: vals(:)
    integer :: i, out

    if (n .le. 1) return

    out = 1
    do i = 2, n
      if (keys(i) .eq. keys(out)) then
        vals(out) = vals(out) + vals(i)
      else
        out = out + 1
        keys(out) = keys(i)
        vals(out) = vals(i)
      end if
    end do
    n = out
  end subroutine compress_sorted_entries

  !> Encode a `(row, col)` pair into one sortable 64-bit key.
  integer(kind=i8) function encode_key(row, col, ncols) result(key)
    integer, intent(in) :: row, col, ncols
    key = int(row - 1, i8) * int(ncols, i8) + int(col, i8)
  end function encode_key

  !> Decode the row index from an encoded `(row, col)` key.
  integer function decode_row(key, ncols) result(row)
    integer(kind=i8), intent(in) :: key
    integer, intent(in) :: ncols
    row = int((key - 1_i8) / int(ncols, i8)) + 1
  end function decode_row

  !> Decode the column index from an encoded `(row, col)` key.
  integer function decode_col(key, ncols) result(col)
    integer(kind=i8), intent(in) :: key
    integer, intent(in) :: ncols
    col = int(mod(key - 1_i8, int(ncols, i8))) + 1
  end function decode_col

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
    do j = 1, size(basis,4)
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
