!> Krylov preconditioner
module identity
  use math
  use utils
  use precon
  use ax_product
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type, public, extends(pc_t) :: ident_t
  contains
     procedure, pass(this) :: solve => ident_solve
  end type ident_t

contains
  !> The (default) naive preconditioner \f$ I z = r \f$
  subroutine ident_solve(this, z, r, n)
       integer, intent(inout) :: n
       class(ident_t), intent(inout) :: this
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), dimension(n), intent(inout) :: r
        call copy(z, r, n)    
  end subroutine ident_solve

end module identity
