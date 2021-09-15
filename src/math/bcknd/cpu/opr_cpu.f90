!> Operators CPU backend
module opr_cpu
  use num_types
  use mxm_wrapper
  use space
  use coefs
  use math
  use mesh
  use field
  use gather_scatter
  use mathops
  implicit none
  
contains

  subroutine opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(in) ::  u, dr, ds, dt
    real(kind=rp) :: drst(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz)
    type(space_t), pointer :: Xh 
    type(mesh_t), pointer :: msh
    integer :: e, k, lxy, lyz, lxyz
    Xh => coef%Xh
    msh => coef%msh 
    lxy  = Xh%lx*Xh%ly
    lyz  = Xh%ly*Xh%lz
    lxyz = Xh%lx*Xh%ly*Xh%lz

    do e = 1,msh%nelv
       if (msh%gdim .eq. 2) then
          call mxm(Xh%dx, Xh%lx, u(1,1,1,e), Xh%lx, du(1,1,1,e), lyz)
          call col2(du(1,1,1,e), dr(1,1,1,e), lxyz)
          call mxm(u(1,1,1,e), Xh%lx, Xh%dyt, Xh%ly, drst,Xh%ly)
          call addcol3(du(1,1,1,e), drst, ds(1,1,1,e), lxyz)
       else
          call mxm(Xh%dx, Xh%lx, u(1,1,1,e), Xh%lx, du(1,1,1,e), lyz)
          call col2(du(1,1,1,e), dr(1,1,1,e), lxyz)
          do k = 1,Xh%lz
             call mxm(u(1,1,k,e), Xh%lx, Xh%dyt, Xh%ly, drst(1,1,k), Xh%ly)
          end do
          call addcol3(du(1,1,1,e), drst, ds(1,1,1,e), lxyz)
          call mxm(u(1,1,1,e), lxy, Xh%dzt, Xh%lz, drst, Xh%lz)
          call addcol3(du(1,1,1,e), drst, dt(1,1,1,e), lxyz)
       end if
    end do
    call col2(du, coef%jacinv, coef%dof%n_dofs)    
  end subroutine opr_cpu_dudxyz

  subroutine opr_cpu_opgrad(ux, uy, uz, u, coef) 
    type(coef_t), intent(in) :: coef  
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: u
    real(kind=rp) :: ur(coef%Xh%lxyz)
    real(kind=rp) :: us(coef%Xh%lxyz)
    real(kind=rp) :: ut(coef%Xh%lxyz)
    integer :: e, i, N

    N = coef%Xh%lx - 1
    do e = 1,coef%msh%nelv
       if(coef%msh%gdim .eq. 3) then
          call local_grad3(ur, us, ut, u(1,e), N, coef%Xh%dx, coef%Xh%dxt)
          do i = 1,coef%Xh%lxyz
             ux(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdx(i,1,1,e) &
                  + us(i)*coef%dsdx(i,1,1,e) &
                  + ut(i)*coef%dtdx(i,1,1,e) )
             uy(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdy(i,1,1,e) &
                  + us(i)*coef%dsdy(i,1,1,e) &
                  + ut(i)*coef%dtdy(i,1,1,e) )
             uz(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdz(i,1,1,e) &
                  + us(i)*coef%dsdz(i,1,1,e) &
                  + ut(i)*coef%dtdz(i,1,1,e) )
          end do
       else

          call local_grad2(ur, us, u(1,e), N, coef%Xh%dx, coef%Xh%dyt)

          do i = 1,coef%Xh%lxyz
             ux(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdx(i,1,1,e) &
                  + us(i)*coef%dsdx(i,1,1,e) )
             uy(i,e) = coef%Xh%w3(i,1,1)*(ur(i)*coef%drdy(i,1,1,e) &
                  + us(i)*coef%dsdy(i,1,1,e) )
          end do
       end if
    end do
  end subroutine opr_cpu_opgrad

  subroutine local_grad3(ur, us, ut, u, n, D, Dt)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(in) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(in) :: D(0:n, 0:n)
    real(kind=rp), intent(in) :: Dt(0:n, 0:n)
    integer :: m1, m2, k

    m1 = n + 1
    m2 = m1*m1

    call mxm(D ,m1,u,m1,ur,m2)
    do k = 0,n
       call mxm(u(0, 0, k), m1, Dt, m1, us(0,0,k), m1)
    end do
    call mxm(u, m2, Dt, m1, ut, m1)

  end subroutine local_grad3

  subroutine local_grad2(ur, us, u, n, D, Dt)
    integer, intent(in) :: n    
    real(kind=rp), intent(inout) :: ur(0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n)
    real(kind=rp), intent(in) :: u(0:n, 0:n)
    real(kind=rp), intent(in) :: D(0:n, 0:n)
    real(kind=rp), intent(in) :: Dt(0:n, 0:n)
    integer :: m1, m2, k

    m1 = n + 1

    call mxm(D, m1, u, m1, ur, m1)
    call mxm(u, m1, Dt, m1, us, m1)

  end subroutine local_grad2

  subroutine opr_cpu_cdtp(dtx, x, dr, ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt
    real(kind=rp) :: wx(coef%Xh%lxyz)
    real(kind=rp) :: ta1(coef%Xh%lxyz)
    real(kind=rp) :: ta2(coef%Xh%lxyz)
    real(kind=rp) :: ta3(coef%Xh%lxyz)
    integer :: e, i1, i2, n1, n2, iz
    type(space_t), pointer :: Xh 

    Xh => coef%Xh
    n1 = Xh%lx*Xh%ly
    n2 = Xh%lx*Xh%ly

    do e = 1,coef%msh%nelv
       call col3(wx, coef%B(1,1,1,e), x(1,e), Xh%lxyz)
       call invcol2(wx, coef%jac(1,1,1,e), Xh%lxyz)
       call col3(ta1, wx, dr(1,e), Xh%lxyz)
       call mxm(Xh%dxt, Xh%lx, ta1, Xh%lx, dtx(1,e), Xh%lyz)
       call col3(ta1, wx, ds(1,e), Xh%lxyz)
       i1 = 1
       i2 = 1
       do iz = 1,Xh%lz
          call mxm(ta1(i2), Xh%lx, Xh%dy, Xh%ly,ta2(i1), Xh%ly)
          i1 = i1 + n1
          i2 = i2 + n2
       end do
       call add2(dtx(1,e), ta2, Xh%lxyz)
       call col3(ta1, wx, dt(1,e), Xh%lxyz)
       call mxm(ta1, Xh%lxy, Xh%dz, Xh%lz, ta2, Xh%lz)
       call add2(dtx(1,e), ta2, Xh%lxyz)
    end do
  end subroutine opr_cpu_cdtp

  subroutine opr_cpu_conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz
    !   Store the inverse jacobian to speed this operation up
    real(kind=rp), dimension(Xh%lx,Xh%ly,Xh%lz) :: dudr, duds, dudt
    integer :: ie, iz, i
    !   Compute vel.grad(u)
    do ie = 1,nelv
       if (gdim .eq. 3) then
          call mxm(Xh%dx, Xh%lx, u(1,1,1,ie), Xh%lx, dudr, Xh%lxy)
          do iz = 1,Xh%lz
             call mxm(u(1,1,iz,ie), Xh%lx, Xh%dyt, Xh%ly, duds(1,1,iz), Xh%ly)
          end do
          call mxm(u(1,1,1,ie), Xh%lxy, Xh%dzt, Xh%lz, dudt, Xh%lz)
          do i = 1,Xh%lxyz
             du(i,ie) = coef%jacinv(i,1,1,ie)*( &
                  vx(i,1,1,ie)*( &
                  coef%drdx(i,1,1,ie)*dudr(i,1,1) &
                  + coef%dsdx(i,1,1,ie)*duds(i,1,1) &
                  + coef%dtdx(i,1,1,ie)*dudt(i,1,1)) &
                  + vy(i,1,1,ie)*( &
                  coef%drdy(i,1,1,ie)*dudr(i,1,1) &
                  + coef%dsdy(i,1,1,ie)*duds(i,1,1) &
                  + coef%dtdy(i,1,1,ie)*dudt(i,1,1)) &
                  + vz(i,1,1,ie)*( &
                  coef%drdz(i,1,1,ie)*dudr(i,1,1) &
                  + coef%dsdz(i,1,1,ie)*duds(i,1,1) &
                  + coef%dtdz(i,1,1,ie)*dudt(i,1,1)))
          end do
       else
          !        2D
          call mxm(Xh%dx, Xh%lx, u(1,1,1,ie), Xh%lx, dudr, Xh%lyz)
          call mxm(u(1,1,1,ie), Xh%lx, Xh%dyt, Xh%ly, duds, Xh%ly)
          do i = 1,Xh%lxyz
             du(i,ie) = coef%jacinv(i,1,1,ie)*( &
                  vx(i,1,1,ie)*( &
                  coef%drdx(i,1,1,ie)*dudr(i,1,1) &
                  + coef%dsdx(i,1,1,ie)*duds(i,1,1)) &
                  + vy(i,1,1,ie)*( &
                  coef%drdy(i,1,1,ie)*dudr(i,1,1) &
                  + coef%dsdy(i,1,1,ie)*duds(i,1,1)))
          end do
       end if
    end do
  end subroutine opr_cpu_conv1

  subroutine opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call opr_cpu_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    end if
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    end if
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_cpu_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_cpu_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine opr_cpu_curl



end module opr_cpu