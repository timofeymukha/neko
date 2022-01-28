!> Krylov preconditioner
module hsmg
  use neko_config
  use math
  use utils
  use precon
  use ax_product
  use ax_helm_fctry  
  use krylov_fctry
  use gather_scatter
  use fast3d
  use interpolation
  use bc
  use dirichlet
  use fdm
  use schwarz
  use tensor
  use jacobi
  use sx_jacobi
  use device_jacobi
  use device
  use device_math
  implicit none
  private

  !Struct to arrange our multigridlevels
  type, private :: multigrid_t
     type(dofmap_t), pointer :: dof
     type(gs_t), pointer  :: gs_h
     type(space_t), pointer :: Xh
     type(coef_t), pointer :: coef
     type(bc_list_t), pointer :: bclst
     type(schwarz_t), pointer :: schwarz
     type(field_t), pointer :: e
  end type multigrid_t
 
  type, public, extends(pc_t) :: hsmg_t
     type(mesh_t), pointer :: msh
     integer :: nlvls !< Number of levels in the multigrid
     type(multigrid_t), allocatable :: grids(:) !< array for multigrids
     type(gs_t) :: gs_crs, gs_mg !< gather scatter for lower levels
     type(space_t) :: Xh_crs, Xh_mg !< spaces for lower levels
     type(dofmap_t) :: dm_crs, dm_mg 
     type(coef_t) :: c_crs, c_mg 
     type(dirichlet_t) :: bc_crs, bc_mg
     type(bc_list_t) :: bclst_crs, bclst_mg
     type(schwarz_t) :: schwarz, schwarz_mg, schwarz_crs !< Schwarz decompostions
     type(field_t) :: e, e_mg, e_crs !< Solve fields
     class(ksp_t), allocatable :: crs_solver !< Solver for course problem
     integer :: niter = 10 !< Number of iter of crs sovlve
     class(pc_t), allocatable :: pc_crs !< Some basic precon for crs
     class(ax_t), allocatable :: ax !< Matrix for crs solve
     real(kind=rp), allocatable :: r(:)!< Residual work array
     type(interpolator_t) :: interp_fine_mid
     type(interpolator_t) :: interp_mid_crs
     real(kind=rp), allocatable :: w(:) !< work array
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
  contains
     procedure, pass(this) :: init => hsmg_init
     procedure, pass(this) :: free => hsmg_free
     procedure, pass(this) :: solve => hsmg_solve
     procedure, pass(this) :: update => hsmg_set_h
  end type hsmg_t
  
contains
  
  !> @note I do not think we actually use the same grids as they do in the original!
  subroutine hsmg_init(this, msh, Xh, coef, dof, gs_h, bclst, crs_pctype)
    class(hsmg_t), intent(inout) :: this
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h 
    type(bc_list_t), intent(inout), target :: bclst
    character(len=*), optional :: crs_pctype
    integer :: n
    
    call this%free()
    if(Xh%lx .lt. 5) then
       call neko_error('Insufficient number of GLL points for hsmg precon. Minimum degree 4 and 5 GLL points required.')
    end if
    this%nlvls = 3 
    this%msh => msh
    allocate(this%grids(this%nlvls))
    allocate(this%w(dof%n_dofs))
    allocate(this%r(dof%n_dofs))

    
    ! Compute all elements as if they are deformed
    call mesh_all_deformed(msh)

    n = dof%n_dofs
    call field_init(this%e, dof,'work array')
    
    call space_init(this%Xh_crs, GLL, 2, 2, 2)
    this%dm_crs = dofmap_t(msh, this%Xh_crs) 
    call gs_init(this%gs_crs, this%dm_crs)
    call field_init(this%e_crs, this%dm_crs,'work crs')
    call coef_init(this%c_crs, this%gs_crs)
    
    ! Create backend specific Ax operator
    call ax_helm_factory(this%ax)


    ! Create a backend specific preconditioner
    ! Note we can't call the pc factory since hsmg is a pc...
    if (NEKO_BCKND_SX .eq. 1) then
       allocate(sx_jacobi_t::this%pc_crs)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       allocate(jacobi_t::this%pc_crs)
    else if (NEKO_BCKND_CUDA .eq. 1 .or. NEKO_BCKND_HIP .eq. 1) then
       allocate(device_jacobi_t::this%pc_crs)
    else
       allocate(jacobi_t::this%pc_crs)
    end if

    ! Create a backend specific krylov solver
    if (present(crs_pctype)) then
       call krylov_solver_factory(this%crs_solver, &
            this%dm_crs%n_dofs, trim(crs_pctype), M = this%pc_crs)
    else
       call krylov_solver_factory(this%crs_solver, &
            this%dm_crs%n_dofs, 'cg', M = this%pc_crs)
    end if

    call this%bc_crs%init(this%dm_crs)
    call this%bc_crs%mark_zone(msh%outlet)
    call this%bc_crs%finalize()
    call this%bc_crs%set_g(real(0d0,rp))
    call bc_list_init(this%bclst_crs)
    call bc_list_add(this%bclst_crs, this%bc_crs)

    call space_init(this%Xh_mg, GLL, 4, 4, 4)
    this%dm_mg = dofmap_t(msh, this%Xh_mg) 
    call gs_init(this%gs_mg, this%dm_mg)
    call field_init(this%e_mg, this%dm_mg,'work midl')
    call coef_init(this%c_mg, this%gs_mg)
    
    call this%bc_mg%init(this%dm_mg)
    call this%bc_mg%mark_zone(msh%outlet)
    call this%bc_mg%finalize()
    call this%bc_mg%set_g(0.0_rp)
    call bc_list_init(this%bclst_mg)
    call bc_list_add(this%bclst_mg, this%bc_mg)

    call this%schwarz%init(Xh, dof, gs_h, bclst, msh)
    call this%schwarz_mg%init(this%Xh_mg, this%dm_mg, this%gs_mg,&
                              this%bclst_mg, msh)

    call this%interp_fine_mid%init(Xh,this%Xh_mg)
    call this%interp_mid_crs%init(this%Xh_mg,this%Xh_crs)

    call hsmg_fill_grid(dof, gs_h, Xh, coef, bclst, this%schwarz, &
                        this%e, this%grids, 3) 
    call hsmg_fill_grid(this%dm_mg, this%gs_mg, this%Xh_mg, this%c_mg, &
                        this%bclst_mg, this%schwarz_mg, this%e_mg, &
                        this%grids, 2) 
    call hsmg_fill_grid(this%dm_crs, this%gs_crs, this%Xh_crs, &
                        this%c_crs, this%bclst_crs, this%schwarz_crs, &
                        this%e_crs, this%grids, 1) 
                 
    call hsmg_set_h(this)
    if (NEKO_BCKND_CUDA .eq. 1 .or. NEKO_BCKND_HIP .eq. 1) then
       call device_map(this%w, this%w_d, n)
       call device_map(this%r, this%r_d, n)
    end if
    select type(pc => this%pc_crs)
    type is (jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    type is (sx_jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    type is (device_jacobi_t)
       call pc%init(this%c_crs, this%dm_crs, this%gs_crs)
    end select



  end subroutine hsmg_init
  
  subroutine hsmg_set_h(this)
   class(hsmg_t), intent(inout) :: this
!    integer :: i
   !Yeah I dont really know what to do here. For incompressible flow not much happens
    this%grids(1)%coef%ifh2 = .false.
    call copy(this%grids(1)%coef%h1, this%grids(3)%coef%h1, this%grids(1)%dof%n_dofs)
    if (NEKO_BCKND_CUDA .eq. 1 .or. NEKO_BCKND_HIP .eq. 1) then
        call device_copy(this%grids(1)%coef%h1_d,this%grids(3)%coef%h1_d,this%grids(1)%dof%n_dofs)
    end if
  end subroutine hsmg_set_h


  subroutine hsmg_fill_grid(dof, gs_h, Xh, coef, bclst, schwarz, e, grids, l) 
    type(dofmap_t), target, intent(in):: dof
    type(gs_t), target, intent(in) :: gs_h
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(bc_list_t), target, intent(in) :: bclst
    type(schwarz_t), target, intent(in) :: schwarz
    type(field_t), target, intent(in) :: e
    integer, intent(in) :: l
    type(multigrid_t), intent(inout), dimension(l) :: grids


    grids(l)%dof => dof
    grids(l)%gs_h => gs_h
    grids(l)%Xh => Xh
    grids(l)%coef => coef
    grids(l)%bclst => bclst
    grids(l)%schwarz => schwarz
    grids(l)%e => e

  end subroutine hsmg_fill_grid

  subroutine hsmg_free(this)
    class(hsmg_t), intent(inout) :: this

    if (allocated(this%ax)) then
       deallocate(this%ax)
    end if
    
    if (allocated(this%grids)) then
       deallocate(this%grids)
    end if
    
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    
    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    call this%schwarz%free()
    call this%schwarz_mg%free()
    call coef_free(this%c_crs)
    call coef_free(this%c_mg)
    call field_free(this%e)
    call field_free(this%e_mg)
    call field_free(this%e_crs)
    call gs_free(this%gs_crs)
    call gs_free(this%gs_mg)
    call this%interp_mid_crs%free()
    call this%interp_fine_mid%free()

    if (allocated(this%crs_solver)) then
       call krylov_solver_destroy(this%crs_solver)
       deallocate(this%crs_solver)
    end if
    
    if (allocated(this%pc_crs)) then 
       select type(pc => this%pc_crs)
       type is (jacobi_t)
          call pc%free()
       type is (sx_jacobi_t)
          call pc%free()
       end select       
       deallocate(this%pc_crs)
    end if
    
  end subroutine hsmg_free

  !> The h1mg preconditioner from Nek5000.
  subroutine hsmg_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(hsmg_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    type(ksp_monitor_t) :: crs_info
     
    if (NEKO_BCKND_CUDA .eq. 1 .or. NEKO_BCKND_HIP .eq. 1) then
       z_d = device_get_ptr(z,n)
       r_d = device_get_ptr(r,n)
       !We should not work with the input 
       call device_copy(this%r_d, r_d, n)

       !OVERLAPPING Schwarz exchange and solve
       call this%grids(3)%schwarz%compute(z, this%r)      
       !! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
       call device_col2(this%r_d, this%grids(3)%coef%mult_d, &
                    this%grids(3)%dof%n_dofs)
       !Restrict to middle level
       call this%interp_fine_mid%map(this%w,this%r,this%msh%nelv,this%grids(2)%Xh)
       call gs_op_vector(this%grids(2)%gs_h, this%w, &
                         this%grids(2)%dof%n_dofs, GS_OP_ADD)
       !OVERLAPPING Schwarz exchange and solve
       call this%grids(2)%schwarz%compute(this%grids(2)%e%x,this%w)  
       call device_col2(this%w_d, this%grids(2)%coef%mult_d, this%grids(2)%dof%n_dofs)
       !restrict residual to crs
       call this%interp_mid_crs%map(this%r,this%w,this%msh%nelv,this%grids(1)%Xh)
       !Crs solve

       call gs_op_vector(this%grids(1)%gs_h, this%r, &
                         this%grids(1)%dof%n_dofs, GS_OP_ADD)
       call bc_list_apply_scalar(this%grids(1)%bclst, this%r, &
                                 this%grids(1)%dof%n_dofs)
       crs_info = this%crs_solver%solve(this%Ax, this%grids(1)%e, this%r, &
                                    this%grids(1)%dof%n_dofs, &
                                    this%grids(1)%coef, &
                                    this%grids(1)%bclst, &
                                    this%grids(1)%gs_h, this%niter)          
       call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x,&
                                 this%grids(1)%dof%n_dofs)


       call this%interp_mid_crs%map(this%w,this%grids(1)%e%x,this%msh%nelv,this%grids(2)%Xh)
       call device_add2(this%grids(2)%e%x_d, this%w_d, this%grids(2)%dof%n_dofs)

       call this%interp_fine_mid%map(this%w,this%grids(2)%e%x,this%msh%nelv,this%grids(3)%Xh)
       call device_add2(z_d, this%w_d, this%grids(3)%dof%n_dofs)
       call gs_op_vector(this%grids(3)%gs_h, z, &
                         this%grids(3)%dof%n_dofs, GS_OP_ADD)
       call device_col2(z_d, this%grids(3)%coef%mult_d, this%grids(3)%dof%n_dofs)
    else
       !We should not work with the input 
       call copy(this%r, r, n)

       !OVERLAPPING Schwarz exchange and solve
       call this%grids(3)%schwarz%compute(z, this%r)      
       ! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
       call col2(this%r, this%grids(3)%coef%mult, &
                    this%grids(3)%dof%n_dofs)
       !Restrict to middle level
       call this%interp_fine_mid%map(this%w,this%r,this%msh%nelv,this%grids(2)%Xh)
       call gs_op_vector(this%grids(2)%gs_h, this%w, &
                         this%grids(2)%dof%n_dofs, GS_OP_ADD)
       !OVERLAPPING Schwarz exchange and solve
       call this%grids(2)%schwarz%compute(this%grids(2)%e%x,this%w)  
       call col2(this%w, this%grids(2)%coef%mult, this%grids(2)%dof%n_dofs)
       !restrict residual to crs
       call this%interp_mid_crs%map(this%r,this%w,this%msh%nelv,this%grids(1)%Xh)
       !Crs solve

       call gs_op_vector(this%grids(1)%gs_h, this%r, &
                         this%grids(1)%dof%n_dofs, GS_OP_ADD)
       call bc_list_apply_scalar(this%grids(1)%bclst, this%r, &
                                 this%grids(1)%dof%n_dofs)
       crs_info = this%crs_solver%solve(this%Ax, this%grids(1)%e, this%r, &
                                    this%grids(1)%dof%n_dofs, &
                                    this%grids(1)%coef, &
                                    this%grids(1)%bclst, &
                                    this%grids(1)%gs_h, this%niter)    
       call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x,&
                                 this%grids(1)%dof%n_dofs)


       call this%interp_mid_crs%map(this%w,this%grids(1)%e%x,this%msh%nelv,this%grids(2)%Xh)
       call add2(this%grids(2)%e%x, this%w, this%grids(2)%dof%n_dofs)

       call this%interp_fine_mid%map(this%w,this%grids(2)%e%x,this%msh%nelv,this%grids(3)%Xh)
       call add2(z, this%w, this%grids(3)%dof%n_dofs)
       call gs_op_vector(this%grids(3)%gs_h, z, &
                         this%grids(3)%dof%n_dofs, GS_OP_ADD)
       call col2(z, this%grids(3)%coef%mult, this%grids(3)%dof%n_dofs)
    end if
  end subroutine hsmg_solve
end module hsmg
