!> NEKTON session data reader
!! @details This module is used to read NEKTON session data in ascii
module rea_file
  use generic_file
  use num_types
  use utils
  use mesh
  use point 
  use map
  use rea
  use re2_file
  use map_file
  use comm
  use datadist
  use htable
  implicit none
  private

  !> Interface for NEKTON ascii files
  type, public, extends(generic_file_t) :: rea_file_t
   contains
     procedure :: read => rea_file_read
     procedure :: write => rea_file_write
  end type rea_file_t

contains

  !> Load NEKTON session data from an ascii file
  subroutine rea_file_read(this, data)
    class(rea_file_t) :: this
    class(*), target, intent(inout) :: data
    type(mesh_t), pointer :: msh
    real(kind=dp), pointer :: params(:)
    character(len=3), pointer :: cbc(:,:)
    character(len=1) :: chtemp
    integer :: ndim, nparam, nskip, nlogic, nbcs
    integer :: nelgs, nelgv, i, j, ierr
    integer :: el_idx, pt_idx
    logical :: read_param, read_bcs, read_map
    real(kind=dp) :: xc(8), yc(8), zc(8)
    type(point_t) :: p(8)
    type(re2_file_t) :: re2_file
    type(map_file_t) :: map_file
    character(len=80) :: re2_fname, map_fname
    integer :: start_el, end_el, nel
    type(linear_dist_t) :: dist
    type(map_t) :: nm
    type(htable_pt_t) :: htp 

    select type(data)
    type is (rea_t)
       call rea_free(data)       
       msh => data%msh
       params => data%params
       cbc => data%cbc
       read_param = .true.
       read_bcs = .true.
    type is (mesh_t)    
       msh => data
       read_param = .false.
       read_bcs = .false.
    class default
       call neko_error('Invalid output data')
    end select

    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Reading NEKTON file ", this%fname
    end if
    
    read(9, *)
    read(9, *)
    read(9, *) ndim
    read(9, *) nparam
    
    if (.not. read_param) then
       ! Skip parameters
       do i = 1, nparam
          read(9, *)
       end do
    else       
       allocate(params(nparam))
       do i = 1, nparam
          read(9, *) params(i)
       end do
    end if
    
    ! Skip passive scalars
    read(9, *) nskip
    do i = 1, nskip
       read(9, *)
    end do
    
    ! Skip logic switches
    read(9, *) nlogic
    do i = 1, nlogic
       read(9, *)
    end do
    
    ! Read mesh info
    read(9, *)
    read(9, *)
    read(9, *) nelgs,ndim, nelgv
    if (nelgs .lt. 0) then
       re2_fname = trim(this%fname(1:scan(trim(this%fname), &
            '.', back=.true.)))//'re2' 
       call re2_file%init(re2_fname)
       call re2_file%read(msh)
    else       
       if (pe_rank .eq. 0) write(*,1) ndim, nelgv
1      format(1x,'ndim = ', i1, ', nelements =', i7)

       call filename_chsuffix(this%fname, map_fname, 'map')
       inquire(file=map_fname, exist=read_map)
       if (read_map) then
          call map_init(nm, nelgv, 2**ndim)
          call map_file%init(map_fname)
          call map_file%read(nm)
       else
          if (pe_rank .eq. 0) call neko_warning('No NEKTON map file found')
       end if

       ! Use a load-balanced linear distribution
       dist = linear_dist_t(nelgv, pe_rank, pe_size, NEKO_COMM)
       nel = dist%num_local()
       start_el = dist%start_idx() + 1
       end_el = dist%end_idx() + 1

       call mesh_init(msh, ndim, dist)

       call htp%init((2**ndim) * nel, ndim)

       el_idx = 1
       pt_idx = 0
       do i = 1, nelgv
          read(9, *)
          if (ndim .eq. 2) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             if (i .ge. start_el .and. i .le. end_el) then
                do j = 1, 4
                   p(j) = point_t(xc(j), yc(j), 0d0)
                   call rea_file_add_point(htp, p(j), pt_idx)
                end do
                call mesh_add_element(msh, el_idx, p(1), p(2), p(3), p(4))
             end if
          else if (ndim .eq. 3) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             read(9, *) (zc(j),j=1,4)
             read(9, *) (xc(j),j=5,8)
             read(9, *) (yc(j),j=5,8)
             read(9, *) (zc(j),j=5,8)
             if (i .ge. start_el .and. i .le. end_el) then
                do j = 1, 8
                   p(j) = point_t(xc(j), yc(j), zc(j))
                   call rea_file_add_point(htp, p(j), pt_idx)
                end do
                call mesh_add_element(msh, el_idx, &
                     p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
             end if
          end if
          if (i .ge. start_el .and. i .le. end_el) then
             el_idx = el_idx + 1
          end if
       end do

       call htp%free()
       
       !> @todo Add support for curved side data
       read(9, *) 
       read(9, *) nskip
       do i = 1, nskip
          read(9, *)
       end do

       ! Read fluid boundary conditions
       read(9,*) 
       read(9,*) 
       if (.not. read_bcs) then
          do i = 1, nelgv
             read(9, *)
          end do
       else 
          allocate(cbc(6,nelgv))
          do i = 1, nelgv
             do j = 1, 2*ndim
                read(9,'(a1, a3)') chtemp, cbc(j, i)
             end do
          end do
       end if

       if (pe_rank .eq. 0) write(*,*) 'Done'       
       close(9)
    endif
    
  end subroutine rea_file_read

  subroutine rea_file_write(this, data)
    class(rea_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
  end subroutine rea_file_write

  subroutine rea_file_add_point(htp, p, idx)
    type(htable_pt_t), intent(inout) :: htp
    type(point_t), intent(inout) :: p
    integer, intent(inout) :: idx
    integer :: tmp
    
    if (htp%get(p, tmp) .gt. 0) then
       idx = idx + 1
       call htp%set(p, idx)
       call p%set_id(idx)
    else
       call p%set_id(tmp)
    end if
    
  end subroutine rea_file_add_point

end module rea_file