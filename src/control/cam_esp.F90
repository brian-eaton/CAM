module cam_esp

! Interfaces for the interaction of CAM with the ESP (external system processing) component.
! The implementation is:
!   1) CAM receives a pause signal from the coupler.
!   2) Before returning control to the coupler CAM writes its state variables (on the physics
!      grid) to a file (the cam.esp file).
!   3) On the next timestep (the resume phase) the file modified by ESP is read, and the tendencies
!      of the state variables are computed and stored in a ptend object.  The physics state
!      is then updated by a call to physics_update just as it would be for any other physics
!      parameterization.   
   

use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
use spmd_utils,       only: masterproc
use cam_instance,     only: inst_suffix
use cam_control_mod,  only: caseid

use ppgrid,           only: begchunk, endchunk, pcols, pver
use constituents,     only: pcnst
use physics_types,    only: physics_state, physics_tend, physics_ptend, &
                            physics_update, physics_ptend_init,         &
                            physics_state_alloc, physics_state_dealloc

use time_manager,     only: get_ref_date, timemgr_get_calendar_cf, get_curr_time

use filenames,        only: interpret_filename_spec

use cam_history_support, only: &
                            write_hist_coord_attrs, lookup_hist_coord_indices, &
                            write_hist_coord_vars, sec2hms, date2yyyymmdd

use cam_grid_support, only: cam_grid_id, cam_grid_header_info_t, &
                            cam_grid_write_attr, cam_grid_write_var, &
                            cam_grid_dimensions, cam_grid_get_decomp, &
                            cam_grid_get_dim_names

use cam_pio_utils,    only: cam_pio_createfile, cam_pio_def_dim, cam_pio_def_var, &
                            cam_pio_handle_error, cam_pio_openfile, cam_pio_closefile

use ncdio_atm,        only: infld

use pio,              only: pio_global, pio_unlimited, pio_offset_kind, pio_double, &
                            pio_seterrorhandling, pio_bcast_error, pio_noerr, pio_nowrite, &
                            file_desc_t, var_desc_t, io_desc_t, &
                            pio_def_dim, pio_def_var, pio_enddef, &
                            pio_put_att, pio_put_var, pio_setframe, pio_write_darray

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   cam_esp_write,  &  ! write state variables to cam.esp file
   cam_esp_resume     ! read cam.esp file and update state variables


character(len=cl) :: fname ! esp filename

integer :: grid_id    ! grid id for the physics grid
integer :: hdimcnt    ! number of dimensions for horizontal grid


! pio variable descriptors
type(var_desc_t) :: time_desc, ps_desc, t_desc, u_desc, v_desc, q_desc

!=========================================================================================
contains
!=========================================================================================

subroutine cam_esp_write(state)

   ! Arguments
   type(physics_state), intent(in), dimension(begchunk:endchunk) :: state

   ! Local workspace
   character(len=cl) :: filename_spec ! filename template for esp file
   type(file_desc_t) :: fh
   real(r8)          :: time
   integer           :: ndcur, nscur
   integer           :: ierr
   !-----------------------------------------------------------------------

   ! Set template for esp filename based on instance suffix
   ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
   filename_spec = '%c.cam' // trim(inst_suffix) //'.esp.%y-%m-%d-%s.nc'

   fname = interpret_filename_spec( filename_spec )

   call cam_pio_createfile(fh, trim(fname), 0)

   ! sets the grid id and variable descriptors in module data
   call define_file(fh)

   ! Write coordinate data:
   ! time
   call get_curr_time(ndcur, nscur)
   time = ndcur + nscur/86400._r8
   ierr = pio_put_var(fh, time_desc, (/1/), (/1/), (/time/))
   ! vertical coordinates
   call write_hist_coord_vars(fh)
   ! horizontal coordinates
   call cam_grid_write_var(fh, grid_id)

   ! write variables
   call write_state(fh, state)

   ! Close file
   call cam_pio_closefile(fh)
      
end subroutine cam_esp_write

!========================================================================================

subroutine cam_esp_resume(ztodt, phys_state, phys_tend)


   ! Arguments
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)

   ! Local workspace
   type(file_desc_t)   :: fh

   integer             :: lchnk, ncol
   logical             :: ls, lu, lv
   logical             :: lq(pcnst)
   logical             :: found

   type(physics_state) :: esp_state(begchunk:endchunk)
   type(physics_ptend) :: esp_ptend(begchunk:endchunk)

   real(r8), allocatable :: chunks2d(:,:)
   real(r8), allocatable :: chunks3d(:,:,:)
   character(len=8)    :: dim1name, dim2name
   

   character(len=*), parameter :: sub = 'cam_esp_resume'
   !---------------------------------------------------------------------------
  
   ! get filehandle pointer to esp file
   write(iulog,*) sub//': open '//trim(fname)
   call cam_pio_openfile(fh, trim(fname), pio_nowrite)

   ! for first cut hardcode the fields expected in the esp file
   ls = .true.
   lu = .true.
   lv = .true.
   lq = .false.
   lq(1) = .true.

   ! allocate storage for file data and tendencies
   do lchnk = begchunk, endchunk
      call physics_state_alloc(esp_state(lchnk), lchnk, pcols)
      call physics_ptend_init(esp_ptend(lchnk), pcols, sub, ls, lu, lv, lq)
   end do

   ! read esp data
   allocate(&
      chunks2d(pcols,begchunk:endchunk), &
      chunks3d(pcols,pver,begchunk:endchunk))

   call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

   call infld('PS', fh, dim1name, dim2name, 1, pcols, &
              begchunk, endchunk, chunks2d, found, gridname='physgrid')
   if (.not. found) call endrun(sub//': PS not found on file')
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      esp_state(lchnk)%ps(:ncol) = chunks2d(:ncol,lchnk)
   end do

   call infld('T', fh, dim1name, dim2name, 'lev', 1, pcols, 1, pver, &
              begchunk, endchunk, chunks3d, found, gridname='physgrid')
   if (.not. found) call endrun(sub//': T not found on file')
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      esp_state(lchnk)%t(:ncol,:) = chunks3d(:ncol,:,lchnk)
   end do
  
   call infld('U', fh, dim1name, dim2name, 'lev', 1, pcols, 1, pver, &
              begchunk, endchunk, chunks3d, found, gridname='physgrid')
   if (.not. found) call endrun(sub//': U not found on file')
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      esp_state(lchnk)%u(:ncol,:) = chunks3d(:ncol,:,lchnk)
   end do
  
   call infld('V', fh, dim1name, dim2name, 'lev', 1, pcols, 1, pver, &
              begchunk, endchunk, chunks3d, found, gridname='physgrid')
   if (.not. found) call endrun(sub//': V not found on file')
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      esp_state(lchnk)%v(:ncol,:) = chunks3d(:ncol,:,lchnk)
   end do
  
   call infld('Q', fh, dim1name, dim2name, 'lev', 1, pcols, 1, pver, &
              begchunk, endchunk, chunks3d, found, gridname='physgrid')
   if (.not. found) call endrun(sub//': Q not found on file')
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      esp_state(lchnk)%q(:ncol,:,1) = chunks3d(:ncol,:,lchnk)
   end do
  
   call cam_pio_closefile(fh)

   ! debug check
   do lchnk = begchunk, endchunk
      ncol = esp_state(lchnk)%ncol
      if (any(esp_state(lchnk)%ps(:ncol) /= phys_state(lchnk)%ps(:ncol))) &
         write(iulog,*) sub//': diffs in PS'
      if (any(esp_state(lchnk)%t(:ncol,:) /= phys_state(lchnk)%t(:ncol,:))) &
         write(iulog,*) sub//': diffs in T'
      if (any(esp_state(lchnk)%u(:ncol,:) /= phys_state(lchnk)%u(:ncol,:))) &
         write(iulog,*) sub//': diffs in U'
      if (any(esp_state(lchnk)%v(:ncol,:) /= phys_state(lchnk)%v(:ncol,:))) &
         write(iulog,*) sub//': diffs in V'
      if (any(esp_state(lchnk)%q(:ncol,:,1) /= phys_state(lchnk)%q(:ncol,:,1))) &
         write(iulog,*) sub//': diffs in Q'
   end do

   ! compute tendencies.  note the temperature tendency is converted to dry
   ! static energy tendency.




   ! update the pressure state variables since this is not done by physics_update


   ! update phys_state



   ! deallocate storage for file data.  ptend already deallocated by phys_update.
   do lchnk = begchunk, endchunk
      call physics_state_dealloc(esp_state(lchnk))
   end do



end subroutine cam_esp_resume

!=========================================================================================
! Private routines
!=========================================================================================

subroutine define_file(fh)

   ! Define dimensions, variables, attributes of rst file.

   ! arguments
   type(file_desc_t), intent(inout) :: fh

   ! local variables
   integer :: i
   integer :: yr, mon, day, nbsec, nbdate
   integer :: time_dimid
   integer :: lev_dimid
   integer :: dimids_2d(3)  ! dimension ids for 2D field w/ time dimension
   integer :: dimids_3d(4)  ! dimension ids for 3D field w/ time dimension
   integer, allocatable :: mdimids(:) ! dim IDs for "middle" dimensions
   integer              :: lev_idx(1)
   character(len=cl)    :: str        ! character temporary

   type(cam_grid_header_info_t) :: hdr

   integer :: ierr, err_handling
   !---------------------------------------------------------------------------

   ! time coordinate
   call cam_pio_def_dim(fh, 'time', pio_unlimited, time_dimid)

   call get_ref_date(yr, mon, day, nbsec)
   nbdate = yr*10000 + mon*100 + day
   str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
   call define_var(fh, 'time', (/time_dimid/), 'time', trim(str), time_desc)
   str = timemgr_get_calendar_cf()
   ierr = pio_put_att (fh, time_desc, 'calendar', trim(str))

   ! vertical coordinate
   call write_hist_coord_attrs(fh, -1, mdimids)
   
   ! find the dimid for lev
   call lookup_hist_coord_indices((/'lev'/), lev_idx)
   lev_dimid = mdimids(lev_idx(1))

   ! horizontal grid
   grid_id = cam_grid_id('physgrid')
   call cam_grid_write_attr(fh, grid_id, hdr)

   ! get dimension ids for the horizontal grid
   ! hdimcnt is 1 for unstructured grid and 2 for rectangular grid
   hdimcnt = hdr%num_hdims()

   ! set the dimension arrays for 2D and 3D time varying fields
   do i = 1, hdimcnt
      dimids_2d(i) = hdr%get_hdimid(i)
      dimids_3d(i) = hdr%get_hdimid(i)
   end do
   dimids_2d(hdimcnt+1) = time_dimid

   dimids_3d(hdimcnt+1) = lev_dimid
   dimids_3d(hdimcnt+2) = time_dimid

   ! field variables
   call define_var(fh, 'PS', dimids_2d(:hdimcnt+1), 'Surface pressure', 'Pa', ps_desc)
   call define_var(fh, 'T',  dimids_3d(:hdimcnt+2), 'Temperature', 'K', t_desc)
   call define_var(fh, 'U',  dimids_3d(:hdimcnt+2), 'Zonal wind', 'm/s', u_desc)
   call define_var(fh, 'V',  dimids_3d(:hdimcnt+2), 'Meridional wind', 'm/s', v_desc)
   call define_var(fh, 'Q',  dimids_3d(:hdimcnt+2), 'Specific humidity', 'kg/kg', q_desc)

   ! global attribute
   ierr = pio_put_att(fh, pio_global, 'caseid', trim(caseid))

   ! End of file definition mode.
   ierr = pio_enddef(fh)

end subroutine define_file

!========================================================================================

subroutine define_var(fh, name, dimids, long_name, units, vdesc)

   ! Define a variable and its attributes

   ! arguments
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: name       ! variable name
   integer,           intent(in)    :: dimids(:)  ! pio dimension IDs of variable
   character(len=*),  intent(in)    :: long_name  ! long name of variable
   character(len=*),  intent(in)    :: units      ! units of variable
   type(var_desc_t),  intent(out)   :: vdesc      ! pio varible descriptor

   ! local variables
   integer :: ierr, err_handling
   character(len=*), parameter :: sub = 'define_var'
   !---------------------------------------------------------------------------

   call cam_pio_def_var(fh, trim(name), pio_double, dimids, vdesc)

   call pio_seterrorhandling(fh, pio_bcast_error, err_handling)

   ierr = pio_put_att(fh, vdesc, 'long_name', trim(long_name))
   call cam_pio_handle_error(ierr, sub//' cannot define long_name for '//trim(name))

   ierr = pio_put_att(fh, vdesc, 'units', trim(units))
   call cam_pio_handle_error(ierr, sub//' cannot define units for '//trim(name))


   call pio_seterrorhandling(fh, err_handling)

end subroutine define_var

!========================================================================================

subroutine write_state(fh, state)

   ! Write the state variables

   ! arguments
   type(file_desc_t),   intent(inout) :: fh
   type(physics_state), intent(in)    :: state(begchunk:endchunk)

   ! local variables
   integer :: i, ncol
   integer :: gdims(3), dims(3)

   type(io_desc_t), pointer :: iodesc

   real(r8), allocatable :: buf2d(:,:), buf3d(:,:,:)

   integer :: ierr, err_handling
   character(len=*), parameter :: sub = 'write_state'
   !---------------------------------------------------------------------------


   ! get horizontal grid dimensions (global)
   call cam_grid_dimensions(grid_id, gdims(1:2))

   ! add dimension for 3D field
   gdims(hdimcnt+1) = pver

   ! write 2D fields

   ! dimensions of local data for 2D field
   dims(1) = pcols
   dims(2) = endchunk - begchunk + 1

   ! PIO decomposition for 2D field
   call cam_grid_get_decomp(grid_id, dims(1:2), gdims(1:hdimcnt), &
                            pio_double, iodesc)

   ! make error handling local to give better diagnostics
   call pio_seterrorhandling(fh, pio_bcast_error, err_handling)

   allocate(buf2d(pcols,begchunk:endchunk))

   ! PS
   do i = begchunk, endchunk
      ncol = state(i)%ncol
      buf2d(:ncol,i) = state(i)%ps(:ncol)
   end do
   call pio_setframe(fh, ps_desc, int(1,kind=pio_offset_kind))
   call pio_write_darray(fh, ps_desc, iodesc, buf2d, ierr)
   call cam_pio_handle_error(ierr, sub//': writing PS')

   nullify(iodesc)
   deallocate(buf2d)

   ! write 3D fields

   ! dimensions of local data for 3D field
   dims(1) = pcols
   dims(2) = pver
   dims(3) = endchunk - begchunk + 1

   ! PIO decomposition for 3D field
   call cam_grid_get_decomp(grid_id, dims, gdims(1:hdimcnt+1), &
                            pio_double, iodesc)

   allocate(buf3d(pcols,pver,begchunk:endchunk))

   ! T
   do i = begchunk, endchunk
      ncol = state(i)%ncol
      buf3d(:ncol,:,i) = state(i)%t(:ncol,:)
   end do
   call pio_setframe(fh, t_desc, int(1,kind=pio_offset_kind))
   call pio_write_darray(fh, t_desc, iodesc, buf3d, ierr)
   call cam_pio_handle_error(ierr, sub//': writing T')

   ! U
   do i = begchunk, endchunk
      ncol = state(i)%ncol
      buf3d(:ncol,:,i) = state(i)%u(:ncol,:)
   end do
   call pio_setframe(fh, u_desc, int(1,kind=pio_offset_kind))
   call pio_write_darray(fh, u_desc, iodesc, buf3d, ierr)
   call cam_pio_handle_error(ierr, sub//': writing U')

   ! V
   do i = begchunk, endchunk
      ncol = state(i)%ncol
      buf3d(:ncol,:,i) = state(i)%v(:ncol,:)
   end do
   call pio_setframe(fh, v_desc, int(1,kind=pio_offset_kind))
   call pio_write_darray(fh, v_desc, iodesc, buf3d, ierr)
   call cam_pio_handle_error(ierr, sub//': writing V')

   ! Q
   do i = begchunk, endchunk
      ncol = state(i)%ncol
      buf3d(:ncol,:,i) = state(i)%q(:ncol,:,1)
   end do
   call pio_setframe(fh, q_desc, int(1,kind=pio_offset_kind))
   call pio_write_darray(fh, q_desc, iodesc, buf3d, ierr)
   call cam_pio_handle_error(ierr, sub//': writing Q')

   nullify(iodesc)
   deallocate(buf3d)

   call pio_seterrorhandling(fh, err_handling)

end subroutine write_state

!========================================================================================

end module cam_esp
