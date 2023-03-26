!> AUTHOR : Ala Bahrami
!> DATE CREATION : 01-25-2023
!> DATES MODIFICATIONS : -
!> DESCRIPTION : The purpose of this module is to write CLASSW variables in a NetCDF format  
!> per time step which are passed to the WATROF subroutine inside the CLASSW module for 4 subarea 
!> categories (CS, GS, C, G) 
!>               
!>
!> source      : 
!> see also    : nc_io
module nc_io_file

#ifdef NETCDF
    use netcdf
	use check_stat
	!use nc_io_watdrn
	!use print_routines
#endif

	implicit none
	
	contains
	
	! Initialize the file  
	subroutine nc4_init_file(FILE_NAME, ILG, level, &
							 ncid, &
							 thlq_varid, time_varid)
	
		!> Input variables.
		character(len = *), intent(in) :: FILE_NAME
		integer, intent(in) :: ILG
		integer, intent(in) :: level
		!integer, intent(in) :: dtype

		!> Output variables.
		integer, intent(out) :: ncid, time_varid, thlq_varid 
		
		! netCDF dimensions.
		integer, parameter :: NDIMS  = 3
		!integer :: NTILES  = ILG
		integer :: tile_dimid, level_dimid, time_dimid

		real :: tiles(ILG)
		integer :: tile_varid, level_varid
		integer :: dimids(NDIMS)

		! We will write THLQ,level and time fields. 
		! Labels  
		character (len = *), parameter :: TILE_NAME  = "tile"
		!character (len = *), parameter :: THLQ_NAME  ="soil_moisture"
		character (len = *), parameter :: THLQ_NAME  ="interflow"
		character (len = *), parameter :: Level_NAME ="level"
		character (len = *), parameter :: Time_NAME  ="time"
		
		! It's good practice for each variable to carry a "units" attribute.
		character (len = *), parameter :: UNITS       = "units"
		character (len = *), parameter :: TILE_UNITS  = "-"
		!character (len = *), parameter :: Time_UNITS  = "minutes"
		character (len = *), parameter :: Time_UNITS  = "hours"
		character (len = *), parameter :: Level_UNITS = "-"
		character (len = *), parameter :: THLQ_UNITS  = "mm"	

		!> Input variables (optional).
		logical :: quiet
		
		! Local variables
		integer :: k

		! create the tile variable 
		do k = 1, ILG
			tiles(k) =  k
		end do
	  
		! open the file 
		!call nc4_open_output_file(file_path, quiet, iun, ierr)
		
		! create netCDF file and enter define mode
		call check(nf90_create(FILE_NAME, NF90_CLOBBER, ncid))
																									
		! define netCDF dimensions
		call check(nf90_def_dim(ncid, TILE_NAME, ILG, tile_dimid))    
		call check(nf90_def_dim(ncid, LEVEL_NAME, level, level_dimid))    
		call check(nf90_def_dim(ncid, Time_NAME, NF90_UNLIMITED, time_dimid))
				
		! define NetCDF variables (time, tile, level , THLQCS)
		call check(nf90_def_var(ncid, Time_NAME, NF90_DOUBLE, time_dimid, time_varid)) ! NF90_DOUBLE
	    call check(nf90_def_var(ncid, LEVEL_NAME, NF90_DOUBLE, level_dimid, level_varid)) ! NF90_DOUBLE
		call check(nf90_def_var(ncid, TILE_NAME, NF90_REAL, tile_dimid, tile_varid) ) ! NF90_REAL
	  
		dimids = (/ tile_dimid, level_dimid, time_dimid/)
		call check(nf90_def_var(ncid, THLQ_NAME, NF90_DOUBLE, dimids, thlq_varid)) ! NF90_DOUBLE	

		! define attributes for time 
		call check(nf90_put_att(ncid, time_varid, UNITS, Time_UNITS))
		call check(nf90_put_att(ncid, tile_varid, UNITS, TILE_UNITS))
		call check(nf90_put_att(ncid, level_varid, UNITS, Level_UNITS))
		call check(nf90_put_att(ncid, thlq_varid, UNITS, THLQ_UNITS))
			 
		! End define mode.
		call check(nf90_enddef(ncid))
	  
		! write netcdf variable 
		call check(nf90_put_var(ncid, tile_varid, tiles))
		
	end subroutine

	!> Add data per timestep 
	subroutine nc4_add_data_file(ncid, ILG, level, &
								 THLQ_varid, time_varid, &
								 THLQCS, count)
		integer, intent(in) :: ncid
		integer, intent(in) :: ILG
		integer, intent(in) :: level
		integer, intent(in) :: THLQ_varid
		integer, intent(in) :: time_varid
		real, intent(in)    :: THLQCS(ILG, level)
		integer, intent(in) :: count
		

		!write (*,*) 'The count number is :', count  
		
		! write time
		call check(nf90_put_var(ncid, time_varid, (count), &
				   start = (/count/)))

		! write soil moisture data 	
		call check(nf90_put_var(ncid, THLQ_varid, THLQCS, &
				   start = (/1,1,count/)))
	
	end subroutine 
	
	!> close the input file  
	subroutine nc4_close_fille(ncid)
		integer, intent(in) :: ncid
		! Close the file.
		call check( nf90_close(ncid))
	end subroutine
	
end module 