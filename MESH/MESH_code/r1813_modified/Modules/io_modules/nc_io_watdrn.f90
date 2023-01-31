!> AUTHOR : Ala Bahrami
!> DATE CREATION : 01-12-2023
!> DATES MODIFICATIONS : -
!> DESCRIPTION : The purpose of this module is to write model inputs which  
!> which are passed to the WATDRN subroutine inside the CLASSW module 
!>               
!>
!> source : https://docs.unidata.ucar.edu/netcdf-fortran/current/f90-variables.html
!>        : nc_io 
module nc_io_watdrn
	
#ifdef NETCDF
    use netcdf
	use check_stat
	use print_routines
    !use parse_utilities
    !use mesh_io_options
    !use model_dates
	!use typesizes
	!use strings
#endif

	implicit none
	
	contains

	!> Open ouput file  
	subroutine  nc4_open_output_file(fpath, quiet, iun, ierr)
			
	    !> Input variables.
        character(len = *), intent(in) :: fpath

        !> Input variables (optional).
        logical, intent(in), optional :: quiet

        !> Output variables.
        integer, intent(out) :: iun, ierr

        !> Local variables.
        character(len = DEFAULT_LINE_LENGTH) line
        character(len = DEFAULT_FIELD_LENGTH) code
        logical :: q = .false.

        !> Print a message (if not quiet).
        if (present(quiet)) q = quiet
        if (.not. q) then
            call reset_tab()
            call print_message("OPENING: " // trim(fpath) // " (for output)")
            call increase_tab()
        end if 	
		
		!> Open the file with write access.
        !ierr = nf90_create(fpath, NF90_NETCDF4, iun)
        ierr = nf90_create(fpath, NF90_CLOBBER, iun)
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            if (q) then
                call print_error("An error occurred opening the file (Code: " // trim(adjustl(code)) // "): " // trim(fpath))
            else
                call print_error("Unable to open the file (Code: " // trim(adjustl(code)) // ").")
            end if
            ierr = 1
            return
        end if
		
		! Open the file with write access
		! ierr = nf90_create(fpath, NF90_NETCDF4, iun)
			! if(ierr /= nf90_noerr) then
			  ! print *, trim(nf90_strerror(ierr))
			  ! stop "Stopped"
			! end if
		 
	end subroutine 
	
	!> Define  the dimensions
	subroutine	nc4_define_dimension_file(iun, dim_name, dim_length, did, ierr)
		
		!> Input variables.
        integer, intent(in) :: iun
        character(len = *), intent(in) :: dim_name

        !> Input variables (optional).
        integer, intent(in), optional :: dim_length

        !> Output variables.
        integer, intent(out) :: did, ierr
		
		!integer :: ierr

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Create dimension.
        if (present(dim_length)) then
            ierr = nf90_def_dim(iun, dim_name, dim_length, did)
        else
            ierr = nf90_def_dim(iun, dim_name, NF90_UNLIMITED, did)
        end if

        !> Check for errors.
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            call print_error( &
                "An error occurred adding the '" // trim(dim_name) // "' dimension to the file (Code " // &
                trim(adjustl(code)) // ").")
            ierr = 1
        else
            ierr = 0
        end if	
		! ierr = nf90_def_dim(iun, SLOPE_NAME, NSLOPES, SLOPE_dimid) ! NSLOPES
			! if(ierr /= nf90_noerr) then
			  ! print *, trim(nf90_strerror(ierr))
			  ! stop "Stopped"
			! end if	
				
	
	end subroutine 
	
	subroutine nc4_define_var(&
		iun, standard_name, dtype, & 
		dim1_id, & !dim2_id, dim3_id, & 
		!dim4_id, dim5_id, & 
		vid, ierr) 
	
		!> Input variables.
        integer, intent(in) :: iun, dtype
        character(len = *), intent(in) :: standard_name

        !> Output variables.
        integer, intent(out) :: ierr

        !> Output variables (optional).
        integer, intent(out), optional :: vid
		
		!> Input variables (optional).
        integer, intent(in), optional :: dim1_id!, dim2_id, dim3_id!, dim4_id, dim5_id
		
        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) dim_name, code
        integer did_c, d, v
	
	
		! if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id) .and. present(dim4_id) .and. present(dim5_id)) &
		    ! then
			! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id, dim3_id, dim4_id, dim5_id/), vid)
		! else if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id) .and. present(dim4_id)) then 
			! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id, dim3_id, dim4_id/), vid)
		! else if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id)) then 
			! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id, dim3_id/), vid)
		! if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id)) then 
			! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id, dim3_id/), vid)
		! else if (present(dim1_id) .and. present(dim2_id)) then 
		    ! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id/), vid)
		if (present(dim1_id)) then 
			ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id/), vid)
		else 
			ierr = nf90_def_var(iun, standard_name, dtype, (/0/), vid)
		end if 
		
		!ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, dim2_id/), vid)
		if (ierr /= NF90_NOERR) then
			call print_error("Unknown data type (" // trim(adjustl(code)) // ") for '" // trim(standard_name) // "' variable.")
            ierr = 1
		else 
			ierr = 0
		end if 
		 	
	end subroutine
	
	! subroutine nc4_define_var ( &
        ! iun, standard_name, dtype, &
        ! dim1_id, dim2_id, dim3_id, dim4_id, dim5_id, &
        ! vid, &
        ! ierr)

        ! !> Input variables.
        ! integer, intent(in) :: iun, dtype
        ! character(len = *), intent(in) :: standard_name

        ! !> Input variables (optional).
        ! !character(len = *), intent(in), optional :: name_dim_char_length
        ! integer, intent(in), optional :: dim1_id, dim2_id, dim3_id, dim4_id, dim5_id

        ! !> Output variables.
        ! integer, intent(out) :: ierr

        ! !> Output variables (optional).
        ! integer, intent(out), optional :: vid

        ! !> Local variables.
        ! character(len = DEFAULT_FIELD_LENGTH) dim_name, code
        ! integer did_c, d, v

        ! !> Create variable (based on type).
        ! select case (dtype)
            ! case (NF90_REAL, NF90_INT)

                ! !> Derive data type.
                ! if (dtype == NF90_REAL .and. kind(1.0) == 8) then
                    ! d = NF90_DOUBLE
                ! else
                    ! d = dtype
                ! end if

                ! !> Add the variable (for known numeric types).
                ! if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id) .and. present(dim4_id) .and. present(dim5_id)) &
                    ! then
                    ! ierr = nf90_def_var(iun, standard_name, d, (/dim1_id, dim2_id, dim3_id, dim4_id, dim5_id/), v)
                ! else if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id) .and. present(dim4_id)) then
                    ! ierr = nf90_def_var(iun, standard_name, d, (/dim1_id, dim2_id, dim3_id, dim4_id/), v)
                ! else if (present(dim1_id) .and. present(dim2_id) .and. present(dim3_id)) then
                    ! ierr = nf90_def_var(iun, standard_name, d, (/dim1_id, dim2_id, dim3_id/), v)
                ! else if (present(dim1_id) .and. present(dim2_id)) then
                    ! ierr = nf90_def_var(iun, standard_name, d, (/dim1_id, dim2_id/), v)
                ! else if (present(dim1_id)) then
                    ! ierr = nf90_def_var(iun, standard_name, d, (/dim1_id/), v)
                ! else
                    ! ierr = nf90_def_var(iun, standard_name, d, (/0/), v)
                ! end if
            ! ! case (NF90_CHAR)

                ! ! !> Prepare to add the variable for type 'char'.
                ! ! if (present(name_dim_char_length)) then

                    ! ! !> Copy the dimension name.
                    ! ! dim_name = trim(name_dim_char_length)
                ! ! else

                    ! ! !> Derive the dimension name.
                    ! ! write(code, FMT_GEN) DEFAULT_FIELD_LENGTH
                    ! ! dim_name = 'string' // trim(adjustl(code))

                    ! ! !> Check to see if the dimension already exists.
                    ! ! if (.not. nc4_inquire_dimension(iun, dim_name, did_c)) then
                        ! ! call nc4_define_dimension(iun, dim_name, DEFAULT_FIELD_LENGTH, did_c, ierr)
                        ! ! if (ierr /= 0) return
                    ! ! end if
                ! ! end if

                ! ! !> Add the variable (for type 'char').
                ! ! if (present(dim1_id)) then
                    ! ! ierr = nf90_def_var(iun, standard_name, dtype, (/dim1_id, did_c/), v)
                ! ! else
                    ! ! ierr = nf90_def_var(iun, standard_name, dtype, (/did_c/), v)
                ! ! end if
            ! case default

                ! !> Unknown type.
                ! write(code, FMT_GEN) dtype
                ! call print_error("Unknown data type (" // trim(adjustl(code)) // ") for '" // trim(standard_name) // "' variable.")
                ! ierr = 1
                ! return
        ! end select
        ! if (present(vid)) vid = v

        ! !> Check for errors.
        ! if (ierr /= NF90_NOERR) then
            ! write(code, FMT_GEN) ierr
            ! call print_error( &
                ! "An error occurred adding the '" // trim(standard_name) // "' variable (Code " // trim(adjustl(code)) // ").")
            ! ierr = 1
        ! else
            ! ierr = 0
        ! end if

    ! end subroutine
	

	subroutine nc4_add_data_1d_reall(iun, standard_name, vid, dat, ierr)

        !> Input variables.
        integer, intent(in) :: iun, vid
        character(len = *), intent(in) :: standard_name
        real, intent(in) :: dat(:)

        !> Input variables (optional).
        !integer, intent(in), optional :: start(:)

        !> Output variables.
        integer, intent(out) :: ierr

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Write data.
        ierr = nf90_put_var(iun, vid, dat)

        !> Check for errors.
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            call print_error( &
                "An error occurred saving the '" // trim(standard_name) // "' variable (Code " // trim(adjustl(code)) // ").")
            ierr = 1
            return
        else
            ierr = 0
        end if

    end subroutine
	
	subroutine nc4_add_data_2d_reall(iun, standard_name, vid, dat, start, ierr)

        !> Input variables.
        integer, intent(in) :: iun, vid
        character(len = *), intent(in) :: standard_name
        real, intent(in) :: dat(:,:)

        !> Input variables (optional).
        integer, intent(in), optional :: start(:)

        !> Output variables.
        integer, intent(out) :: ierr

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code

        !+ define variables 
		
		!+ set attributes 
		
		!> Write data.
        ierr = nf90_put_var(iun, vid, dat, start)

        !> Check for errors.
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            call print_error( &
                "An error occurred saving the '" // trim(standard_name) // "' variable (Code " // trim(adjustl(code)) // ").")
            ierr = 1
            return
        else
            ierr = 0
        end if

    end subroutine
	
	
	! added for define_var
	subroutine nc4_define_dimension(iun, dim_name, dim_length, did, ierr)

        !> Input variables.
        integer, intent(in) :: iun
        character(len = *), intent(in) :: dim_name

        !> Input variables (optional).
        integer, intent(in), optional :: dim_length

        !> Output variables.
        integer, intent(out) :: did, ierr

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Create dimension.
        if (present(dim_length)) then
            ierr = nf90_def_dim(iun, dim_name, dim_length, did)
        else
            ierr = nf90_def_dim(iun, dim_name, NF90_UNLIMITED, did)
        end if

        !> Check for errors.
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            call print_error( &
                "An error occurred adding the '" // trim(dim_name) // "' dimension to the file (Code " // &
                trim(adjustl(code)) // ").")
            ierr = 1
        else
            ierr = 0
        end if

    end subroutine
	
	! added for define_var
	logical function nc4_inquire_dimension(iun, dimension_name, did)

        !> Input variables.
        integer, intent(in) :: iun
        character(len = *), intent(in) :: dimension_name

        !> Output variables (optional).
        integer, intent(out), optional :: did

        !> Local variables.
        integer d

        !> Check if the dimension exists.
        nc4_inquire_dimension = (nf90_inq_dimid(iun, dimension_name, d) == NF90_NOERR)
        if (present(did)) did = d

    end function
	
	subroutine nc4_close_file(iun, fpath, quiet, ierr)

        !> Input variables.
        integer, intent(in) :: iun
        character(len = *), intent(in) :: fpath

        !> Input variables (optional).
        logical, intent(in), optional :: quiet

        !> Output variables.
        integer, intent(out) :: ierr

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        logical :: q = .true.

        !> Close the file.
        if (present(quiet)) q = quiet
        ierr = nf90_close(iun)
        if (ierr /= NF90_NOERR) then
            write(code, FMT_GEN) ierr
            if (q) then
                call print_error("An error occurred closing the file (Code: " // trim(adjustl(code)) // "): " // trim(fpath))
            else
                call print_error("Unable to close the file (Code: " // trim(adjustl(code)) // ").")
            end if
            ierr = 1
        else
            ierr = 0
        end if

    end subroutine
	
end module 