!> AUTHOR : Ala Bahrami
!> DATE CREATION : 01-11-2023
!> DATES MODIFICATIONS : -
!> DESCRIPTION : Call the check statment for the errors
!>               
!>
module check_stat

#ifdef NETCDF
    use netcdf
#endif
	implicit none
	
	contains
	
	subroutine check(status)

	   !use netcdf 
	   !implicit none 

	   integer, intent (in) :: status

	   if(status /= nf90_noerr) then 
			print *, trim(nf90_strerror(status))
			stop "Stopped"
	   end if
	end subroutine ! check 

end module 