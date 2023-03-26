!> AUTHOR : Ala Bahrami
!> DATE CREATION : 03-23-2023
!> DATES MODIFICATIONS : -
!> DESCRIPTION : The purpose of this script is to run the stand-alone version of WATDRN 
!>               and do some experiments based on paper of WATDRN critique. 
!> see also : WATDRN_SA_simulation.py 
program WATDRN_SA_simulation

	use netcdf
	use nc_io_file
    use check_stat
	implicit none	 

	INTERFACE 
	  FUNCTION exav(x)
		! REAL, exav
		! REAL, INTENT(IN) :: x
	  END FUNCTION exav
	END INTERFACE

	! constant parameters for synthetic simulation 
	integer, parameter :: MM_PER_M     = 1000  ! [-]
	integer, parameter :: S_PER_HOUR   = 3600  ! [-]
	
	integer, parameter :: TOPMODEL_exp = 3.0   ! [-]
	real   , parameter :: surfHydCond  = 1.0   ! [m hour^-1]
	 
	real   , parameter :: tan_slope    = 0.3   ! [-]
	real   , parameter :: soilDepth    = 1.5   ! [m] 
	real   , parameter :: porosity     = 0.25  ! [-]
	real   , parameter :: qRain        = 0.002 ! [m hour^-1] 

	integer, parameter :: totalLength1 = 50    ! [m]
	integer, parameter :: totalLength2 = 10    ! [m]

	integer, parameter :: hillWidth    = 100   ! [m]
	
	integer, parameter :: ILG = 1
	integer, parameter :: IG  = 1
	integer, parameter :: IL1 = 1
	integer, parameter :: IL2 = 1
	integer, parameter :: IWF = 1
	
	! save interflow output 
    integer :: iun = 20	
    
	character(len=*), parameter :: fpath = 'interflow.nc'
	integer :: time_varid = 21, interflow_varid = 22 
	
	! time duration of simulation and time step  
	integer, parameter :: duration   = 20*24           ![hours]
	real, parameter    :: delt       = 1          	   ![hours]  
        
	! input forcing 
	real:: precip (duration) 
	
	! work array and internal variables 
	real :: kH, Dd1, Dd2, bij, xlambda, delzw, ksat, xslope, &
		    h0, ktop, kl, grkeff1, grkeff2, thpor_avail, ztop, &
			sWATDRN1(duration), qWATDRN1(duration), &
			sWATDRN2(duration), qWATDRN2(duration)
	
	integer :: i, ierr
	
	! output variables 
	real :: asat_t1(duration), subflw1(duration), & 
			basflw1(duration), satfc1(duration), &
			asat_t2(duration), subflw2(duration), & 
			basflw2(duration), satfc2(duration)
    
	!-------------------------------------------------------------	
	! construct variables and initialize them  
	! construct precipitation
	precip(1:duration/2) = qRain                      ![m hour^-1]
	precip(duration/2+1:duration) = 0.0

	! Horizontal conductivity 
	kH = surfHydCond * tan_slope              		  ![m hour^-1]

    !-------------------------------------------------------------
	! call WATDRN_storage subroutine 
	call WATDRN_storage(precip, kH, totalLength1, soilDepth, &
						porosity, TOPMODEL_exp, delt,duration , &
						sWATDRN1, qWATDRN1)
	
	call WATDRN_storage(precip, kH, totalLength2, soilDepth, &
						porosity, TOPMODEL_exp, delt,duration , &
						sWATDRN2, qWATDRN2)					
    
	!-------------------------------------------------------------
	! initialize variables and parameters requried for running WATDRN 
	Dd1 = 1.0/(2.0*totalLength1)                           ! [1/m]
	Dd2 = 1.0/(2.0*totalLength2)                           ! [1/m]

	! obtain b coefficient  based on Clapp-Hornberger connectivity index  
	bij = (TOPMODEL_exp -3)/2                          ! [-]

	! calculating xlambda explicitly by equating eq.(21) from Martyn's paper and eq.(24) from Mekonnen's draft 
	! NB: in MESH implementation, the xlambda decay factor is derived from xdrainh
	xlambda = - log(1+tan_slope**2)/(soilDepth)        ![1/m]
	xlambda = - xlambda                                ! to be consistent with fortran code 		


	delzw  = soilDepth                                 ! [m]
	ksat   = surfHydCond                               ! [m hour^-1]
	xslope = tan_slope                                 ! [-]

	! parameter - will be used to compute xdrainh (the fractional change in horizontal
	! conductivity in a depth change h0) in Vincent's new formula.
	h0 = 1.0 

	! Integration of saturated conductivity across the layer -> kl
	! Important note :xlambda is derived based on function of slope and soil Depth instead of xdrainh 
	! xlambda        = -np.log(xdrainh[i])/h0 #  
	ktop           = ksat * exp(xlambda*ztop)          !  [m hour^-1]

	! Important note : I commented calculation of kl using exav functions as it 
	! reduces the precision of calculation of grkeff
	!kl             = ktop * exav(xlambda*delzw)       !  [m hour^-1]
	kl             = ktop * exp(xlambda*delzw)         !  [m hour^-1]

	! calculate grkeff for two hillslope lengths 
	grkeff1         = kl*xslope*2.0*Dd1/(1+xslope**2)  ! [hour^-1]
	grkeff2         = kl*xslope*2.0*Dd2/(1+xslope**2)  ! [hour^-1]

	! thpor_avail[i] = np.max((thlmin[i,j],thpor_avail[i]))
	thpor_avail = porosity 
	
	!-------------------------------------------------------------
	! call WATDRN for entire time window of simulation 
	asat_t1 = 0
	subflw1 = 0
    basflw1 = 0
	satfc1  = 0 
	asat_t2 = 0
	subflw2 = 0 
	basflw2 = 0
	satfc2  = 0

    !  compute interflow from the layer (subflow). Baseflow from the layer (basflw) is
	!  also computed but is not used at present.
	
	
	
	do i = 1, duration 	
	   
	 call watdrn (delzw,bij,thpor_avail, ksat, grkeff1, sWATDRN1(i),iwf, &
				  asat_t1(i),subflw1(i),basflw1(i),satfc1(i), &  
				  ilg,il1,il2,1,delt)
				  
	 ! convert interflow components to  [mm hour^-1]
	 subflw1(i) = MM_PER_M * subflw1(i)/delt                
	 basflw1(i) = MM_PER_M * basflw1(i)/delt
	 
     
	 call watdrn (delzw,bij,thpor_avail, ksat, grkeff2, sWATDRN2(i),iwf, &
				  asat_t2(i),subflw2(i),basflw2(i),satfc2(i), &  
				  ilg,il1,il2,1,delt)
				  
	 ! convert interflow components to  [mm hour^-1]
	 subflw2(i) = MM_PER_M * subflw2(i)/delt                
	 basflw2(i) = MM_PER_M * basflw2(i)/delt     
	
	 
	end do 
	
	!-------------------------------------------------------------
	! save output 
	! intialize saving the output 
	call nc4_init_file(fpath, duration, 8, &
					  iun, &
					  interflow_varid, time_varid)
					  
	! Write and append data each time step  
	 call nc4_add_data_file(iun, duration, 8,                     &
							 interflow_varid, time_varid,         & 
							 (/asat_t1,subflw1,basflw1,satfc1,    &
							  asat_t2,subflw2,basflw2,satfc2/) , 1)				  
	
	! close output file 
	call nc4_close_fille(iun)
	
	
end program 

! **********************************************************************
! 
!      function exav(x)
!      finds average of an exponential function exp(-x)
!      over the interval [0,x], which is (1-exp(-x))/x
!      deals with limit values of 1-x/2 when x->0 and 1/x when x->inf
! 
! **********************************************************************
!
function exav(x)
	implicit none
	!
	real, intent(in) :: x
	real :: exphuge, expbig, expsmall, exav
	data exphuge/1.0e+9/,expbig/1.0e+8/,expsmall/1.0e-8/
    !
	if (x .gt. exphuge) then
	   exav = 0.0
	elseif (x .gt. expbig) then
	   exav = 1.0/x
	elseif (x .gt. expsmall) then
	   exav = (1.0-exp(-x))/x
	else
	   exav = 1.0-x/2.0
	endif
	return
end	