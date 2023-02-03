!> AUTHOR : Ala Bahrami
!> DATE CREATION : 02-01-2023
!> DATES MODIFICATIONS : -
!> DESCRIPTION : The purpose of this program is to run WATDRN as stand-alone 
!> and by calling both WATROF and WATDRN subroutines. The input variables and                
!> and parameters are obtained from MESH Fraser setup which was run for one month
!> source      : 
!>------------------------------------------------------------------------------- 
!
! todo : 1) Add unity and variable description   
!		
!>------------------------------------------------------------------------------- 
program WATDRN_SA

	use netcdf
	use nc_io_file
    use check_stat
	implicit none	

	! input data names  
	character(len=*), parameter :: fpath_cs = 'THLQCS_2000_09_01.nc'
	character(len=*), parameter :: fpath_gs = 'THLQGS_2000_09_01.nc'
	character(len=*), parameter :: fpath_c  = 'THLQC_2000_09_01.nc'
	character(len=*), parameter :: fpath_g  = 'THLQG_2000_09_01.nc'
	character(len=*), parameter :: fpath_p  = 'parameter.nc'
	   
	! variables dimensiosn 
	character (len = *), parameter :: THLQ_NAME  ="soil_moisture"
	character (len = *), parameter :: Level_NAME ="level"
	character (len = *), parameter :: Time_NAME  ="time"
		
	! constant parameters
	integer, parameter :: ILG = 15236, IG = 4, IL1 = 1, IL2 = ILG
	integer, parameter :: levelp = 30, level = 34, t_steps = 48 !(2*24 for one day)
	
	! local variables 
	integer :: i	
    
	! Diagnostic array 
	real :: OVRFLW    (ILG), TOVRFL    (ILG), &
            SUBFLW (ILG,IG), TSUBFL (ILG,IG)
	
	! work array for WATROF
	real :: DODRN(ILG), DOVER(ILG), DIDRN (ILG,IG)
	
	! canopy over snow (CS)
	real :: THLQCS(ILG,IG), THICCS(ILG, IG), ZPNDCS   (ILG), &
			TPNDCS   (ILG), RUNFCS    (ILG), TRNFCS   (ILG), &
			FCS      (ILG), ZPLMCS    (ILG), TBRWCS(ILG,IG)
            
	! snow covered ground (GS)
	real :: THLQGS(ILG,IG), THICGS(ILG, IG), ZPNDGS   (ILG), &
			TPNDGS   (ILG), RUNFGS    (ILG), TRNFGS   (ILG), &
			FGS      (ILG), ZPLMGS    (ILG), TBRWGS(ILG,IG)
	 
	! canopy over bare ground (C)
	real :: THLQCO(ILG,IG), THICCO(ILG, IG), ZPONDC   (ILG), &
			TPONDC   (ILG), RUNFC     (ILG), TRUNFC   (ILG), &
			FC       (ILG), ZPLIMC    (ILG), TBARWC(ILG,IG)
	
	! bare ground (G)
	real :: THLQGO(ILG,IG), THICGO(ILG, IG), ZPONDG   (ILG), &
			TPONDG   (ILG), RUNFG     (ILG), TRUNFG   (ILG), &
			FG       (ILG), ZPLIMG    (ILG), TBARWG(ILG,IG)
	
	! time-invariant variables 
	real :: XSLOPE     (ILG), XDRAINH  (ILG), MANNING_N(ILG), &
			DD         (ILG), KSAT     (ILG), DELZW (ILG,IG), & 
			THPOR   (ILG,IG), THLMIN(ILG,IG), BI    (ILG,IG), &
			BULK_FC(ILG, IG)
	integer :: ISAND(ILG, IG), IWF(ILG)
	
	!> Read data 
	!---------------------------------------------------------------------------
	! store data into these local variables 
	real :: datp (ILG, levelp, 1)
	real :: datcs(ILG, level, t_steps-1)
	real :: datgs(ILG, level, t_steps-1)
	! NB: the Canopy and Ground subareas have one more time step
	real :: datc (ILG, level, t_steps)
	real :: datg (ILG, level, t_steps)
	
	! Read parameter file  
	call nc4_read_data_file(fpath_p , datp , THLQ_NAME)
	call nc4_read_data_file(fpath_cs, datcs, THLQ_NAME)
	call nc4_read_data_file(fpath_gs, datgs, THLQ_NAME)
	call nc4_read_data_file(fpath_c , datc , THLQ_NAME)
	call nc4_read_data_file(fpath_g , datg , THLQ_NAME)
	!> Assign variables before calling WATROF  
	!---------------------------------------------------------------------------
	! NB: The order of coloums corresponds the way data are written in CLASSW 
	! parameters
	XSLOPE    = datp(:,1,1) 	   
	XDRAINH   = datp(:,2,1) 
	MANNING_N = datp(:,3,1)       
	DD 	      = datp(:,4,1) 
	KSAT      = datp(:,5,1)
	DELZW 	  = datp(:,6:9,1)  
	THPOR     = datp(:,10:13,1)
	THLMIN    = datp(:,14:17,1)
	BI 	      = datp(:,18:21,1)
	ISAND     = int(datp(:,22:25,1))
	IWF       = int(datp(:,26,1))
	BULK_FC   = datp(:,27:30,1)  
	!--------------------------------------------------------------------------
	! loop over time steps and call WATROF for each subareas
	!--------------------------------------------------------------------------
	do  i = 1, t_steps
		!--------------------------------------------------------------------------- 
		!> C: Calculations for Canopy over bare ground  
		!--------------------------------------------------------------------------- 
		THLQCO = datc(:,1:4,i)    
		THICCO = datc(:,5:8,i) 
		ZPONDC = datc(:,9,i)     
		TPONDC = datc(:,10,i)
		OVRFLW = datc(:,11,i)
		TOVRFL = datc(:,12,i)
		SUBFLW = datc(:,13:16,i)
		TSUBFL = datc(:,17:20,i) 
		RUNFC  = datc(:,21,i)
		TRUNFC = datc(:,22,i)
		FC     = datc(:,23,i)
		ZPLIMC = datc(:,24,i)
		TBARWC = datc(:,25:28,i)
		
		call WATROF(THLQCO, THICCO, ZPONDC, TPONDC, OVRFLW, TOVRFL, &
					SUBFLW, TSUBFL, RUNFC, TRUNFC, FC, ZPLIMC,      &
					XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBARWC,   &
					DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN,  &
					ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
					
		write (*,*) 'SUBFLWC(90,:):', SUBFLW(90,:)
		!--------------------------------------------------------------------------- 
		!> G: Calculations for bare ground 
		!--------------------------------------------------------------------------- 
		THLQGO = datg(:,1:4,i)    
		THICGO = datg(:,5:8,i) 
		ZPONDG = datg(:,9,i)     
		TPONDG = datg(:,10,i)
		OVRFLW = datg(:,11,i)
		TOVRFL = datg(:,12,i)
		SUBFLW = datg(:,13:16,i)
		TSUBFL = datg(:,17:20,i) 
		RUNFG  = datg(:,21,i)
		TRUNFG = datg(:,22,i)
		FG     = datg(:,23,i)
		ZPLIMG = datg(:,24,i)
		TBARWG = datg(:,25:28,i)
		
		call WATROF(THLQGO, THICGO, ZPONDG, TPONDG, OVRFLW, TOVRFL,	&
					SUBFLW, TSUBFL, RUNFG, TRUNFG, FG, ZPLIMG,      &
					XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBARWG,   & 
					DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN,  &
					ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
		
		write (*,*) 'SUBFLWG(90,:):', SUBFLW(90,:)
		if (i < t_steps-1) then 
			!--------------------------------------------------------------------------- 
			!> CS: Calculation for canopy over snow 
			!--------------------------------------------------------------------------- 
			THLQCS = datcs(:,1:4,i)    
			THICCS = datcs(:,5:8,i) 
			ZPNDCS = datcs(:,9,i)     
			TPNDCS = datcs(:,10,i)
			OVRFLW = datcs(:,11,i)
			TOVRFL = datcs(:,12,i)
			SUBFLW = datcs(:,13:16,i)
			TSUBFL = datcs(:,17:20,i) 
			RUNFCS = datcs(:,21,i)
			TRNFCS = datcs(:,22,i)
			FCS    = datcs(:,23,i)
			ZPLMCS = datcs(:,24,i)
			TBRWCS = datcs(:,25:28,i)
			
			call WATROF(THLQCS, THICCS, ZPNDCS, TPNDCS, OVRFLW, TOVRFL,  &
						 SUBFLW, TSUBFL, RUNFCS, TRNFCS, FCS, ZPLMCS,    &
						 XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBRWCS,   &
						 DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN,  &
						 ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
			
			write (*,*) 'SUBFLWCS(90,:):', SUBFLW(90,:)
			!--------------------------------------------------------------------------- 
			!> GS: Calculations for snow-covered ground 
			!---------------------------------------------------------------------------
			THLQGS = datgs(:,1:4,i)    
			THICGS = datgs(:,5:8,i) 
			ZPNDGS = datgs(:,9,i)     
			TPNDGS = datgs(:,10,i)
			OVRFLW = datgs(:,11,i)
			TOVRFL = datgs(:,12,i)
			SUBFLW = datgs(:,13:16,i)
			TSUBFL = datgs(:,17:20,i) 
			RUNFGS = datgs(:,21,i)
			TRNFGS = datgs(:,22,i)
			FGS    = datgs(:,23,i)
			ZPLMGS = datgs(:,24,i)
			TBRWGS = datgs(:,25:28,i)
			
			call WATROF(THLQGS, THICGS, ZPNDGS, TPNDGS, OVRFLW, TOVRFL,	 &
						SUBFLW, TSUBFL, RUNFGS, TRNFGS, FGS, ZPLMGS,	 & 
						XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBRWGS,	 &
						DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN,   &
						ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
			
			write (*,*) 'SUBFLWGS(90,:):', SUBFLW(90,:)
		end if 
	enddo
end program 