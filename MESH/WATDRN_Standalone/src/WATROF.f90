subroutine WATROF(THLIQ,THICE,ZPOND,TPOND,OVRFLW,TOVRFL,   &
                  SUBFLW,TSUBFL,RUNOFF,TRUNOF,FI,ZPLIM,    &
                  XSLOPE,XDRAINH,MANNING_N,DD,KSAT,TBARW,  & 
                  DELZW,THPOR,THLMIN,BI,DODRN,DOVER,DIDRN, &
                  ISAND,IWF,IG,ILG,IL1,IL2,BULK_FC)

!     * AUG 26/20 - D.PRINCZ.   CHANGED THE DIMENSIONS OF SUBFLW/TSUBFL
!                               TO PRESERVE THE PER-LAYER VALUES FOR
!                               INTERFLOW. THE TILE TOTALS ARE THE SUMS
!                               OF THE VALUES ALONG THE 2ND DIMENSION.
!     * MAR 03/10 - M.A.MEKONNEN/B.DAVISON/M.MACDONALD
!     *             RE-WRITTEN FOR TWO REASONS:
!     *             -TO USE VINCENT'S VERSION OF WATDRN;
!     *             -TO INCLUDE MORE COMMENTS.
!     * SEP 16/06 - R.SOULIS/F.SEGLENIEKS/A.PIETRONIRO/B.DAVISON.
!     *             MODIFICATIONS TO OVERLAND FLOW.
!     * SEP 15/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3.
!     * MAR 30/05 - D.VERSEGHY. ADDITIONAL FIELDS.
!     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
!     * AUG 02/02 - R.SOULIS/D.VERSEGHY. UPDATES DEVELOPED AT
!     *             WATERLOO.
!     * DEC 10/01 - R.SOULIS/K.SNELGROVE/T.WHIDDEN/D.VERSEGHY
!     *             WATFLOOD ROUTINE TO CALCULATE OVERLAND FLOW AND
!     *             INTERFLOW COMPONENTS OF SURFACE RUNOFF.
!
!
! ------------------ADDITIONAL CLARIFICATION-----------------------------------------------
! 
!      FOR A SLOPING ELEMENT, WATROF COMPUTES:
!          1. OVERLAND FLOW AT THE SURFACE BASED ON THE AVAILABLE PONDING DEPTH;
!          2. LATERAL FLOW FROM EACH LAYER BASED ON THE AVAILABLE WATER IN EACH SOIL LAYER.
! 
!      THE BASE FLOW FROM EACH LAYER IS ALSO COMPUTED BUT IT IS NOT USED AT PRESENT. SO,
!      THE TOTAL BASEFLOW REMAINS THE SAME AS THE ONE CALCULATED IN CLASS.
! 
! -----------------------------------------------------------------------------------------
!      DEFINITIONS
!      IWF         - FLAG GOVERNING OVERLAND AND LATERAL FLOW CALCULATIONS
!                    0 REPRESENTS FLAT ELEMENT - CLASS CALCULATES OVERLAND AND LATERAL FLOW
!                    NE 0 REPRESENTS SLOPING ELEMENT - WATROF CALCULATES OVERLAND AND LATERAL FLOW
!      IG          - TOTAL NUMBER OF SOIL LAYERS
!      ILG         - TOTAL NUMBER OF ELEMENTS
!      IL1         - STARTING INDEX OF ACTIVE ELEMENT
!      IL2         - FINAL INDEX OF ACTIVE ELEMENT
!      THLIQ       - VOLUMETRIC LIQUID WATER CONTENT OF SOIL LAYERS
!      THICE       - VOLUMETRIC FROZEN WATER CONTENT OF SOIL LAYERS
!      FI          - FRACTIONAL COVERAGE OF SUBAREA IN QUESTION ON MODELLED AREA
!      ZPLIM       - SUBAREA MAXIMUM PONDING DEPTH
!      XSLOPE      - SURFACE SLOPE
!      GRKFAC      - WATROF PARAMETER USED WHEN RUNNING MESH CODE ? NEEDS MORE CLARIFICATION
!      WFCINT      - WATROF PARAMETER USED WHEN RUNNING MESH CODE ? NEEDS MORE CLARIFICATION
!      TBARW       - TEMPERATURE OF WATER IN SOIL LAYER
!      ZPOND       - DEPTH OF PONDED WATER ON SURFACE
!      TPOND       - SUBAREA TEMPERATURE OF SURFACE PONDED WATER
!      OVRFLW      - OVERLAND FLOW FROM TOP OF SOIL COLUMN
!      TOVRFL      - TEMPERATURE OF OVERLAND FLOW
!      SUBFLW      - INTERFLOW FROM SIDES OF SOIL COLUMN
!      TSUBFL      - TEMPERATURE OF INTERFLOW FROM SIDES OF SOIL COLUMN
!      RUNOFF      - TOTAL RUNOFF
!      TRUNOF      - TEMPERATURE OF TOTAL RUNOFF
!      DELZW       - PERMEABLE THICKNESS OF SOIL LAYER
!      THPOR       - PORE VOLUME IN SOIL LAYER
!      THLMIN      - RESIDUAL SOIL LIQUID WATER CONTENT REMAINING AFTER FREEZING OR EVAPORATION
!      PSISAT      - SOIL MOISTURE SUCTION AT SATURATION
!      BI          - CLAPP AND HORNBERGER EMPIRICAL “B” PARAMETER
!      ISAND       - SAND CONTENT FLAG ? NEEDS MORE CLARIFICATION FOR THE VALUES
!      BULK_FC     - BULK FIELD CAPACITY
!      DELT        - TIME STEP
!      TFREZ       - FREEZING POINT OF WATER
!      MANNING_N   - MANNING'S ROUGHNESS COEFFICIENT
!      DD          - DRAINAGE DENSITY
!      ASAT_T0     - BULK SATURATION AT INITIAL TIME
!      ASAT_T1     - BULK SATURATION AT FINAL TIME
! 
!--------------------------------------------------------------
!+ todo 1: work on explanation of varibales and reformating them 
!       2: data precision to be real64 
!---------------------------------------------------------------
  implicit none
! 
! 
! * INPUT/OUTPUT ARRAYS.
! 
  real :: THLIQ (ILG,IG),  THICE (ILG,IG)

  real :: ZPOND (ILG),     TPOND (ILG),     OVRFLW(ILG), &
          TOVRFL(ILG),     SUBFLW(ILG,IG),               &
          RUNOFF(ILG),     TRUNOF(ILG)
! 
! * INPUT ARRAYS.
! 
  real :: FI    (ILG),    ZPLIM (ILG),     XSLOPE(ILG),   &
          xdrainh(ILG),    ksat(ILG),      TBARW (ILG,IG)
! 
! * SOIL INFORMATION ARRAYS.
! 
  real :: DELZW (ILG,IG), THPOR (ILG,IG),  THLMIN(ILG,IG), &
		  BI    (ILG,IG)

  integer :: ISAND (ILG,IG)
! 
! * WORK ARRAYS.
! 
  real :: DODRN (ILG),    DOVER (ILG), &
          DIDRN (ILG,IG), BULK_FC(ILG,IG)
! 
!  * COMMON BLOCK PARAMETERS.
! 
  
  ! NB : These two variables are called from CLASS1. Here it should set
  ! modified here as in the standalone these two varaibles are not called 
  ! from CLASS
  ! real :: DELT,TFREZ
!
  !COMMON /CLASS1/ DELT,TFREZ
  real, parameter :: TFREZ = 273.16			!< Freezing point of water (K)
  real, parameter :: DELT  = 1800.			!< Time step (s)
  

!  * INTERNAL SCALARS AND VECTORS
  real:: VEL_T0(ILG),NUC_DOVER(ILG),MANNING_N(ILG),DD(ILG), 		&
		 GRKEFF(ILG),ASAT_T0(ILG),ASAT_T1(ILG),DELZWJ(ILG),  		&
         BIJ(ILG),THPORJ(ILG),ASAT0(ILG),ASAT1(ILG),SATFC(ILG),     &
         DAVAIL,DTOT,SUBFLWJ(ILG),TSUBFL(ILG,IG),THLIQ_AVAIL(ILG),  &
         THPOR_AVAIL(ilg),BASFLWJ(ILG),XLAMBDA,ktop,kl,h0,c1,c2,    &
         ztop(ilg,ig)

  integer :: IG,ILG,IL1,IL2,I,J
  integer :: IWF(ILG)

  real :: exav

! -----------------------------------------------------------------------------------------
! return if no nml is expected to run in this cycle
  if(.not. any(iwf == 1)) return

! -----------------------------------------------------------------------------------------
!  coefficients
  c1 = 2.0/3.0
  c2 = 1.5 !3.0/2.0

! -----------------------------------------------------------------------------------------
! parameter - will be used to compute xdrainh (the fractional change in horizontal
! conductivity in a depth change h0) in Vincent's new formula.
  h0 = 1.0 
! -----------------------------------------------------------------------------------------
! loop through each element
	  	
  do i = il1,il2
! -----------------------------------------------------------------------------------------
!    skip if using flat class
     if (iwf(i) /= 1) cycle
!    ---------------------------------------------------------------------------------
!    compute overland flow and add to runoff and to the overall overland flow
!    ---------------------------------------------------------------------------------
     if(fi(i) .gt. 0.0 .and. zpond(i) .gt. zplim(i))then
!       ------------------------------------------------------------------------------
!       calculate the depth of water available for overland flow
!       ------------------------------------------------------------------------------
		dover(i) = zpond(i)-zplim(i)
!       ------------------------------------------------------------------------------
!       calculate the flow velocity at the beginning of the timestep
!       (based on kinematic wave velocity) - eqn (1) in notes on overland flow
!       ------------------------------------------------------------------------------
		vel_t0(i) = (dover(i)**c1)*sqrt(xslope(i))/(manning_n(i))
!       ------------------------------------------------------------------------------
!       calculate a normalized unconstrained overland flow to avoid numerical
!       problems with a division of small dover(i) values.
!       eqn (29) in notes on overland flow
!       ------------------------------------------------------------------------------
		nuc_dover(i) = -2*dd(i)*vel_t0(i)*delt
!       ------------------------------------------------------------------------------
!       constrained overland flow - limited by physically possible flow.
!       eqn (30) in notes on overland flow
!       ------------------------------------------------------------------------------
		dodrn(i) = dover(i)*(1.0-1./((1.0-c1*nuc_dover(i))**c2))
!       ------------------------------------------------------------------------------
!       add overland flow to runoff and to the overall overland flow
!       ------------------------------------------------------------------------------
		if(runoff(i) .gt. 1.0e-08) then
		   trunof(i) = (trunof(i)*runoff(i)+(tpond(i)+tfrez)*  &
						dodrn(i))/(runoff(i)+dodrn(i))
		endif
		runoff(i)    = runoff(i) + dodrn(i)
		if(dodrn(i) .gt. 1.0e-08)then
		   tovrfl(i) = (tovrfl(i)*ovrflw(i)+(tpond(i)+tfrez)*  &
						fi(i)*dodrn(i))/(ovrflw(i)+fi(i)*dodrn(i))
		   ovrflw(i) = ovrflw(i) + fi(i)*dodrn(i)
!       ---------------------------------------------------------------------------
!       subtract overland flow depth from the ponding depth
!       ---------------------------------------------------------------------------
		   zpond(i)  = zpond(i)- dodrn(i)
		endif
     endif
  enddo
! -----------------------------------------------------------------------------------------
!      compute interflow flow from each layer
! -----------------------------------------------------------------------------------------
  thliq_avail = 0.0
  thpor_avail = 0.0
  asat_t0     = 0.0
  ztop        = 0.0

! -----------------------------------------------------------------------------------------
! loop through each soil layer
  do j = 1,ig
!    ---------------------------------------------------------------------------------
!    loop through each element
	 do i = il1,il2
!    ---------------------------------------------------------------------------------
!    skip if not using watrof
	   if (iwf(i) /= 1) cycle
!      ---------------------------------------------------------------------------------
!      form vecotors for the layer - to be compatible with WATDRN arguments
	   delzwj(i)   = delzw(i,j)
	   bij(i)      = bi(i,j)
	   thporj(i)   = thpor(i,j)
!      ---------------------------------------------------------------------------------
!      Find the top of each soil layer for the calculation of grkeff
	   if(j .lt. ig)ztop(i,j+1) = ztop(i,j) - delzw(i,j)
	   if(fi(i) .gt. 0.0 .and. isand(i,j) .ge. -2 .and. &
		  delzw(i,j) .gt. 0.0)then
!          ---------------------------------------------------------------------------
!          determine available liquidwater in layer
!          ---------------------------------------------------------------------------
		   thliq_avail(i) = max(0.0,thliq(i,j)-thlmin(i,j))
!          ---------------------------------------------------------------------------
!          determine available porosity
!          ---------------------------------------------------------------------------
		   thpor_avail(i)    = max(thliq(i,j),thlmin(i,j), &
							   thpor(i,j)-thice(i,j))
!          ---------------------------------------------------------------------------
!          saturation defined as liquid water content over available porosity
!          ---------------------------------------------------------------------------
		   asat_t0(i)     = thliq_avail(i)/thpor_avail(i)
		endif
!           ------------------------------------------------------------------------------
!           grkeff - average value of the parameter controlling the time scale of
!           interflow process - kl * (tile slope / tile length) (1/s)
!           Note: this formula is not the same as the one in Fhydro2_VF_20100226.f
!           and needs to be confirmed by Vincent.
!           ------------------------------------------------------------------------------
!           Integration of k across the layer -> kl
		xlambda   = -log(xdrainh(i))/h0
		ktop      = ksat(i)*exp(xlambda*ztop(i,j))
		kl        = ktop * exav(xlambda*delzw(i,j))
		grkeff(i) = kl*xslope(i)*2.0*dd(i)/(1+xslope(i)**2)
		thpor_avail(i) = max(thlmin(i,j),thpor_avail(i))
	 enddo
!    ---------------------------------------------------------------------------------
!    compute interflow from the layer (subflowj). Baseflow from the layer (basflwj) is
!    also computed but is not used at present.
!    ---------------------------------------------------------------------------------
	 call watdrn (delzwj,bij,thpor_avail,ksat,grkeff,asat_t0,iwf, &
				  asat_t1,subflwj,basflwj,satfc, &  
				  ilg,il1,il2,1,delt)
!    ---------------------------------------------------------------------------------
!    loop through each element
	 do i = il1,il2
!       ------------------------------------------------------------------------------
!       skip if not using watrof
		if (iwf(i) /= 1) cycle
!       -----------------------------------------------------------------------------
!       allow lateral flow if liquid water content is greater than
!       bulk field capacity.
!       -----------------------------------------------------------------------------
		if(thliq_avail(i).gt.0.0.and.thliq(i,j).ge.bulk_fc(i,j))then
		   didrn(i,j) = subflwj(i)
!          ---------------------------------------------------------------------------
!          compute davail: volume of available water in a soil layer per land
!                              area [m^3/m^2]
!          ---------------------------------------------------------------------------
		   davail = thliq_avail(i)*delzw(i,j)
!          ---------------------------------------------------------------------------
!          limit the lateral flow not to exceed the available water in the layer
!          ---------------------------------------------------------------------------
		   didrn(i,j) = max(0.0,min(davail,didrn(i,j)))
!          ---------------------------------------------------------------------------
!          add the lateral flow to the runoff and to the subflow
!          ---------------------------------------------------------------------------
		   if(didrn(i,j).gt.1.0e-8)then
			  trunof(i)  = (trunof(i)*runoff(i)+tbarw(i,j)* &
						   didrn(i,j))/(runoff(i)+didrn(i,j))
			  runoff(i)  = runoff(i)+didrn(i,j)
			  subflw(i,j)= subflw(i,j)+fi(i)*didrn(i,j)
!             ------------------------------------------------------------------------
!             remove the lateral flow from the layer
!             ------------------------------------------------------------------------
			  thliq(i,j) = thliq(i,j)-didrn(i,j)/delzw(i,j)
		   endif
		endif
	 enddo
  enddo
  
  return
end subroutine WATROF

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
