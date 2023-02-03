SUBROUTINE WATDRN(delzw,bcoef,thpora,grksat,grkeff,asat0,iwf, &
     			  asat1,subflw,basflw,satsfc, &
                  ilg,il1,il2,wfk,delt)
!
!     * December 4, 2009, Vincent Fortin
!
!       Initial WATDRN code derived from WAT_DRAIN by Ric Soulis
!       Code entirely rewritten for three main reasons:
!       - simplify the parameterization by working with average
!         hydraulic conductivity within a layer (runs faster, but
!         more importantly is easier to understand for users)
!       - meet CCC and RPN source code standards
!       - improve readability
!
!     * December 11, 2009, Vincent Fortin
!
!       - Define grkeff as ksat * (tile slope / tile length) (1/s)
!         (take out thpor, as thpor changes with ice content)
!       - Modify computation of critical time accordingly:
!         tc = thpor / (c*grkeff)
!
!     * March 22, 2010, Vincent Fortin
!
!       - Stop computing t0 and t1 to avoid overflow, compute tc/t
!         directly when t>tc and t/tc when t<tc (where t=t0 or t1)
!       - Avoid dividing by tc to make the code more robust
!         when tc is small
!       - Simplify calculation of baseflow rate
!       - Bugfix: move satsfc calculation after computation of tc/t0
!       - Change variable name from thpor to thpora
!
!     Summary:
!
!     This routine calculates the outflow from a tilted landscape
!     element (tile). Rather than using average saturation for the 
!     element(asat0), it finds flow using the saturation at the
!     seepage face (for interflow) at at the bottom of the layer
!     (for baseflow). The routine returns:
!     1) the saturation at the end of the time step (asat1)
!     2) the interflow (subflw) and baseflow (basflw) integrated 
!        over the duration of the time step (in m)
!     3) the fraction of the surface which is saturated at the
!        surface according to the model (satsfc)
!
!     Interflow amount during the time step in m is determined from:
!     subflw = (asat0 - asat1) * thpor * delzw
!     (implicit estimation)
!
!     Baseflow amount during the time step in m is determined from:
!     basflw = grksat * asatb0**(2*bcoef+3) * delt
!     where asatb0 is the average saturation across the bottom of the
!     element at initial time (explicit estimation, 1st order)
!
!     The fraction of the surface that is saturated at the initial
!     time (satsfc) is estimated as MIN(1-t0/tc,0) where t0 is the 
!     (theoretical) time elapsed since saturation and tc is the
!     critical time at which the saturation front reaches the
!     seepage face.
!
!     Background info on WATDRN:
!
!     WATDRN is a parameterization for sub-grid scale interflow,
!     i.e. subsurface flow, which is thought to be an important flow
!     generation process following a rain event at the basin scale
!
!     The underlying principles behind WATDRN are described in
!     a paper by Soulis et al. (2000), Atmosphere-Ocean.
!     However, this code presents a simpler solution to the problem
!     which assumes that hydraulic conductivity is constant within
!     a soil layer. This is a departure from previous versions of
!     WATDRN aimed at making the code easier to understand
!     and faster to run.
!
!     Here is the basic idea: interflow is especially important
!     after a rain event which leaves the soil close to saturation.
!     Following such a rain event, water flows from the hillslope
!     to the river which creates a gradient in water content along
!     the hillslope. But land-surface models have generally as their
!     prognostic variable the mean water content of the grid box.
!     Under the hypothesis of Darcian flow along the hillslope
!     following a rain event which left the soil saturated, WATDRN
!     recovers the saturation distribution along the hillslope
!     from the bulk saturation, estimates from this the interflow
!     rate at the seepage face and integrates this rate over the
!     duration of the time step.
!
!     Given a bulk saturation value asat_t0 at the start of the
!     time step, and assuming that this bulk saturation is the
!     result of Darcian flow in the direction of the hillslope
!     starting from a saturated soil (with no rain after that),
!     there is a one-to-one relationship between bulk saturation
!     and the time elapsed since the rain stopped. So we can figure
!     out:
!
!     1. t0, the (theoretical) time elapsed since the soil was
!        saturated, knowing the bulk saturation asat0.
!     2. t1, the (theoretical) time at the end of the time step
!        (that's just t0+delt)
!     3. asat1, the bulk saturation at the end of the time step
!
!     Then interflow is proportional to the difference between
!     asat0 and asat1.
!
!     From bulk saturation, it is also possible to figure out
!     average saturation across the bottom of the element, which
!     is used to estimate baseflow, and the fraction of the surface
!     that is saturated, which can be used to separate runoff and
!     infiltration.

  implicit none

! Input parameters
  integer ilg         ! Size of arrays
  integer il1         ! index of first grid point to process
  integer il2         ! index of last grid point to process
  integer wfk         ! IWF of parent call (e.g., WATROF, LATFLOW)
  real    delt        ! duration of the time step (s)

! Input arrays
  real    delzw(ilg)  ! layer thickness (m)
  real    bcoef(ilg)  ! slope of retention curve, related to
					  ! Clapp-Hornberger connextivity index by c=2b+3
  real    thpora(ilg) ! Available porosity of the soil layer
					  ! (total porosity - ice content)
  real    grksat(ilg) ! Vertical hydraulic conductivity at saturation
					  ! at the bottom of the layer (m/s)
  real    grkeff(ilg) ! average value of the parameter
					  ! controlling the time scale of interflow process
					  ! ksat * (tile slope / tile length) (1/s)
  real    asat0(ilg)  ! bulk saturation at initial time
  integer iwf(ilg)    ! IWF flag (runs if == wfk)

! Output arrays
  real    asat1(ilg)  ! bulk saturation at the end of the time step
  real    subflw(ilg) ! interflow amount during the time step (m)
  real    basflw(ilg) ! baseflow rate during the time step (m)
  real    satsfc(ilg) ! saturated fraction of the surface (0 to 1)

! Work arrays
  real    c(ilg)      ! Clapp-Hornberger connectivity index (c>1)
  real    cm1(ilg)    ! c-1
  real    c2m1(ilg)   ! 2*c-1
  real    asatc(ilg)  ! bulk saturation at the critical time tc
  real    asat00(ilg) ! MIN(1,asat0) because we don't handle supersaturated soils
  real    tc(ilg)     ! critical time at which the seepage face becomes unsaturated
  real    ratiot(ilg) ! ratio tc/t (t0 or t1) if t>tc and t/tc if t<=tc
  logical satspf(ilg) ! indicates if seepage face is saturated
					  ! equivalent to knowing if t<=tc

! Local variables
  integer i           ! Array index

! return if no nml is expected to run in this cycle
  if (.NOT. any(iwf == wfk)) return
!  
! **********************************************************************
!      STEP 0: Initialize a few things before we start
!              - output variables
!              - functions of the input parameters
! **********************************************************************
  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk .or. asat0(i) == 0.0) cycle
!    c and c factors
	 c(i)    = 2.*bcoef(i)+3.
	 cm1(i)  = c(i)-1.
	 c2m1(i) = 2.*c(i)-1.
!    bulk saturation at critical time
!    (just before the seepage face becomes unsaturated)
	 asatc(i) = 1.-1./c(i)
!    layer average saturation asat0 may be greater than 1
!    e.g. frost heave but it is not possible for wat_drain
	 asat00(i) = MIN(1.,asat0(i))
!    assess if seepage face is saturated at initial time
	 satspf(i) = asat00(i) .GE. asatc(i)
  enddo
! **********************************************************************
!      STEP 1: Find theoretical time t0 elapsed since element was last
!              saturated and estimate baseflow rate at initial time
!              Also estimate fraction of surface that is saturated
! ********************************************************************** 
  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
!    determine time at which seepage face becomes unsaturated
	 tc(i) = thpora(i)/(c(i)*grkeff(i))
  enddo

  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
!    find theoretical start of recession (t0) from bulk saturation
!    and at the same time estimate baseflow based on rate at t0
	 if (satspf(i)) then 
!       saturated seepage face at initial time:
!       compute t0/tc
		ratiot(i) = c(i)*(1.-asat00(i))
!       normalized baseflow rate
		basflw(i) = 1.-c(i)*c(i)/c2m1(i)*(1.-asat00(i))
!       the fraction of the surface that is saturated at t0
!       varies linearly with t0/tc
		satsfc(i) = 1.-ratiot(i)
	 else
!       unsaturated seepage face at initial time:
!       calculate tc/t0 instead of t0 to avoid overflow
		ratiot(i) = (asat00(i)/asatc(i))**cm1(i)
!       normalized baseflow rate
		basflw(i) = cm1(i)/c2m1(i)*ratiot(i)*asat00(i)/asatc(i)
!       the fraction of the surface that is saturated at t0 is zero
		satsfc(i) = 0.
	 endif
  enddo

  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
!    Compute baseflow in m from normalized baseflow rate
	 basflw(i) = grksat(i)*basflw(i)*delt
  enddo
! **********************************************************************
!      STEP 2: Find theoretical time t1 at the end of the time step
! **********************************************************************
  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
	 if (satspf(i)) then
!       Assess if seepage face will still be saturated at the
!       end of the time step
		satspf(i) = tc(i)*ratiot(i)+delt .LE. tc(i)
		if (satspf(i)) then
!          Seepage face still saturated, compute t1/tc from t0/tc
		   ratiot(i) = (tc(i)*ratiot(i)+delt)/tc(i)
		else
!          Seepage face not saturated anymore, compute tc/t1
		   ratiot(i) = tc(i)/(tc(i)*ratiot(i)+delt)
		end if
	 else 
!       If seepage face was not saturated initially, we compute
!       tc/t1=tc/(t0+delt) from tc/t0
		ratiot(i) = tc(i)*ratiot(i)/(tc(i)+delt*ratiot(i))
	 end if
  enddo
! **********************************************************************
!      STEP 3: Obtain bulk saturation at the end of the time step
!              and interflow amount
! **********************************************************************
  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
	 if (satspf(i)) then
!       saturated seepage face at the end of the time step
		asat1(i) = 1.-ratiot(i)/c(i)
	 else
!       unsaturated seepage face at the end of the time step
		asat1(i) = asatc(i)*ratiot(i)**(1./cm1(i))
	 endif
  enddo

  do i=il1,il2
!    cycle if not using watrof or latflow
	 if (iwf(i) /= wfk  .or. asat0(i) == 0.0) cycle
!    Sanity check: bulk saturation should not increase with time
	 asat1(i) = MIN(asat00(i),asat1(i))
!    Obtain interflow from the difference in bulk saturation
	 subflw(i) = (asat00(i)-asat1(i))*thpora(i)*delzw(i)
  enddo

  return
end subroutine 
