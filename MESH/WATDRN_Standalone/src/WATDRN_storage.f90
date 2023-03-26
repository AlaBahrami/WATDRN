subroutine WATDRN_storage(qForce, kH, xL, zSoil, phi, eta, dt,tsteps , &
						  sWATDRN, qWATDRN)
	
  ! the purpose of this subroutine is to produce relativbe storage required 
  ! for running WATDRN for synthetic test
	
  ! The subroutine is the translation of the WATDRN_SA_simulation.py script
  ! which is based on native WATDRN algorithm 

  implicit none
   
  ! input parameters 
  integer eta
  integer xL 
  integer tsteps 
  
  real kH
  real zSoil
  real phi
  real dt 
  
  ! input array
  real qForce(tsteps)
  
  ! output arrays 
  real sWATDRN(tsteps)
  real qWATDRN(tsteps)
  
  ! work arrays 
  real :: S0, S1, Sc,Sn, &
		  tf,  tr, t1,  & 
		  qx 
  
  ! local variable 
  integer i 
  
  ! intialize variables 
  !--------------------------------------------------------
  !initialize (mx points distributed across the full space)
  S0 = 1.e-8
    
  !save variables
  sWATDRN = 0
  qWATDRN = 0

  ! get temporally constant variables
  ! time that the drying front reaches the bottom of the hillslope
  tf = phi*xL/(kH*eta) 
  ! spatially-averaged relatve storage at time tf
  Sc = 1.0 - 1.0/eta   
  
  ! loop through points on the characteristic curve
  do i = 1,tsteps 
	 
	 ! compute the time required for spatially-averaged storage to reach the state from the host land model
     ! - only consider case 2: the toe of the hillslope is not fully saturated
     tr = tf*(S0/Sc)**(1.0 - eta)
	 
	 ! compute the storage at the time tr + Delta t
     t1 = tr + dt
     S1 = Sc*(tf/t1)**(1.0/(eta - 1.0))
        
     ! compute the lateral flow from mass balance
     qx = phi*zSoil*(S0 - S1)/dt
     
        
     !update states
     Sn = S0 + (qForce(i) - qx)/(zSoil*phi)
     S0 = Sn
        
     ! save data
     sWATDRN(i)   = S0
     qWATDRN(i)   = qx
		
  end do 
  
end subroutine 