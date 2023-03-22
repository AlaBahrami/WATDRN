# -*- coding: utf-8 -*-
"""
Purpose: To produce graphds based on MoC for WATDRN's paper

Created on Thu Feb 16 23:39:11 2023

See also: https://github.com/CH-Earth/laughTests/tree/master/lt4_wigmosta1994

@author: Ala Bahrami

last modified: 03/14/2023
              1) include WATDRN stand-alone routines 
              2) add WATDRN simulations to MoC and Soulis' upscaling method 
                
              03/21/2023
              1) Add WATDRN storage required for running stand-alone WATDRN for 
                synthetic test based on section 2.3.2 and 2.4 paper of Martyn et al. (2023)
"""
#%% import modules 
import matplotlib
import matplotlib.pyplot as plt
#import xarray as xr
import numpy as np

#%% I/O directory 
outdir           = 'output/'

#%% input parameters and setting for synthetic simulations 
S_PER_HOUR = 3600
MM_PER_M   = 1000

TOPMODEL_exp = 3.0
surfHydCond = 1.0
 
tan_slope = 0.3
soilDepth = 1.5                            # [m] 
porosity  = 0.25                           # [-]
qRain     = 0.002                          # [m hour^-1] 

totalLength1 = 50.                         # [m]
totalLength2 = 10.                         # [m]

hillWidth = 100.                           # [m]
#domain_area = hillWidth * totalLength

hill_dist1 = np.arange(1, 1+totalLength1)  
hill_dist2 = np.arange(1, 1+totalLength2)

kinematic_dist  = np.linspace(0, totalLength1, num=1000)
kinematic_dist2 = np.linspace(0, totalLength2, num=1000)

# construct precipitation
precip= np.zeros(480)
precip[0:240]= 0.002                       #[m hour^-1]

# NB: horizontal conductivity I am not sure how to produce depth averaged hs, then Kh
# but I guess here Martyn in his simulations instead of depth averaged used the surface
# value instead of depth averaged value
kH = surfHydCond * tan_slope               #[m hour^-1]

# time step settings 
#deltat    = 1/(24)                                              # time intervals [days]  
duration  = 20                                                  # time period [days]
# ntsteps   = np.int32((duration/deltat) + 1)                     # number of time steps [-] 
# time      = np.linspace(0,duration, ntsteps)

#%% set variable and parameters required for WATDRN 
# parameters required for WATDRN 
# drainage density 
# this is baed on Soulis et al (2000) and Mekonnen's draft 
Dd1 = 1/(2*totalLength1)                           # [1/m]
Dd2 = 1/(2*totalLength2)                           # [1/m]

# obtain b coefficient  based on Clapp-Hornberger connectivity index  
bij = (TOPMODEL_exp -3)/2                            # [-]

# calculating xlambda explicitly by equating eq.(21) from Martyn's paper and eq.(24) from Mekonnen's draft 
# NB: in MESH implementation, the xlambda decay factor is derived from xdrainh
lambdda = - np.log(1+tan_slope**2)/(soilDepth)     #[1/m]
lambdda = - lambdda                                # to be consistent with fortran code 

# time step 
delt = 1.0                                           # [hour]
#---------------------------------------
# initialize variables
# NB: I assumed soil as one single layer 
ILG = 1                                     # [-] number of tiles 
IG  = 1                                     # [-] number of soil layers
IL1 = 1                                     # [-] index of first tile
IL2 = 1                                     # [-] index of last tile  
IWF = 1                                     # [-] FLAG governing flat or slope CLASS. IWF = 0 activates flat CLASS
ztop   = 0                                  # [m]
delzw  = soilDepth                          # [m]
ksat   = surfHydCond                        # [m hour^-1]
xslope = tan_slope                          # [-]

# Integration of saturated conductivity across the layer -> kl
#xlambda        = -np.log(xdrainh[i])/h0
ktop           = ksat * np.exp(lambdda*ztop)      #  [m hour^-1]

# Important note : I commented calculation of kl using exav functions as it 
# reduces the precision of calculation of grkeff
#kl             = ktop * exav(lambdda*delzw)       #  [m hour^-1]
kl             = ktop * np.exp(lambdda*delzw)       #  [m hour^-1]

# calculate grkeff for two hillslope lengths 
grkeff1         = kl*xslope*2.0*Dd1/(1+xslope**2)  # [hour^-1]
grkeff2         = kl*xslope*2.0*Dd2/(1+xslope**2)  # [hour^-1]

#thpor_avail[i] = np.max((thlmin[i,j],thpor_avail[i]))
thpor_avail = porosity 

time_sim = np.linspace(0, duration, len(precip))

#%% # Set font labels
font = {'family' : 'Times New Roman',
         'weight' : 'bold',
         'size'   : 18}
matplotlib.rc('font', **font)

#%% relative storage based on MoC
def relative_storage(itime, qRain, porosity, soilDepth, kinematic_dist, kH, TOPMODEL_exp):
    # eq.(18)
    Smax = (itime * qRain)/(porosity * soilDepth)
    # rewrite eq.(18) based on equation 20
    S = ((kinematic_dist * qRain)/(kH * soilDepth))**(1/TOPMODEL_exp)
    S[S>Smax]=Smax
    return S

#%% kinematic runoff section 2 (relative storage)
def kinematic_runoff_RS(surfHydCondH, soilDepth, TOPMODEL_exp, porosity,
                                      domain_area, hillWidth, totalLength, hill_dist, qRain):
    
    # NB: This function is written based on m 
    t0 = kH * soilDepth
    qKin = np.ones_like(qRain) * -9999
    
    for i, q in enumerate(qRain):
        if q > 0:
            tc = i+1
            Smax = (tc * q)/(porosity * soilDepth)
            # rewrite eq.(18) based on equation 20
            S = ((hill_dist * q)/(kH * soilDepth))**(1/TOPMODEL_exp)
            S[S>Smax]=Smax
            # domain averaged value 
            qKin[i] = (1/domain_area) * hillWidth * t0 * ((np.max(S)) ** TOPMODEL_exp)
          
    # # reconstruct based on kinematic_runoff formulation        
    # # I dont know whether the derivation of celerity is correct (based on section 2.2.2)
    dxdt = (surfHydCondH * TOPMODEL_exp) * (S) ** (TOPMODEL_exp - 1)/ porosity
    # comes from Martyn's experiment 
    td = 240.
    delX = totalLength - hill_dist
    tExit = td + delX / dxdt
    qFall = (1/domain_area) * hillWidth * t0 * ((S) ** TOPMODEL_exp)
    qKin = qKin[qKin > 0]
    return qKin , qFall, tExit

#%% kinematic runoff recession section 2
def kinematic_runoff_RS_recession(tf, time_res, TOPMODEL_exp, 
                                  totalLength1, soilDepth, surfHydCondH):
    # eq. (22)
    S = (tf/time_res)**(1/(TOPMODEL_exp-1)) 
    
    qs = (1/totalLength1) * soilDepth * surfHydCondH # [mm h^-1]
    qx = qs * S ** TOPMODEL_exp
    return qx, qs

#%% kinematic storage function section 3
def kinematic_storage(surfHydCond, soilDepth, TOPMODEL_exp, porosity, qrain, itime, hill_dist):
    t0 = surfHydCond * soilDepth / TOPMODEL_exp
    # i don't know source of x1, x2. It seems comes from solving eqs 50-53  
    # x1 is obtained when zwat = 0
    x1 = t0 * tan_slope / qRain
    # maybe it comes form dx/dt
    x2 = itime * qRain / (porosity * soilDepth)
    # NB: xL is a threshold here to limit the storage. I guess it is a distance from the top of hillslope 
    # where water table is parallel to slope
    xL = min(x1 * x2 ** TOPMODEL_exp, hill_dist[-1])
    
    # calculate hs from inversion of equation 56 where transmisivity is calculated from eq 59
    hS = soilDepth * ((qRain * hill_dist) / (t0 * tan_slope)) ** (1./TOPMODEL_exp)
    hL = soilDepth * ((qRain * xL) / (t0 * tan_slope)) ** (1./TOPMODEL_exp)
    hS[hS>hL] = hL
    return hS

#%% kinematic_runoff section 3
def kinematic_runoff(surfHydCond, soilDepth, TOPMODEL_exp, porosity, domain_area, hillWidth, totalLength, hill_dist, qRain):
    t0 = surfHydCond * soilDepth / TOPMODEL_exp
    hS = soilDepth * (qRain[0] * hill_dist / (t0 * tan_slope)) ** (1/TOPMODEL_exp)
    qKin = np.ones_like(qRain) * -9999
    # xL2 = np.zeros(len(qRain))
    # hL2= np.zeros(len(qRain))
    for i, q in enumerate(qRain):
        if q > 0:
            tc = i+1
            x1 = t0 * tan_slope / q
            # tc here is time 
            x2 = tc * q / (porosity * soilDepth)
            xL = min(x1*x2 ** TOPMODEL_exp, hill_dist[-1])
            hL = soilDepth * (q * xL / (t0 * tan_slope)) ** (1/TOPMODEL_exp)
            # xL2[i] = xL
            # hL2[i] = hL
            # does xlen equals to (1/domain_area) * hillWidth 
            # when hL is calculated, then from eq.(56) you can calculate the outflow
            qKin[i] = (1/domain_area) * hillWidth * tan_slope * t0 * ((hL / soilDepth) ** TOPMODEL_exp)
    
    dtdz = (tan_slope * t0 / soilDepth) * TOPMODEL_exp * (hS / soilDepth) ** (TOPMODEL_exp - 1)
    dxdt = dtdz / porosity
    # comes from Martyn's experiment 
    td = 240.
    delX = totalLength - hill_dist
    tExit = td + delX / dxdt
    qFall = (1/domain_area) * hillWidth * tan_slope * t0 * ((hS/soilDepth) ** TOPMODEL_exp)
    qKin = qKin[qKin > 0]
    return qKin, qFall, tExit 

#%% relative storage recession function 
def relative_storage_recession(kinematic_dist, porosity, itime, kH, TOPMODEL_exp):
    # eq.(22)
    S = ((kinematic_dist*porosity)/(itime*kH*TOPMODEL_exp))**(1/(TOPMODEL_exp-1))
    # constrain values more than 1
    S[S>1] = 1
    return S

#%% calculate fractional outflow section 2.5
def fractional_outflow (S, Sc , TOPMODEL_exp):  
    qx_qs = (S/Sc) ** TOPMODEL_exp
    # constrain values when soil is saturated 
    qx_qs[qx_qs>1] = 1
    
    return qx_qs

#%% WATDRN area-average storage-discharge relationship for saturated soil 
# section 3.4 Martyn's paper   
def WATDRN_saturated_runoff(q0s, zwt, nm, zwt_crit, TOPMODEL_exp):
    # [m], eq.(67)
    qx           =  q0s * ((nm -zwt)/(nm - zwt_crit)) ** TOPMODEL_exp  
    # nm -zwt equal to water table height
    # constrain the qx for values larger than the critical value 
    qx[qx>q0s]   =  q0s                                         
    return qx 

#%% TOPMODEL area-averaged storage-discharge relationship  for saturated soil 
def TOPMODEL_saturated_runoff(totalLength, zwt, 
                              kH, m, TOPMODEL_exp, nm):
    # eq.(73)
    lambdda = ((totalLength)/(kH * m)) ** (1/TOPMODEL_exp) * (TOPMODEL_exp/(TOPMODEL_exp +1)) 
    # eq.(A.15)
    qx = MM_PER_M * 1/(lambdda) ** TOPMODEL_exp * (1 - zwt/nm) ** TOPMODEL_exp
    return qx

#%% Soulis's upscaling outflow (section 2.5)
# NB: This function is based on Martyn's solution    
def Soulis_upscaling_outflow (qForce, kH, xL, Scrit, zSoil, phi, eta):
    
    #get the hillslope outflow when the seepage face is fully saturated
    q0 = kH*zSoil/xL
    
    #define the time step
    nForce = len(qForce)
    dt     = 1
    nSub   = 100
    xTime  = 0
    tVec   = np.zeros(nForce)
    
    #initialize (mx points distributed across the full space)
    S = 0
    
    #save variables
    sSave = np.zeros(nForce)
    qSave = np.zeros(nForce)
    
    # loop through points on the characteristic curve 
    for iForce in range(nForce):
        
        #get the length of the substep
        dtSub = dt/nSub
        
        # numerically integrate over time
        for iSub in range(nSub):
            # preliminaries
            xTime = xTime + dtSub  # increment time
    
            # solution
            # eq. (35)
            qb  = q0*(S/Scrit)**eta
            #  eq. (16)
            rhs = (qForce[iForce] - qb)/(zSoil*phi)
            S = S + rhs*dtSub
        
        # save data
        sSave[iForce]   = S
        qSave[iForce]   = qb
        tVec[iForce]    = xTime 
    
    return sSave, qSave, tVec

#%% Soulis's upscaling outflow for saturated flow problems 
def Soulis_upscaling_outflow_saturated (qForce, kH, xL, Zcrit, zSoil, phi, nm, m, n):
    
    #get the hillslope outflow when the seepage face is fully saturated
    q0s = kH*m/xL
    
    #define the time step
    nForce = len(qForce)
    dt     = 1
    nSub   = 100
    xTime  = 0
    tVec   = np.zeros(nForce)
    
    #initialize (mx points distributed across the full space)
    Z = zSoil
    
    #save variables
    zSave = np.zeros(nForce)
    qSave = np.zeros(nForce)
    
    # loop through points on the characteristic curve 
    for iForce in range(nForce):
        
        #get the length of the substep
        dtSub = dt/nSub
        
        # numerically integrate over time
        for iSub in range(nSub):
            # preliminaries
            xTime = xTime + dtSub  # increment time
    
            # solution
            # eq.(67)
            qb  = q0s*((nm -Z)/(nm - Zcrit))**n
            # eq.(45) 
            rhs = -(qForce[iForce] - qb)/(phi)
            Z = Z + rhs*dtSub
        
        # save data
        zSave[iForce]   = Z
        qSave[iForce]   = qb
        tVec[iForce]    = xTime 
    
    return zSave, qSave, tVec

#%% TOPMODEL upscaling outflow for saturated flow problems 
def TOPMODEL_upscaling_outflow_saturated (qForce, kH, xL, zSoil, phi, nm, m, n):
    
    # eq.(73)
    lambdda = ((xL)/(kH * m)) ** (1/n) * (n/(n +1)) 
     
    #define the time step
    nForce = len(qForce)
    dt     = 1
    nSub   = 100
    xTime  = 0
    tVec   = np.zeros(nForce)
    
    #initialize (mx points distributed across the full space)
    Z = zSoil
    
    #save variables
    zSave = np.zeros(nForce)
    qSave = np.zeros(nForce)
    
    # loop through points on the characteristic curve 
    for iForce in range(nForce):
        
        #get the length of the substep
        dtSub = dt/nSub
        
        # numerically integrate over time
        for iSub in range(nSub):
            # preliminaries
            xTime = xTime + dtSub  # increment time
    
            # solution
            # # eq.(A.15)
            qb  = 1/(lambdda) ** n * (1 - Z/nm) ** n
            # eq.(45) 
            rhs = -(qForce[iForce] - qb)/(phi)
            Z = Z + rhs*dtSub
        
        # save data
        zSave[iForce]   = Z
        qSave[iForce]   = qb
        tVec[iForce]    = xTime 
    
    return zSave, qSave, tVec
#%% WATDRN storage 
def WATDRN_storage (qForce, kH, xL, zSoil, phi, eta, dt):
    
    # this is based on native WATDRN algorithm 
    #define the time step
    nForce = len(qForce)
    
    #initialize (mx points distributed across the full space)
    S0 = 1.e-8
    
    #save variables
    sWATDRN = np.zeros(nForce)
    qWATDRN = np.zeros(nForce)
    
    # get temporally constant variables
    tf = phi*xL/(kH*eta) # time that the drying front reaches the bottom of the hillslope
    Sc = 1.0 - 1.0/eta   # spatially-averaged relatve storage at time tf
    
    # loop through points on the characteristic curve 
    for iForce in range(nForce):
        
        # compute the time required for spatially-averaged storage to reach the state from the host land model
        # - only consider case 2: the toe of the hillslope is not fully saturated
        tr = tf*(S0/Sc)**(1.0 - eta)
        
        # compute the storage at the time tr + Delta t
        t1 = tr + dt
        S1 = Sc*(tf/t1)**(1.0/(eta - 1.0))
        
        # compute the lateral flow from mass balance
        qx = phi*zSoil*(S0 - S1)/dt
        #print, iForce, S0, S1, [qSave[iForce], qx]*1000.d
        
        #update states
        Sn = S0 + (qForce[iForce] - qx)/(zSoil*phi)
        S0 = Sn
        
        # save data
        sWATDRN[iForce]   = S0
        qWATDRN[iForce]   = qx
    
    return sWATDRN, qWATDRN  

#%% WATDRN routine 
# NB: as this experiment is run over one tile, the loop over tiles has been removed  
def WATDRN(delzw,bcoef,thpora,grksat,grkeff,asat0,iwf, \
     		ilg,il1,il2,wfk,delt):
                   
   # output variables 
   # asat1(ilg)  ! bulk saturation at the end of the time step
   # subflw(ilg) ! interflow amount during the time step (m)
   # basflw(ilg) ! baseflow rate during the time step (m)
   # satsfc(ilg) ! saturated fraction of the surface (0 to 1)
    
    if (not(np.all(iwf == wfk))):
        exit
    
    # **********************************************************************
    #      STEP 0: Initialize a few things before we start
    #              - output variables
    #              - functions of the input parameters
    # **********************************************************************
    # c      = np.zeros(ilg)            # Clapp-Hornberger connectivity index (c>1)
    # cm1    = np.zeros(ilg)            # c-1
    # c2m1   = np.zeros(ilg)            # 2*c-1
    # asatc  = np.zeros(ilg)            # bulk saturation at the critical time tc
    # asat00 = np.zeros(ilg)            # MIN(1,asat0) because we don't handle supersaturated soils
    # tc     = np.zeros(ilg)            # critical time at which the seepage face becomes unsaturated
    # ratiot = np.zeros(ilg)            # ratio tc/t (t0 or t1) if t>tc and t/tc if t<=tc
    # satspf = np.ones(ilg, dtype=bool) # indicates if seepage face is saturated
    #             					  # equivalent to knowing if t<=tc 
    
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    
    #c and c factors
    c    = 2.*bcoef+3.
    cm1  = c-1.
    c2m1 = 2.*c-1.
    
    # bulk saturation at critical time
    #(just before the seepage face becomes unsaturated)
    asatc = 1.-1./c

    # layer average saturation asat0 may be greater than 1
    #e.g. frost heave but it is not possible for wat_drain
    asat00 = np.min((1., asat0))
    #assess if seepage face is saturated at initial time
    satspf = asat00 >= asatc
        
    # **********************************************************************
    #      STEP 1: Find theoretical time t0 elapsed since element was last
    #              saturated and estimate baseflow rate at initial time
    #              Also estimate fraction of surface that is saturated
    # ********************************************************************** 
    # tc = np.zeros(ilg)     # critical time at which the seepage face becomes unsaturated
    # ratiot = np.zeros(ilg) # ratio tc/t (t0 or t1) if t>tc and t/tc if t<=tc
    # basflw = np.zeros(ilg) # baseflow rate during the time step (m)
    # satsfc = np.zeros(ilg) # saturated fraction of the surface (0 to 1)
    
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
        # determine time at which seepage face becomes unsaturated
    tc = thpora/(c*grkeff)
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    #find theoretical start of recession (t0) from bulk saturation
    # and at the same time estimate baseflow based on rate at t0
    if (satspf):
      #saturated seepage face at initial time:
      # compute t0/tc
      ratiot = c*(1.-asat00)
      # normalized baseflow rate
      basflw = 1.-c*c/c2m1*(1.-asat00)
      #the fraction of the surface that is saturated at t0
      #varies linearly with t0/tc
      satsfc= 1.-ratiot
    else:
      #unsaturated seepage face at initial time:
      #calculate tc/t0 instead of t0 to avoid overflow 
      ratiot = (asat00/asatc)**cm1
      #normalized baseflow rate
      basflw = cm1/c2m1*ratiot*asat00/asatc 
      # the fraction of the surface that is saturated at t0 is zero
      satsfc = 0.
            
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    #Compute baseflow in m from normalized baseflow rate
    basflw = grksat*basflw*delt
        
    # **********************************************************************
    #      STEP 2: Find theoretical time t1 at the end of the time step
    # **********************************************************************
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    if (satspf):
        # Assess if seepage face will still be saturated at the
        #  end of the time step
        satspf = tc * ratiot + delt <= tc
        if (satspf):
            # Seepage face still saturated, compute t1/tc from t0/tc
            ratiot= (tc*ratiot+delt)/tc
        else:
            # Seepage face not saturated anymore, compute tc/t1
            ratiot = tc/(tc*ratiot+delt)
    else:
        # If seepage face was not saturated initially, we compute
        # tc/t1=tc/(t0+delt) from tc/t0
        ratiot = tc*ratiot/(tc+delt*ratiot)
            
    # **********************************************************************
    #      STEP 3: Obtain bulk saturation at the end of the time step
    #              and interflow amount
    # **********************************************************************  
    # asat1  = np.zeros(ilg)
    # subflw = np.zeros(ilg) 
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    if (satspf):
       # saturated seepage face at the end of the time step
       asat1 = 1.-ratiot/c
    else:
       # unsaturated seepage face at the end of the time step
       asat1= asatc*ratiot**(1./cm1)
    
    # for i in range(il1-1,il2):
    #     # cycle if not using watrof or latflow
    #     if (iwf[i] != wfk or asat0[i] == 0.0):
    #         continue
    # Sanity check: bulk saturation should not increase with time
    asat1 = np.min((asat00,asat1))
    # Obtain interflow from the difference in bulk saturation
    subflw = (asat00-asat1)*thpora*delzw
    
    return asat1,subflw,basflw,satsfc

#%% exav
#**********************************************************************
# 
#      function exav(x)
#      finds average of an exponential function exp(-x)
#      over the interval [0,x], which is (1-exp(-x))/x
#      deals with limit values of 1-x/2 when x->0 and 1/x when x->inf
# 
# **********************************************************************
def exav(x):
    # initialize variables 
    exphuge  = 1.0e+9
    expbig   = 1.0e+8
    expsmall = 1.0e-8
    
    if (x > exphuge):
        exav = 0.0
    elif (x > expbig):
        exav = 1.0/x  
    elif (x > expsmall):
        exav = (1.0-np.exp(-x))/x
    else: 
        exav = 1.0-x/2.0
    return exav

#%% Figure 3- Rainfall pulse 
itime1 = 24
itime2 = 2*24
itime3 = 3*24
itime4 = 4*24

S1 = relative_storage(itime1, qRain, porosity, soilDepth, kinematic_dist, kH, TOPMODEL_exp)
S2 = relative_storage(itime2, qRain, porosity, soilDepth, kinematic_dist, kH, TOPMODEL_exp)
S3 = relative_storage(itime3, qRain, porosity, soilDepth, kinematic_dist, kH, TOPMODEL_exp)
S4 = relative_storage(itime4, qRain, porosity, soilDepth, kinematic_dist, kH, TOPMODEL_exp)

# calling kinematic_runoff based on relative storage for rising and falling limb 
qKin, qFall, tExit=  kinematic_runoff_RS(kH, soilDepth, TOPMODEL_exp, porosity, hillWidth * totalLength1, 
                                          hillWidth, totalLength1, hill_dist1,  precip[0:240])


fig,axs = plt.subplots(1,2, figsize=(20, 20))

axs[0].plot(kinematic_dist, S1, '#00007fe6', alpha = 0.2, label='1 day')
axs[0].plot(kinematic_dist, S2, '#00007fe6', alpha = 0.4, label='2 day')
axs[0].plot(kinematic_dist, S3, '#00007fe6', alpha = 0.6, label='3 day')
axs[0].plot(kinematic_dist, S4, '#00007fe6', alpha = 0.8, label='4 day')

axs[0].legend(fontsize = 14, loc = 'upper left',frameon=False)
axs[0].grid(alpha = 0.5)


# set axis 
axs[0].set_xlabel('x (m)')
axs[0].set_ylabel('Relative storage (-)')

axs[0].set_xlim(0,np.max(kinematic_dist))
axs[0].set_ylim(0,1)

# presenting rising and falling hydrograph 
time_rise  = np.linspace(0,10,240) # [days]
axs[1].plot(time_rise, MM_PER_M * qKin, '#00007fe6', label='Method of characteristics')
axs[1].plot(tExit/24, MM_PER_M * qFall, '#00007fe6') 
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)
# set axis 
# set axis 
axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')

axs[1].set_xlim(0,20)
axs[1].set_ylim(0,2.5)
axs[1].grid(alpha=0.5)

# save figure 
plt.savefig(outdir+'Figure3.png', format='png', dpi=300)
plt.close()

#%% reproducing figure 4
# NB: time should be presented in hours as the unity of 
# hydraulic conductivity is m h^-1  

fig,axs = plt.subplots(1,2, figsize=(20, 20))

itime1 = 0.5*24  # [hours]
itime2 = 2*24    
itime3 = 5*24  
itime4 = 10*24

# obtain relative sotrage for different time steps (left side of plot)
S1 = relative_storage_recession(kinematic_dist, porosity, itime1, kH, TOPMODEL_exp)
S2 = relative_storage_recession(kinematic_dist, porosity, itime2, kH, TOPMODEL_exp)
S3 = relative_storage_recession(kinematic_dist, porosity, itime3, kH, TOPMODEL_exp)
S4 = relative_storage_recession(kinematic_dist, porosity, itime4, kH, TOPMODEL_exp)

# calculating the time that the drying front reaches the bottom of hillslope 
# eq. (21)
tf = np.round((porosity*totalLength1)/(kH*TOPMODEL_exp),2)   # tf[hours]  
# calculate relative storage and discharge during the recession period 
time_dry = np.linspace(0, tf ,10)     # [hours]
time_res = np.linspace(tf, 5*24 ,100) # [hours] 

# calling recession runoff function  
qx,qs = kinematic_runoff_RS_recession(tf, time_res, TOPMODEL_exp, totalLength1, soilDepth, kH)

# plot results 
axs[0].plot(kinematic_dist, S1, color = '#00007fe6', alpha = 0.9, label='Method of characteristics')
axs[0].plot(kinematic_dist, S2, color = '#00007fe6', alpha = 0.7)
axs[0].plot(kinematic_dist, S3, color = '#00007fe6', alpha = 0.5)
axs[0].plot(kinematic_dist, S4, color = '#00007fe6', alpha = 0.3)
axs[0].legend(fontsize = 14, loc = 'upper left',frameon=False)

# set annotation
axs[0].text(4.0,0.50, '0.5 days')
axs[0].text(40,0.52, '2 days')
axs[0].text(43,0.34, '5 days')
axs[0].text(43,0.24, '10 days')

# set axis 
axs[0].set_xlabel('x (m)')
axs[0].set_ylabel('Relative storage (-)')

axs[0].set_xlim(0,np.max(kinematic_dist))
axs[0].set_ylim(0,1)
axs[0].grid(alpha=0.5)

axs[1].plot(time_dry/24, MM_PER_M * qs*np.ones(len(time_dry)), color = '#00007fe6',label='Method of characteristics')
axs[1].plot(time_res/24, MM_PER_M * qx, color = '#00007fe6')
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)

# set axis 
axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')

axs[1].set_xlim(0,5)
axs[1].set_ylim(0,10)
axs[1].grid(alpha=0.5)

# save the output 
plt.savefig(outdir+'Figure4.png', format='png', dpi=300)
plt.close()

#%% Figure 5
# calculate the spatially-averaged relative saturation at the critical time 

# intialize spatially-averaged relative saturation 
S = np.linspace(0, 1, num=100)

# critical value of spatially-averaged relative saturation based on eq.(27)
Sc = 1 - 1/TOPMODEL_exp 

qx_qs = fractional_outflow (S, Sc , TOPMODEL_exp)

#----------------------------------------------------------------------------
# calling kinematic_runoff based on relative storage for rising and falling limb 
# condition 1
qKin1, qFall1, tExit1=  kinematic_runoff_RS(kH, soilDepth, TOPMODEL_exp, porosity, hillWidth * totalLength1, 
                                          hillWidth, totalLength1, kinematic_dist,  precip[0:240])

# condition 2
qKin2, qFall2, tExit2=  kinematic_runoff_RS(kH, soilDepth, TOPMODEL_exp, porosity, hillWidth * totalLength2, 
                                          hillWidth, totalLength2, kinematic_dist2,  precip[0:240])

# simulating outflow based on Soulis's upscaling method
S_Soulis1, qx_soulis1, time_sim = Soulis_upscaling_outflow (precip, kH, totalLength1, Sc, soilDepth, porosity, TOPMODEL_exp)

S_Soulis2, qx_soulis2, time_sim = Soulis_upscaling_outflow (precip, kH, totalLength2, Sc, soilDepth, porosity, TOPMODEL_exp)


# calling simulated WATDRN storage required for running Stand-alone WATDRN 
sWATDRN1, qWATDRN1 = WATDRN_storage (precip, kH, totalLength1, soilDepth, porosity, TOPMODEL_exp, delt)
sWATDRN2, qWATDRN2 = WATDRN_storage (precip, kH, totalLength2, soilDepth, porosity, TOPMODEL_exp, delt)

#---------------------------------------
# calling WATDRN routine  
# variables related to WATDRN routine 
p = len(S_Soulis1)

asat_t1  = np.zeros(p) 
subflw1  = np.zeros(p) 
basflw1  = np.zeros(p)
satfc1   = np.zeros(p)

asat_t2  = np.zeros(p) 
subflw2  = np.zeros(p) 
basflw2  = np.zeros(p)
satfc2   = np.zeros(p)

# compute interflow from the layer (subflowj). Baseflow from the layer (basflwj) is 
# also computed but is not used at present.

# loop over entire time window of simulation for two hillslope lengths 
for i in range(p):
    # hillslope with length of totalLength1
    [asat_t1[i], subflw1[i], basflw1[i], satfc1[i]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff1,sWATDRN1[i], IWF,\
                                                             ILG,IL1,IL2,1,delt)
    
    # NB: in the source code of WATDRN the lateral flow is not devided by time step
    subflw1[i] = MM_PER_M * subflw1[i]/delt                # convert to [mm hour^-1]

    
    # hillslope with length of totalLength2    
    [asat_t2[i], subflw2[i],basflw2[i],satfc2[i]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff2,sWATDRN2[i],IWF,\
                                                            ILG,IL1,IL2,1,delt)
    
    subflw2[i] = MM_PER_M * subflw2[i]/delt               # convert to [mm hour^-1]

#---------------------------------------
# plot results 
# to be in same range of MoC
time_sim = (time_sim-1)/len(precip)*duration

# plot the results Figure 5.1 
figure,axs = plt.subplots(1,3, figsize= (20,20))

axs[0].plot(S, qx_qs, linewidth=3)

# set axes, labels, grids  
axs[0].set_xlabel('Relative storage (-)')
axs[0].set_ylabel('Fractional outflow (-)')

axs[0].set_xlim(np.min(S), np.max(S))
axs[0].set_ylim(np.min(qx_qs), np.max(qx_qs)) 

axs[0].grid(alpha=0.5)

# presenting rising and falling hydrograph 
time_rise  = np.linspace(0,10,240) # [days]
axs[1].plot(time_rise, MM_PER_M * qKin1, '#00007fe6', label='Method of characteristics')
axs[1].plot(tExit1/24, MM_PER_M * qFall1, '#00007fe6') 

axs[1].plot(time_sim, MM_PER_M * qx_soulis1, '#47b1e2', label='Soulis upscaling method')

axs[1].plot(time_sim, subflw1, '#ff0000', linestyle = '--',label='WATDRN Stand-alone')

axs[1].text(14.8, 1.05, '$X_{L}$ = 50 m')
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)

axs[2].plot(time_rise, MM_PER_M * qKin2, '#00007fe6', label='Method of characteristics')
axs[2].plot(tExit2/24, MM_PER_M * qFall2, '#00007fe6') 

axs[2].plot(time_sim, MM_PER_M * qx_soulis2, '#47b1e2', label='Soulis upscaling method')

axs[2].plot(time_sim, subflw2, '#ff0000', linestyle = '--', label='WATDRN Stand-alone')

axs[2].text(14.8, 1.05, '$X_{L}$ = 10 m')
axs[2].legend(fontsize = 14, loc = 'upper left',frameon=False)
 
# set axis 
axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')
axs[2].set_xlabel('Time (days)')

axs[1].set_xlim(0,duration)
axs[1].set_ylim(0,2.5)
axs[1].grid(alpha=0.5)

axs[2].set_xlim(0,duration)
axs[2].set_ylim(0,2.5)
axs[2].grid(alpha=0.5)

plt.savefig(outdir+'Figure5_WATDRN.png', format='png', dpi=300)
plt.close()
       
#%% figure 6 Behavior of Soulis’ upscaling method applied to saturated flow problems  
# calculate hillslope outflow based on water table depth 
nm         = soilDepth                                        #[m]
m          =  soilDepth/TOPMODEL_exp                          #[m]
q0s1       = MM_PER_M * kH * m / totalLength1    #[mm hour^-1], eq. (55) 
q0s2       = MM_PER_M * kH * m / totalLength2    #[mm hour^-1], eq. (55) 
zwt_crit   = m                                                #[m], eq.(60)
zwt        = np.linspace(1.5, 0,100)    

qx1 = WATDRN_saturated_runoff(q0s1, zwt, nm, zwt_crit, TOPMODEL_exp)
qx2 = WATDRN_saturated_runoff(q0s2, zwt, nm, zwt_crit, TOPMODEL_exp)

# NB: MoC figures 6 & 7 are same
#calling kinematic_runoff for generating otflow for MOC 
qKin1, qFall1, texit1 = kinematic_runoff(surfHydCond, soilDepth, TOPMODEL_exp, porosity, hillWidth * totalLength1, 
                                         hillWidth, totalLength1, kinematic_dist,  precip)
 
qKin2, qFall2, texit2 = kinematic_runoff(surfHydCond, soilDepth, TOPMODEL_exp, porosity, hillWidth * totalLength2, 
                                         hillWidth, totalLength2, kinematic_dist2,  precip)

# simulating outflow based on Soulis's upscaling method for saturated flow 
Zwt_Soulis1, qx_soulis1, time_sim = Soulis_upscaling_outflow_saturated (precip, kH, totalLength1, zwt_crit, 
                                                                        soilDepth, porosity, nm, m, TOPMODEL_exp)

Zwt_Soulis2, qx_soulis2, time_sim = Soulis_upscaling_outflow_saturated (precip, kH, totalLength2, zwt_crit, 
                                                                        soilDepth, porosity, nm, m, TOPMODEL_exp)

# to be in same range of MoC
time_sim = (time_sim-1)/len(precip)*duration

fig,axs = plt.subplots(1,3, figsize=(20, 20))

axs[0].plot(zwt, qx1, color='#00007fe6', label='$X_{L}$ = 50 m')
axs[0].plot(zwt, qx2, color='#47b1e2', label='$X_{L}$ = 10 m')
axs[0].legend(fontsize = 14, loc = 'upper left', frameon=False) 
axs[0].grid(alpha = 0.5)

axs[1].plot(time_rise, MM_PER_M * qKin1, color='#00007fe6', label='Method of characteristics')
axs[1].plot(texit1/24, MM_PER_M * qFall1, color='#00007fe6')

axs[1].plot(time_sim, MM_PER_M * qx_soulis1, '#47b1e2', label='Soulis upscaling method')

axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False) 
axs[1].text(14.7, 1.05, '$X_{L}$ = 50 m')
axs[1].grid(alpha = 0.5)

axs[2].plot(time_rise, MM_PER_M * qKin2, color='#00007fe6', label='Method of characteristics')
axs[2].plot(texit2/24, MM_PER_M * qFall2, color='#00007fe6')

axs[2].plot(time_sim, MM_PER_M * qx_soulis2, '#47b1e2', label='Soulis upscaling method')

axs[2].legend(fontsize = 14, loc = 'upper left',frameon=False) 
axs[2].text(14.7, 1.05, '$X_{L}$ = 10 m')
axs[2].grid(alpha = 0.5)

# set axes and labels 
axs[0].set_ylabel('Hillslope outflow [mm/hour]')
axs[1].set_ylabel('Hillslope outflow [mm/hour]')
#axs[2].set_ylabel('Runoff [mm/hour]')
axs[0].set_xlabel('Water table depth [m]')
axs[1].set_xlabel('Time (days)')
axs[2].set_xlabel('Time (days)')

axs[0].set_xlim(1.5,0)
axs[0].set_ylim(0,20)

axs[1].set_xlim(0,duration)
axs[1].set_ylim(0,2.5)

axs[2].set_xlim(0,duration)
axs[2].set_ylim(0,2.5)

plt.savefig(outdir+'Figure6.png', format='png', dpi=300)
plt.close()

#%% Figure 7 TOPMODEL spatially-averaged storage-discharge relationship 

qx1 = TOPMODEL_saturated_runoff(totalLength1, zwt, kH, m, TOPMODEL_exp, nm)
qx2 = TOPMODEL_saturated_runoff(totalLength2, zwt, kH, m, TOPMODEL_exp, nm)

# simulating outflow based on Soulis's upscaling method for saturated flow
Zwt_TOPMODEL1, qx_TOPMODEL1, time_sim = TOPMODEL_upscaling_outflow_saturated (precip, kH, totalLength1, 
                                                                               soilDepth, porosity, nm, m, TOPMODEL_exp)
Zwt_TOPMODEL2, qx_TOPMODEL2, time_sim = TOPMODEL_upscaling_outflow_saturated (precip, kH, totalLength2, 
                                                                               soilDepth, porosity, nm, m, TOPMODEL_exp)
# to be in same range of MoC
time_sim = (time_sim-1)/len(precip)*duration

fig,axs = plt.subplots(1,3, figsize=(20, 20))

axs[0].plot(zwt, qx1, color='#00007fe6', label='$X_{L}$ = 50 m')
axs[0].plot(zwt, qx2, color='#47b1e2', label='$X_{L}$ = 10 m')
axs[0].legend(fontsize = 14, loc = 'upper left', frameon=False) 
axs[0].grid(alpha = 0.5)

axs[1].plot(time_rise, MM_PER_M * qKin1, color='#00007fe6', label='Method of characteristics')
axs[1].plot(texit1/24, MM_PER_M * qFall1, color='#00007fe6')

axs[1].plot(time_sim, MM_PER_M * qx_TOPMODEL1, '#47b1e2', label='TOPMODEL')

axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False) 
axs[1].text(14.7, 1.05, '$X_{L}$ = 50 m')
axs[1].grid(alpha = 0.5)

axs[2].plot(time_rise, MM_PER_M * qKin2, color='#00007fe6', label='Method of characteristics')
axs[2].plot(texit2/24, MM_PER_M * qFall2, color='#00007fe6')

axs[2].plot(time_sim, MM_PER_M * qx_TOPMODEL2, '#47b1e2', label='TOPMODEL')

axs[2].legend(fontsize = 14, loc = 'upper left',frameon=False) 
axs[2].text(14.7, 1.05, '$X_{L}$ = 10 m')
axs[2].grid(alpha = 0.5)

# set axes and labels 
axs[0].set_ylabel('Hillslope outflow [mm/hour]')
axs[1].set_ylabel('Hillslope outflow [mm/hour]')
#axs[2].set_ylabel('Runoff [mm/hour]')
axs[0].set_xlabel('Water table depth [m]')
axs[1].set_xlabel('Time (days)')
axs[2].set_xlabel('Time (days)')

axs[0].set_xlim(1.5,0)
axs[0].set_ylim(0,20)

axs[1].set_xlim(0,duration)
axs[1].set_ylim(0,2.5)

axs[2].set_xlim(0,duration)
axs[2].set_ylim(0,2.5)

plt.savefig(outdir+'Figure7.png', format='png', dpi=300)
plt.close()