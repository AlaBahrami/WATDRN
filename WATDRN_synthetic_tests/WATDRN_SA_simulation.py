# -*- coding: utf-8 -*-
"""
The purpose of this script is to run the stand-alone version of WATDRN and do some
experiments based on paper of WATDRN critique. 
Created on Mon Mar 13 09:17:22 2023

@author: Ala Bahrami

See also : 
          WATDRN_SA.py and watdrn_simulation.py 
"""
#%% import modules 
import matplotlib
import matplotlib.pyplot as plt
#import xarray as xs 
import numpy as np

#%% I/O directory 
outdir           = 'output/'
 
#%% input parameters and setting for synthetic simulations 
MM_PER_M     = 1000                        # [-]
S_PER_HOUR   = 3600                        # [-]

TOPMODEL_exp = 3.0                         # [-]
surfHydCond  = 1.0                         # [m hour^-1]
 
tan_slope    = 0.3                         # [-]
soilDepth    = 1.5                         # [m] 
porosity     = 0.25                        # [-]
qRain        = 0.002                       # [m hour^-1] 

totalLength1 = 50.                         # [m]
totalLength2 = 10.                         # [m]

hillWidth    = 100.                        # [m]
#domain_area = hillWidth * totalLength

hill_dist1 = np.arange(1, 1+totalLength1)  
hill_dist2 = np.arange(1, 1+totalLength2)

kinematic_dist  = np.linspace(0, totalLength1, num=1000)
kinematic_dist2 = np.linspace(0, totalLength2, num=1000)

# construct precipitation
precip= np.zeros(480)
precip[0:240]= 0.002                       #[m hour^-1]

# Horizontal conductivity 
kH = surfHydCond * tan_slope               #[m hour^-1]

# time duration of simulation  
duration  = 20                              #[days]            

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

#%% simulate relative storage and outflow based on Soulis's upscaling method

# critical value of spatially-averaged relative saturation based on eq.(27)
Sc = 1 - 1/TOPMODEL_exp

S_Soulis1, qx_soulis1, time_sim = Soulis_upscaling_outflow (precip, kH, totalLength1, Sc, soilDepth, porosity, TOPMODEL_exp)

S_Soulis2, qx_soulis2, time_sim = Soulis_upscaling_outflow (precip, kH, totalLength2, Sc, soilDepth, porosity, TOPMODEL_exp)

#%% # Set font labels
font = {'family' : 'Times New Roman',
         'weight' : 'bold',
         'size'   : 18}
matplotlib.rc('font', **font)

#%% calculate interflow based on WATDRN routine 
# NB: the purpose of this script is produce lateral flow simulation. The relative
# storage is simulated from synthetic test.  
#---------------------------------------
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

# parameter - will be used to compute xdrainh (the fractional change in horizontal
# conductivity in a depth change h0) in Vincent's new formula.
h0 = 1.0 

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

#---------------------------------------
# compute interflow from the layer (subflowj). Baseflow from the layer (basflwj) is 
# also computed but is not used at present.

# loop over entire time window of simulation for two hillslope lengths 
for i in range(p):
    # hillslope with length of totalLength1
    [asat_t1[i], subflw1[i], basflw1[i], satfc1[i]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff1,S_Soulis1[i], IWF,\
                                                             ILG,IL1,IL2,1,delt)
    
    # NB: in the source code of WATDRN the lateral flow is not devided by time step
    subflw1[i] = MM_PER_M * subflw1[i]/delt                # convert to [mm hour^-1]

    
    # hillslope with length of totalLength2    
    [asat_t2[i], subflw2[i],basflw2[i],satfc2[i]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff2,S_Soulis2[i],IWF,\
                                                            ILG,IL1,IL2,1,delt)
    
    subflw2[i] = MM_PER_M * subflw2[i]/delt               # convert to [mm hour^-1]

# plot the results referri
time_sim = np.linspace(0, duration, len(subflw1))    

figure,axs = plt.subplots(1,2, figsize= (20,20))

axs[0].plot(time_sim, subflw1, '#ff0000', label='(-lambda) = %f'%lambdda)
axs[0].text(14.8, 1.05, '$X_{L}$ = 50 m')
axs[0].legend(fontsize = 14, loc = 'upper left',frameon=False)

axs[1].plot(time_sim, subflw2, '#ff0000', label='(-lambda) = %f'%lambdda)
axs[1].text(14.8, 1.05, '$X_{L}$ = 10 m')
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)

# set axes and title 
axs[0].set_title('synthetic WATDRN simulaiton')
axs[1].set_title('synthetic WATDRN simulaiton')

axs[0].set_xlabel('Time (days)')
axs[1].set_xlabel('Time (days)')

axs[0].set_ylabel('Hillslope outflow (mm $h^-1$)')
axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')

axs[0].set_xlim(0,duration)
axs[0].set_ylim(0,2.5)
axs[0].grid(alpha=0.5)

axs[1].set_xlim(0,duration)
axs[1].set_ylim(0,2.5)
axs[1].grid(alpha=0.5)

plt.savefig(outdir+'WATDRN_synthetic.png', format='png', dpi=300)
plt.close()

#%% compare WATDRN stand-alone with WATDRN derived from eqs(32) and eqs(33) of paper
qx1 = np.zeros(p)
qx2 = np.zeros(p)

for i in range(p):
    # for Hillslope with totalLength1 
    # eq.(33)
    tr = ((porosity*totalLength1)/(kH*TOPMODEL_exp))*(S_Soulis1[i]/Sc)**(1-TOPMODEL_exp)
    tq = tr + delt
    # eq.(32)
    S2 = Sc * ((porosity*totalLength1)/(kH*TOPMODEL_exp*tq))**(1/(TOPMODEL_exp-1))
    qx1[i] = MM_PER_M * np.abs(S2 - S_Soulis1[i]) * porosity * soilDepth/delt
    
    # for Hillslope with totalLength2
    # eq.(33)
    tr = ((porosity*totalLength2)/(kH*TOPMODEL_exp))*(S_Soulis2[i]/Sc)**(1-TOPMODEL_exp)
    tq = tr + delt
    # eq.(32)
    S2 = Sc * ((porosity*totalLength2)/(kH*TOPMODEL_exp*tq))**(1/(TOPMODEL_exp-1))
    qx2[i] = MM_PER_M * np.abs(S2 - S_Soulis2[i]) * porosity * soilDepth/delt
    
figure,axs = plt.subplots(1,2, figsize= (20,20))

axs[0].plot(time_sim, subflw1, '#ff0000', label='Stand-alone WATDRN')
axs[0].text(14.8, 1.05, '$X_{L}$ = 50 m')
axs[0].legend(fontsize = 14, loc = 'upper left',frameon=False)

axs[0].plot(time_sim, qx1, '#00007fe6', label='Explicit WATDRN')
axs[0].legend(fontsize = 14, loc = 'upper left',frameon=False)

axs[1].plot(time_sim, subflw2, '#ff0000', label='Stand-alone WATDRN')
axs[1].text(14.8, 1.05, '$X_{L}$ = 10 m')
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)

axs[1].plot(time_sim, qx2, '#00007fe6', label='Explicit WATDRN')
axs[1].legend(fontsize = 14, loc = 'upper left',frameon=False)

# set axes and title 
axs[0].set_title('WATDRN simulaiton')
axs[1].set_title('WATDRN simulaiton')

axs[0].set_xlabel('Time (days)')
axs[1].set_xlabel('Time (days)')

axs[0].set_ylabel('Hillslope outflow (mm $h^-1$)')
axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')

axs[0].set_xlim(0,duration)
axs[0].set_ylim(0,2.5)
axs[0].grid(alpha=0.5)

axs[1].set_xlim(0,duration)
axs[1].set_ylim(0,2.5)
axs[1].grid(alpha=0.5)

plt.savefig(outdir+'WATDRN_SA_explicit.png', format='png', dpi=300)
plt.close()    
        
#%% simulate WATDRN based on different xdrainh values 
# the purpose is to find out how the outflow is responds based on changing values of xdrainh

# initialize variables  

# produce  a range for xdrainh 
xdrainh = np.linspace(0.05,1.05,11)

k = len(xdrainh)

asat_t1  = np.zeros((p,k)) 
subflw1  = np.zeros((p,k)) 
basflw1  = np.zeros((p,k))
satfc1   = np.zeros((p,k))

asat_t2  = np.zeros((p,k)) 
subflw2  = np.zeros((p,k)) 
basflw2  = np.zeros((p,k))
satfc2   = np.zeros((p,k))

for j in range(k):
    
    lambdda        = -np.log(xdrainh[j])/h0
    ktop           = ksat * np.exp(lambdda*ztop)      #  [m hour^-1]
    kl             = ktop * np.exp(lambdda*delzw)       #  [m hour^-1]
    
    # calculate grkeff for two hillslope lengths 
    grkeff1         = kl*xslope*2.0*Dd1/(1+xslope**2)  # [hour^-1]
    grkeff2         = kl*xslope*2.0*Dd2/(1+xslope**2)
    
    # loop over entire time window of simulation for two hillslope lengths 
    for i in range(p):
        # hillslope with length of totalLength1
        [asat_t1[i,j], subflw1[i,j], basflw1[i,j], satfc1[i,j]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff1,S_Soulis1[i], IWF,\
                                                                          ILG,IL1,IL2,1,delt)
        
        # NB: in the source code of WATDRN the lateral flow is not devided by time step
        subflw1[i,j] = MM_PER_M * subflw1[i,j]/delt                # convert to [mm hour^-1]
    
        
        # hillslope with length of totalLength2    
        [asat_t2[i,j], subflw2[i,j],basflw2[i,j],satfc2[i,j]] = WATDRN (delzw,bij,thpor_avail,ksat,grkeff2,S_Soulis2[i],IWF,\
                                                                ILG,IL1,IL2,1,delt)
        
        subflw2[i,j] = MM_PER_M * subflw2[i,j]/delt               # convert to [mm hour^-1]


# plot results 
figure,axs = plt.subplots(1,2, figsize= (20,20))

for j in range(k):
    lambdda        = -np.log(xdrainh[j])/h0
    axs[0].plot(time_sim, subflw1[:,j], label='XDRAINH = %f'%xdrainh[j] + ' or (-lambda) = %f'%lambdda)
    axs[0].text(1.5, 65, '$X_{L}$ = 50 m')
    axs[0].legend(fontsize = 10, loc = 'upper right',frameon=False)
    
    axs[1].plot(time_sim, subflw2[:,j], label='XDRAINH = %f'%xdrainh[j] + ' or (-lambda) = %f'%lambdda)
    axs[1].text(1.5, 65, '$X_{L}$ = 10 m')
    axs[1].legend(fontsize = 10, loc = 'upper right',frameon=False)
    
    # set axes
    axs[0].set_xlabel('Time (days)')
    axs[1].set_xlabel('Time (days)')
    
    axs[0].set_ylabel('Hillslope outflow (mm $h^-1$)')
    axs[1].set_ylabel('Hillslope outflow (mm $h^-1$)')
    
    axs[0].set_xlim(0,duration)
    axs[0].set_ylim(0,70)
    axs[0].grid(alpha=0.5)
    
    axs[1].set_xlim(0,duration)
    axs[1].set_ylim(0,70)
    axs[1].grid(alpha=0.5)

plt.savefig(outdir+'WATDRN_synthetic_xdrainh.png', format='png', dpi=300) 
plt.close()
   