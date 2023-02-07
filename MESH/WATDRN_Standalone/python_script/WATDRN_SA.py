# -*- coding: utf-8 -*-
"""
NAME 
    WATDRN_SA
    
Input 

        saved CLASSW output variables and parameter 
Output     
    
    
PURPOSE    
        The purpose of this program is to run WATDRN as stand-alone 
        and by calling both WATROF and WATDRN subroutines. The input variables and                
        and parameters are obtained from MESH Fraser setup which was run for one month
    
PROGRAMMER: Ala Bahrami

REVISION HISTORY
    20230205 -- Initial version created and posted online 
    20230206 -- WATROF and WATDRN subroutines have been added based on fortran codes 
    
    
    
REFERENCE
               WATDRN_SA.f90, WATROF.f90, WATDRN.f90 
    
See also: 
"""
#%% import modules 
import xarray as xs
import numpy as np

#%% input files and local variables  
# files 
fpath_cs = 'THLQCS_2000_09_01.nc'
fpath_gs = 'THLQGS_2000_09_01.nc'
fpath_c  = 'THLQC_2000_09_01.nc'
fpath_g  = 'THLQG_2000_09_01.nc'
fpath_p  = 'parameter.nc'
	  
#variables dimensiosn 
THLQ_NAME  ="soil_moisture"
Level_NAME ="level"
Time_NAME  ="time"
		
#constant parameters
ILG = 15236 
IG = 4 
IL1 = 1 
IL2 = ILG
levelp = 30 
level = 34 
t_steps = 48 

#%% read input netcdf files  
data = xs.open_dataset(fpath_p)
data.close()
datp = data[THLQ_NAME].values 

data = xs.open_dataset(fpath_cs)
data.close()
datcs = data[THLQ_NAME].values 

#%% assign variables 
# NB: matrix order is inversed in  python compared to fortran
XSLOPE    = datp[0,0,:] 	   
XDRAINH   = datp[0,1,:] 
MANNING_N = datp[0,2,:]       
DD 	      = datp[0,3,:] 
KSAT      = datp[0,4,:]
DELZW 	  = datp[0,5:9,:].transpose()  
THPOR     = datp[0,9:13,:].transpose()
THLMIN    = datp[0,13:17,:].transpose()
BI 	      = datp[0,17:21,:].transpose()
ISAND     = datp[0,21:25,:].astype(int).transpose()
IWF       = datp[0,25,:].astype(int)
BULK_FC   = datp[0,26:30,:].transpose()  


#%% WATDRN 
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
    c      = np.zeros(ilg)            # Clapp-Hornberger connectivity index (c>1)
    cm1    = np.zeros(ilg)            # c-1
    c2m1   = np.zeros(ilg)            # 2*c-1
    asatc  = np.zeros(ilg)            # bulk saturation at the critical time tc
    asat00 = np.zeros(ilg)            # MIN(1,asat0) because we don't handle supersaturated soils
    tc     = np.zeros(ilg)            # critical time at which the seepage face becomes unsaturated
    ratiot = np.zeros(ilg)            # ratio tc/t (t0 or t1) if t>tc and t/tc if t<=tc
    satspf = np.ones(ilg, dtype=bool) # indicates if seepage face is saturated
                					  # equivalent to knowing if t<=tc 
    
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
    
        #c and c factors
        c[i]    = 2.*bcoef[i]+3.
        cm1[i]  = c[i]-1.
        c2m1[i] = 2.*c[i]-1.
        
        # bulk saturation at critical time
        #(just before the seepage face becomes unsaturated)
        asatc[i] = 1.-1./c[i]
    
        # layer average saturation asat0 may be greater than 1
        #e.g. frost heave but it is not possible for wat_drain
        asat00[i] = np.min((1., asat0[i]))
        #assess if seepage face is saturated at initial time
        satspf[i] = asat00[i] >= asatc[i]
        
    # **********************************************************************
    #      STEP 1: Find theoretical time t0 elapsed since element was last
    #              saturated and estimate baseflow rate at initial time
    #              Also estimate fraction of surface that is saturated
    # ********************************************************************** 
    tc = np.zeros(ilg)     # critical time at which the seepage face becomes unsaturated
    ratiot = np.zeros(ilg) # ratio tc/t (t0 or t1) if t>tc and t/tc if t<=tc
    basflw = np.zeros(ilg) # baseflow rate during the time step (m)
    satsfc = np.zeros(ilg) # saturated fraction of the surface (0 to 1)
    
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        # determine time at which seepage face becomes unsaturated
        tc[i] = thpora[i]/(c[i]*grkeff[i])
    
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        #find theoretical start of recession (t0) from bulk saturation
        # and at the same time estimate baseflow based on rate at t0
        if (satspf[i]):
          #saturated seepage face at initial time:
          # compute t0/tc
          ratiot[i] = c[i]*(1.-asat00[i])
          # normalized baseflow rate
          basflw[i] = 1.-c[i]*c[i]/c2m1[i]*(1.-asat00[i])
          #the fraction of the surface that is saturated at t0
          #varies linearly with t0/tc
          satsfc[i]= 1.-ratiot[i]
        else:
          #unsaturated seepage face at initial time:
          #calculate tc/t0 instead of t0 to avoid overflow 
          ratiot[i] = (asat00[i]/asatc[i])**cm1[i]
          #normalized baseflow rate
          basflw[i] = cm1[i]/c2m1[i]*ratiot[i]*asat00[i]/asatc[i] 
          # the fraction of the surface that is saturated at t0 is zero
          satsfc[i] = 0.
            
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        #Compute baseflow in m from normalized baseflow rate
        basflw[i] = grksat[i]*basflw[i]*delt
        
    # **********************************************************************
    #      STEP 2: Find theoretical time t1 at the end of the time step
    # **********************************************************************
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        if (satspf[i]):
            # Assess if seepage face will still be saturated at the
            #  end of the time step
            satspf[i] = tc[i] * ratiot[i] + delt <= tc[i]
            if (satspf[i]):
                # Seepage face still saturated, compute t1/tc from t0/tc
                ratiot[i]= (tc[i]*ratiot[i]+delt)/tc[i]
            else:
                # Seepage face not saturated anymore, compute tc/t1
                ratiot[i] = tc[i]/(tc[i]*ratiot[i]+delt)
        else:
            # If seepage face was not saturated initially, we compute
            # tc/t1=tc/(t0+delt) from tc/t0
            ratiot[i] = tc[i]*ratiot[i]/(tc[i]+delt*ratiot[i])
            
    # **********************************************************************
    #      STEP 3: Obtain bulk saturation at the end of the time step
    #              and interflow amount
    # **********************************************************************  
    asat1  = np.zeros(ilg)
    subflw = np.zeros(ilg) 
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        if (satspf[i]):
           # saturated seepage face at the end of the time step
           asat1[i] = 1.-ratiot[i]/c[i]
        else:
           # unsaturated seepage face at the end of the time step
           asat1[i]= asatc[i]*ratiot[i]**(1./cm1[i])
    
    for i in range(il1-1,il2):
        # cycle if not using watrof or latflow
        if (iwf[i] != wfk or asat0[i] == 0.0):
            continue
        # Sanity check: bulk saturation should not increase with time
        asat1[i] = np.min((asat00[i],asat1[i]))
        # Obtain interflow from the difference in bulk saturation
        subflw[i] = (asat00[i]-asat1[i])*thpora[i]*delzw[i]
    
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

#%% WATROF 
def WATROF(thliq,thice,zpond,tpond,ovrflw,tovrfl,      
           subflw,tsubfl,runoff,trunof,fi,zplim,        
           xslope,xdrainh,manning_n,dd,ksat,tbarw,     
           delzw,thpor,thlmin,bi,dodrn,dover,didrn,  
           isand,iwf,ig,ilg,il1,il2,bulk_fc):

    # input and output arrays 
    # THLIQ (ILG,IG)  
    # THICE (ILG,IG)
    # ZPOND (ILG)
    # TPOND (ILG)
    # OVRFLW(ILG)
    # TOVRFL(ILG)
    # SUBFLW(ILG,IG)
    # RUNOFF(ILG)
    # TRUNOF(ILG)
    
    
    # work arrays
    dodrn = np.zeros(ilg)
    dover = np.zeros(ilg)
    didrn = np.zeros((ilg, ig))
    vel_t0 = np.zeros(ilg)
    nuc_dover = np.zeros(ilg)
    
    delzwj = np.zeros(ilg)
    bij    = np.zeros(ilg)
    thporj = np.zeros(ilg)
    
    # constant values 
    tfrez = 273.16
    delt  = 1800.
    
    #return if no nml is expected to run in this cycle
    if (not(np.all(iwf == 1))):
        exit 
        
    # coefficients
    c1 = 2.0/3.0
    c2 = 1.5

    # parameter - will be used to compute xdrainh (the fractional change in horizontal
    # conductivity in a depth change h0) in Vincent's new formula.
    h0 = 1.0 

    #loop through each element
    for i in range(il1-1,il2):
        #skip if using flat class
        if (iwf[i] != 1):
            continue
         
        #compute overland flow and add to runoff and to the overall overland flow
        if ((fi[i] > 0.0) and (zpond[i] > zplim[i])):
            
            # calculate the depth of water available for overland flow
            dover[i] = zpond[i]-zplim[i]
            
            #calculate the flow velocity at the beginning of the timestep
            #(based on kinematic wave velocity) - eqn (1) in notes on overland flow
            vel_t0[i] = (dover[i]**c1) * np.sqrt(xslope[i])/(manning_n[i])
            
            #calculate a normalized unconstrained overland flow to avoid numerical
            #problems with a division of small dover(i) values.
            #eqn (29) in notes on overland flow
            nuc_dover[i] = -2 * dd[i] * vel_t0[i] * delt
            
            #constrained overland flow - limited by physically possible flow.
            #eqn (30) in notes on overland flow
            dodrn[i] = dover[i] * (1.0 - 1./((1.0-c1 * nuc_dover[i])**c2))

            # add overland flow to runoff and to the overall overland flow
            if(runoff[i] > 1.0e-08):
                trunof[i] = (trunof[i]*runoff[i]+(tpond[i]+tfrez)*  \
						    dodrn[i])/(runoff[i]+dodrn[i])
            
            runoff[i]    = runoff[i] + dodrn[i]
                        
            if(dodrn[i] > 1.0e-08):
                tovrfl[i] = (tovrfl[i]*ovrflw[i]+(tpond[i]+tfrez)*  \
						     fi[i]*dodrn[i])/(ovrflw[i]+fi[i]*dodrn[i])
                 
                ovrflw[i] = ovrflw[i] + fi[i] * dodrn[i]
                
                # subtract overland flow depth from the ponding depth
                zpond[i]  = zpond[i]- dodrn[i]
                             
    #-------------------------------------------------------------------------
    # compute interflow flow from each layer
    thliq_avail = np.zeros(ilg)
    thpor_avail = np.zeros(ilg)
    asat_t0     = np.zeros(ilg)
    ztop        = np.zeros((ilg, ig))
    
    grkeff      = np.zeros(ilg)
    
    # loop through each soil layer
    for j in range(ig):
        # loop through each element 
        for i in range(il1-1,il2):
            # skip if not using watrof
            if (iwf[i] != 1):
                continue
            #form vecotors for the layer - to be compatible with WATDRN arguments
            delzwj[i]   = delzw[i,j]
            bij[i]      = bi[i,j]
            thporj[i]   = thpor[i,j]
            
            # Find the top of each soil layer for the calculation of grkeff
            # modi here 
            if(j < (ig-1)):
                ztop[i,j+1] = ztop[i,j] - delzw[i,j]
            
            if (fi[i] > 0.0 and isand[i,j] >= -2 and \
                delzw[i,j] > 0.0):
                
                # determine available liquidwater in layer
                thliq_avail[i] = np.max((0.0,thliq[i,j]-thlmin[i,j]))
                
                # determine available porosity
                thpor_avail[i]    = np.max((thliq[i,j],thlmin[i,j], \
                                            thpor[i,j]-thice[i,j]))
							 
                # saturation defined as liquid water content over available porosity
                asat_t0[i]     = thliq_avail[i]/thpor_avail[i]
                
            # grkeff - average value of the parameter controlling the time scale of
            # interflow process - kl * (tile slope / tile length) (1/s)
            # Note: this formula is not the same as the one in Fhydro2_VF_20100226.f
            # and needs to be confirmed by Vincent
            
            # Integration of k across the layer -> kl
            xlambda        = -np.log(xdrainh[i])/h0
            ktop           = ksat[i]* np.exp(xlambda*ztop[i,j])
            kl             = ktop * exav(xlambda*delzw[i,j])
            grkeff[i]      = kl*xslope[i]*2.0*dd[i]/(1+xslope[i]**2)
            thpor_avail[i] = np.max((thlmin[i,j],thpor_avail[i]))
            
        # compute interflow from the layer (subflowj). Baseflow from the layer (basflwj) is 
        # also computed but is not used at present.
        [asat_t1, subflwj,basflwj,satfc] = WATDRN (delzwj,bij,thpor_avail,ksat,grkeff,asat_t0,iwf,\
                                                   ilg,il1,il2,1,delt)
        
        #loop through each element
        for i in range(il1-1,il2):
            #skip if not using watrof
            if (iwf[i] != 1):
                continue
            
            # allow lateral flow if liquid water content is greater than
            #bulk field capacity.
            if(thliq_avail[i] > 0.0 and thliq[i,j] >= bulk_fc[i,j]):
                didrn[i,j] = subflwj[i]
                
            #---------------------------------------------------------------------------
            #compute davail: volume of available water in a soil layer per land
            #                   area [m^3/m^2]
            
            davail = thliq_avail[i]*delzw[i,j]
            
            # limit the lateral flow not to exceed the available water in the layer
            # NB: take care about paranthesis here 
            didrn[i,j] = np.max((0.0,np.min((davail,didrn[i,j]))))
            
            # add the lateral flow to the runoff and to the subflow
            if(didrn[i,j] > 1.0e-8):
               trunof[i]   = (trunof[i]*runoff[i]+tbarw[i,j]* \
						    didrn[i,j])/(runoff[i]+didrn[i,j]) 
               runoff[i]   = runoff[i] + didrn[i,j]
               subflw[i,j] = subflw[i,j] + fi[i] * didrn[i,j]
               # remove the lateral flow from the layer
               thliq[i,j] = thliq[i,j]-didrn[i,j]/delzw[i,j]
        
    # NB: I dont see tsubfl where is called    
    return  thliq,thice,zpond,tpond,ovrflw,tovrfl,subflw,tsubfl,runoff,trunof,\
            dodrn,dover,didrn


#%% loop over time steps and call WATROF for each subareas 
for k in range(1):

    THLQCS = datcs[k,0:4,:].transpose()    
    THICCS = datcs[k,4:8,:].transpose() 
    ZPNDCS = datcs[k,8,:]     
    TPNDCS = datcs[k,9,:]
    OVRFLW = datcs[k,10,:]
    TOVRFL = datcs[k,11,:]
    SUBFLW = datcs[k,12:16,:].transpose()
    TSUBFL = datcs[k,16:20,:].transpose() 
    RUNFCS = datcs[k,20,:]
    TRNFCS = datcs[k,21,:]
    FCS    = datcs[k,22,:]
    ZPLMCS = datcs[k,23,:]
    TBRWCS = datcs[k,24:28,:].transpose()
    DODRN  = datcs[k,28,:]
    DOVER  = datcs[k,29,:]
    DIDRN  = datcs[k,30:34,:].transpose()    	

    # check to see whether the algorithm works     
    #print('Before calling WATROF')
    #print(SUBFLW[89,:])	

    # NB: I dont know whether I should intoduce DODRN, DOVER, DIDRN as inputs or 
    # as they are local variables I can put them as outputs only  	
    [THLQCS, THICCS, ZPNDCS, TPNDCS, OVRFLW, TOVRFL,\
      SUBFLW, TSUBFL, RUNFCS, TRNFCS,\
      DODRN, DOVER, DIDRN] =   WATROF(THLQCS, THICCS, ZPNDCS, TPNDCS, OVRFLW, TOVRFL, \
    	                                                         SUBFLW, TSUBFL, RUNFCS, TRNFCS, FCS, ZPLMCS,    \
    	                                                         XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBRWCS,   \
    	                                                         DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN,  \
    	                                                         ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
    #check to see whether the algorithm works
    #print('After calling WATROF')
    #print(SUBFLW[89,:])                                 