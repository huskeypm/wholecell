
# This function computes the effective substrate flux 
# over a given time interval. 
# The results contain the concentration flux, the average flux due to the
# individual 'J' terms in the ode mode, as well as the averaged 'J' terms
# themselves
# not sure if this is the correct/efficient way of doing this
# runs iterator, grabs flux values at each time point 

from scipy.integrate import odeint
import numpy as np
class empty:pass 

print "WARNING: ode idx is hardcoded" 
print "WARNING: This model is almost certinaly wrong (dbl check fluxes)" 
Cai_ode_idx = 37 


## THIS IS A TEMPLATE 
def separateFluxes(model,s,pi,tsteps,warn=False):

  # run model to get states vs time  
  states = odeint(model.rhs,s,tsteps,(pi,))
  dt = tsteps[-1] - tsteps[0] # TODO check 

  # namely, need to determine which fluxes are going into the PDE domain 
  fluxes = empty()

  ## TEMPLATE 
  #$ grep monitor shannon_2004.ode | perl -ane '$n++; chomp $_; ($nm)=($_=~/\((.*)\)/);printf("${nm}_idx=%d # $_\n",$n-1)'
  Vol_SR_idx=0 # monitor(Vol_SR)  
  Vol_myo_idx=1 # monitor(Vol_myo) 
  fCa_SL_idx=2 # monitor(fCa_SL) 
  fCa_jct_idx=3 # monitor(fCa_jct) 
  i_NaCa_idx=4 # monitor(i_NaCa)   
  j_rel_SR_idx=5 # monitor(j_rel_SR)
  j_pump_SR_idx=6 # monitor(j_pump_SR)
  J_Ca_SL_myo_idx=7 # monitor(J_Ca_SL_myo)
  dCa_TroponinC_idx=8 # monitor(dCa_TroponinC)
  i_Stim_idx=9 # monitor(i_Stim)   
  
  jSums= np.zeros(10) # total num of monitored fluxes 

 
  ## DONT TOUCH 
  # TODO put in its own function 
  ## grab js from each time step 
  # Idea here is that 
  # [c(Tf) - c(T0)]/(Tf-T0) = \int_T0^Tf j
  # --> delC ~ delT Sum[ j*dt] 
  # get dt between time steps
  # ignore first time step 
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::]
  for i,t in enumerate(tstepsm1):
    # get current state
    si = states[i,:]
    # extract minitored fluxes 
    jis = model.monitor(si, t, pi)

    ## sum fluxes 
    jSums += jis*dtt
    #print "m3/v",m3[7:9]/V_cyt

  ## get average jAvg over the entire time interval 
  dT =tsteps[-1] - tsteps[0]
  jAvg = jSums/dT

  ## determine concentration flux from first/last state  
  f = states[-1,]
  s = states[0,]
  #dc= f[ATP_idx] - states[0,ATP_idx]
  dc= f - s
  fluxes.dcdt = dc/dT


  ## END LEAVE MEE 


  ## TEMPLATE
  # from ODE: 
  # dCai_dt = -j_pump_SR*Vol_SR/Vol_myo + J_Ca_SL_myo/Vol_myo - one*dCa_cytosol_tot_bound
  fluxes.jBoundary = np.zeros(np.shape(s)[0]) 
  Vol_myo = jis[Vol_myo_idx] # unfortunate abuse of variable names  
  Vol_SR  = jis[Vol_SR_idx] # unfortunate abuse of variable names  
  fluxes.jBoundary[Cai_ode_idx] = jAvg[J_Ca_SL_myo_idx]/Vol_myo

  fluxes.jVol = np.zeros(np.shape(s)[0]) 
  # this is quirkly but I think we need to substract SR component, but 
  # add in TnC component (so that we are left with just the contribution 
  # from the SL (see SB expression for dCai) 
  fluxes.jVol[Cai_ode_idx] = -jAvg[j_pump_SR_idx]*Vol_SR/Vol_myo
  fluxes.jVol[Cai_ode_idx] += 0 # TODO add TnC component here

  #print Vol_myo
  #print jAvg[J_Ca_SL_myo_idx]
  #print jAvg[j_pump_SR_idx]
   
  
  return (states,fluxes)



