from separate import * # This imports separate fluxes
import shannon_2004 as model
#import shannon_pkh as model

## MONITORS
# New gotran 
#monitors = ["fCa_SL","fCa_jct","i_NaCa","j_rel_SR","j_pump_SR","i_Stim"]
#units = ["unk","unk","unk","unk","unk","unk"]
#idxMonitors = np.zeros(np.shape(monitors)[0],dtype=int)
#for i,monitor in enumerate(monitors):
#  idxMonitors[i] = model.monitor_indices(monitor)


## STATE VAR
Cai_idx =  model.state_indices( "Cai" )
Ca_SR_idx =  model.state_indices( "Ca_SR" )
V_idx =  model.state_indices( "V" )


## PARAMS
# stim_period=(121
#stim_period_pIdx=121
#V_max_Jpump_pIdx = 71 # SERCA
#V_max_pIdx = 45 # NCX

# old gotran 
stim_period_pIdx =  model.param_indices( "stim_period" )
V_max_pIdx =  model.param_indices( "V_max" )

# new gotran
#stim_period_pIdx =  model.parameter_indices( "stim_period" )
#V_max_pIdx =  model.parameter_indices( "V_max" )



## Misc
mM_to_uM = 1e3

## Monitors 
#huskeypm@huskeypm-ubuntu12:~/sources/wholecell$ grep monitor shannon_2004.ode 
#monitor(fCa_SL) 
#monitor(fCa_jct) 
#monitor(i_NaCa)   
#monitor(j_rel_SR)
#monitor(j_pump_SR)
#monitor(i_Stim)
fCa_SL  =0
fCa_jct =1
i_NaCa  =2 
j_rel_SR=3
j_pump_SR=4
i_CaL=5
i_Stim = 6
totMonitors = 7

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# get monitored values 
def monitorstepper(model,states,pi,tsteps):
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::] 
  jall = np.zeros((np.shape(tstepsm1)[0],totMonitors)) 

  jSums = np.zeros(totMonitors)
  for i,t in enumerate(tstepsm1):
    # get current state
    si = states[i,:]
    # extract monitored fluxes 
    jis = model.monitor(si, t, pi) 
    #print np.shape(jis)
    jall[i,] = jis
    
    # sum fluxes 
    jSums += jis*dtt   
    
    
  return (tstepsm1,jall)  

def init():
  # old gotran 
  s=model.init_values()
  p=model.default_parameters()

  # new gotran 
  #s=model.init_state_values()
  #p=model.init_parameter_values()
  t=0; #dt=1000; dtn=5;

  model.s = s; model.t =t; model.p = p


# run simulation 
# dtn <1, otherwise action potential isn't correct
def runner(dt=1000,dtn=1.,\
           stim_period=1000,\
  #         V_max_Jpump = 0.0053114, # SERCA, [mM/ms]
  #         V_max = 9, # NCX, [uA/uF]
           mxsteps=500):

  s = model.s; t =model.t; p = model.p
  
  # Important variables, states
  p[stim_period_pIdx]=stim_period
  #p[V_max_Jpump_pIdx]= V_max_Jpump 
  #p[V_max_pIdx] = V_max


  
  # Basic run and grab outcomes
  tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
  
  states = odeint(model.rhs,s,tsteps,(p,),mxstep=mxsteps)
  
  # get monitored variables
  (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)
  
  return (p,states,ts,js)


def plotting(p,states,ts,js,case="default"):
  tsteps = np.zeros(np.shape(ts)[0]+1)
  tsteps[1:]=ts
  state_indices = model.state_indices    


  ## Ca transients 
  plt.figure()
  plt.subplot(2,1,1)
  plt.plot(tsteps,mM_to_uM*states[:,Cai_idx],label=Cai_idx)
  plt.title("Cytosolic Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplot(2,1,2)
  plt.plot(tsteps,mM_to_uM*states[:,Ca_SR_idx],label=Ca_SR_idx)
  plt.title("SR Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplots_adjust(hspace=0.5) 
  plt.gcf().savefig(case+"_calcium.png",dpi=300)
  
  ## voltage 
  plt.figure()
  plt.plot(tsteps,states[:,V_idx],label=V_idx)
  plt.ylabel("V [mV]")
  plt.xlabel("t [ms]") 
  plt.gcf().savefig(case+"_potential.png",dpi=300)
  
  ## fluxes
  (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)

  #
  plt.figure(figsize=(10,5))
  plt.subplot(1,2,1)
  plt.plot(ts,js[:,fCa_SL_idx],label="fCa_SL_idx")
  plt.plot(ts,js[:,fCa_jct_idx],label="fCa_jct_idx")
  plt.legend(loc=0)
  
  plt.subplot(1,2,2)
  plt.plot(ts,js[:,i_NaCa_idx],label="i_NaCa_idx")
  plt.legend(loc=0)

  #
  plt.figure(figsize=(10,5))
  plt.subplot(1,2,1)
  plt.plot(ts,js[:,j_rel_SR_idx],label="j_rel_SR_idx")
  plt.legend(loc=0)
  
  plt.subplot(1,2,2)
  plt.plot(ts,js[:,j_pump_SR_idx],label="j_pump_SR_idx")
  plt.legend(loc=0)

  # w
  plt.figure(figsize=(10,5))
  plt.subplot(1,2,1)
  plt.plot(ts,js[:,i_CaL_idx],label="i_CaL")               
  plt.legend(loc=0)
  
  plt.subplot(1,2,2)
  plt.plot(ts,js[:,i_Stim_idx],label="i_Stim_idx")
  plt.legend(loc=0)

  #plt.gcf().savefig(case+"_fluxes.png",dpi=300)

  # returning only 
  return js[:,j_rel_SR_idx]

# <codecell>



#!/usr/bin/env python
import sys
#
# Revisions
#       10.08.10 inception
#

def doit():
  init()
  (p,si,tsi,jsi) = runner(dt=1000)
  plotting(p,si,tsi,jsi)
  raise RuntimeError("Not actually a validation yet")

def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-validation"):
      #arg1=sys.argv[i+1] 
      doit()
      quit()
  





  raise RuntimeError("Arguments not understood")




