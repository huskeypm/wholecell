
from separate import * # This imports separate fluxes

import shannon_2004 as model

## STATE VAR
#Cai=37, V=38
#Cai_idx=37; 
#Ca_SR_idx = 25
#V_idx=38
Cai_idx =  model.state_indices( "Cai" )
Ca_SR_idx =  model.state_indices( "Ca_SR" )
V_idx =  model.state_indices( "V" )


## PARAMS
# stim_period=(121
#stim_period_pIdx=121
#V_max_Jpump_pIdx = 71 # SERCA
#V_max_pIdx = 45 # NCX
stim_period_pIdx =  model.param_indices( "stim_period" )
V_max_Jpump_pIdx =  model.param_indices( "V_max_Jpump" )
V_max_pIdx =  model.param_indices( "V_max" )



## Misc
mM_to_uM = 1e3

## Monitors 
# WARNING: defined in  separate.py

import shannon_2004 as model
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# get monitored values 
def monitorstepper(model,states,pi,tsteps):
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::] 
  jall = np.zeros((np.shape(tstepsm1)[0],totMonitored)) 

  jSums = np.zeros(totMonitored)
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
  s=model.init_values()
  p=model.default_parameters()
  t=0; #dt=1000; dtn=5;

  model.s = s; model.t =t; model.p = p


# run simulation 
def runner(dt=1000,dtn=5,\
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
  plt.plot(tsteps,mM_to_uM*states[:,state_indices("Cai")],label=Cai_idx)
  plt.title("Cytosolic Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplot(2,1,2)
  plt.plot(tsteps,mM_to_uM*states[:,state_indices("Ca_SR")],label=Ca_SR_idx)
  plt.title("SR Ca")
  plt.ylabel("Ca [uM]")
  plt.xlabel("t [ms]") 
  
  plt.subplots_adjust(hspace=0.5) 
  plt.gcf().savefig(case+"_calcium.png",dpi=300)
  
  ## voltage 
  plt.figure()
  plt.plot(tsteps,states[:,state_indices("V")],label=V_idx)
  plt.ylabel("V [mV]")
  plt.xlabel("t [ms]") 
  plt.gcf().savefig(case+"_potential.png",dpi=300)
  
  ## fluxes
  (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)

  if 0:
    plt.figure(figsize=(10,10))
    plt.subplot(2,2,1)
    plt.plot(ts,js[:,fCa_SL],label="fCa_SL")
    plt.plot(ts,js[:,fCa_jct],label="fCa_jct")
    plt.legend(loc=0)
  
    plt.subplot(2,2,2)
    plt.plot(ts,js[:,i_NaCa],label="i_NaCa")
    plt.legend(loc=0)
  
  
    plt.subplot(2,2,3)
    plt.plot(ts,js[:,j_rel_SR],label="j_rel_SR")
    plt.legend(loc=0)
  
  
    plt.subplot(2,2,4)
    plt.plot(ts,js[:,j_pump_SR],label="j_pump_SR")
    plt.legend(loc=0)

  plt.gcf().savefig(case+"_fluxes.png",dpi=300)

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




