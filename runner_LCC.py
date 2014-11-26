## STATE VAR
## Misc
mM_to_uM = 1e3

## Monitors 
#huskeypm@huskeypm-ubuntu12:~/sources/wholecell$ grep monitor shannon_2004.ode 
i_CaL = 0  
totMonitors = 1

import shannon_LCC as model
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
    #print np.shape(jall)
    #print np.shape(jis)  
    
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
  V = -78,
  #         V_max_Jpump = 0.0053114, # SERCA, [mM/ms]
  #         V_max = 9, # NCX, [uA/uF]
  paramName=None,
  paramVal=0.,
           mxsteps=500):

  s = model.s; t =model.t; p = model.p
  
  # Important variables, states
  #p[stim_period_pIdx]=stim_period
  #p[V_max_Jpump_pIdx]= V_max_Jpump 
  param_indices = model.param_indices    

  p[param_indices("V")] = V
  if paramName!=None:
    p[param_indices(paramName)] = paramVal

  
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
  label="d"
  plt.plot(tsteps,mM_to_uM*states[:,state_indices(label)],label=label)       
  label="f"
  plt.plot(tsteps,mM_to_uM*states[:,state_indices(label)],label=label)       
  plt.legend(loc=0)
  plt.xlabel("t [ms]") 
  
  plt.subplots_adjust(hspace=0.5) 
  #plt.gcf().savefig(case+"_calcium.png",dpi=300)
  
  ## fluxes
  (ts,js)=monitorstepper(model,states,np.copy(p),tsteps)

  if 1:
    plt.figure(figsize=(10,10))
    plt.subplot(2,2,1)
    plt.plot(ts,js[:,i_CaL],label="I_CaL")  
    plt.legend(loc=0)
  
  #plt.gcf().savefig(case+"_fluxes.png",dpi=300)

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




