from scipy.integrate import odeint
import numpy as np      


#monitors = ["fCa_SL","fCa_jct","i_NaCa","j_rel_SR","j_pump_SR","i_Stim"]
#units = ["unk","unk","unk","unk","unk","unk"]
def RunnerMonitor(model,tsteps,s,p,monitors):
    # run 
    states = odeint(model.rhs,s,tsteps,(p,),mxstep=5000)

    idxMonitors = np.zeros(np.shape(monitors)[0],dtype=int)
    for i,monitor in enumerate(monitors):
      idxMonitors[i] = model.monitor_indices(monitor)

    # pull out fluxes 
    dtt= tsteps[1]-tsteps[0]
    tstepsm1 = tsteps[1::]

    monitored = np.zeros([np.shape(tstepsm1)[0],np.shape(idxMonitors)[0]])
    for i,t in enumerate(tstepsm1):
        # get current state
        si = states[i,:]
        # extract monitored fluxes 
        allJs = model.monitor(si, t, p)
        monitored[i,:] = allJs[idxMonitors]
    return states,monitored

    
def GetMonitored(case,model,tsteps,mlabel="J_Ca_SL_myo"):
  case.js = np.zeros(np.shape(tsteps)[0])
  for i, ti in enumerate(tsteps):    
    r=model.monitor(case.states[i,:],ti,case.p)
    case.js[i]=r[model.monitor_indices(mlabel)]


##
##
##
import gotranJIT
from gotran.model.loadmodel import load_ode
def PrintStates(odeName="shannon_2004.ode"):
  ode = load_ode(odeName)
  
  for state in ode.full_states:
    print state.name,state.value

def PrintParams(odeName="shannon_2004.ode"):
  ode = load_ode(odeName)

  for param in ode.parameters:    
    print param.name,param.value


 

    



#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#
def doit(fileIn):
  1


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -print <states/params>" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-print"):
      arg1=sys.argv[i+1] 
      if arg1=="states":
        PrintStates()
      elif arg1=="params":
        PrintParams()
      else:
        raise RuntimeError("unknown") 
      quit()
  





  raise RuntimeError("Arguments not understood")




