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

    



