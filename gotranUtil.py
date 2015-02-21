from scipy.integrate import odeint
import numpy as np      





# Runner

# In[ ]:


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


