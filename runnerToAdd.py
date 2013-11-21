isubstrate_idx=7

import shannon_2004 as model
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
class empty:pass

Cai_idx = 37
getFluxes=True

#def propagator():
if 0:

  initvals=model.init_values()    
  params=model.default_parameters()
  t0=0; tss=1000; dtn=5; tstepsSS = np.linspace(t0, tss, (tss-t0)/dtn+1)



  ## get steady state first 
  statesSS= odeint(model.rhs,initvals,tstepsSS,(params,))
  finalStateSS = statesSS[-1::,]
  # plot steady state
  pSS, = plt.plot(tstepsSS,statesSS[:,Cai_idx])

  # init w prev
  t0 = tss
  finalStateSS= np.ndarray.flatten(finalStateSS)
  statesPrev = finalStateSS
    
 
# bigerator     
if 1:  
  ## overwriting
  statesPrev = initvals
  finalStateSS= np.ndarray.flatten(statesSS[0,])
  t0=0  
    
  # this is where the PDE 'concentration' values will be stored
  # I'm saving all variables here for simplicity, but really we'd 
  # only be concerned w ATP/ADP
  localVals = np.copy(finalStateSS)    
    
  ## big loop 
  TF = 1000; dT = 10  # validated
  noiseScale = 0.  # to simulate errors in PDE 
    
  jAvgs=[]
  tAvgs=[]
  iters = TF/dT
  for i in np.arange(iters):
    #print "iter ", i
    ## get forward solution
    tf = t0 + dT
    tsteps = np.linspace(t0, tf, (tf-t0)/dtn+1)
    initvalsi = np.ndarray.flatten(statesPrev)
    statesi= odeint(model.rhs,initvalsi,tsteps,(params,))
    # plot next
    pEst,=plt.plot(tsteps,statesi[:,Cai_idx],'k--')

    # get flux info 
    pi = np.copy(params)  
    jAvg = separateFluxes(model,statesF,pi,tsteps)
    jAvgs.append(jAvg)
    tAvgs.append(t0+dT/2.)
    
    
    # get init/final states 
    states0    = np.ndarray.flatten(statesi[0,]) # matched 'initvalsi'
    statesF    = np.ndarray.flatten(statesi[-1::,])

    ## estimate change in local variables based on ode model  
    #plt.plot(tf,statesPrev[ATPcyt_idx],'.')                        
    dStatedt = (statesF-states0) / dT
    localValsi = localVals + dStatedt*dT
    # when using PDE, the updated localVals might not exactly match the state values, so we update here
    noise = noiseScale*np.random.randn(2)  # basically to add some uncertainity to mimic PDE errors 
    noise = 0.
    # NOT WORKING statesF[Cai_idx] = localValsi[Cai_idx] + noise[0]
    localVals = localValsi
    pC,=plt.plot(tf,localVals[Cai_idx],'k.')


    # update
    t0 = tf
    statesPrev = statesF   
    
  # return (jAvgs,tAvgs)  

def separateFluxes(model,s,pi,tsteps,warn=False):
  #specific to model
  jidx =0
  jSums= np.zeros(1)
  
  # run model to get states vs time  
  states = odeint(model.rhs,s,tsteps,(pi,))
    
  dtt= tsteps[1]-tsteps[0]
  tstepsm1 = tsteps[1::] 
    
  for i,t in enumerate(tstepsm1):
    # get current state
    si = states[i,:]
    # extract monitored fluxes 
    jis = model.monitor(si, t, pi)    
    
    
    # sum fluxes 
    jSums += jis*dtt    
    
  ## get average jAvg over the entire time interval 
  dT =tsteps[-1] - tsteps[0]
  jAvg = jSums/dT

  ## determine concentration flux from first/last state  
  f = states[-1,]
  s = states[0,]
  #dc= f[ATP_idx] - states[0,ATP_idx]
  dc= f - s
  dcdt = dc/dT    
    
    
  return jAvg  

# test
s=model.init_values()
pdef=model.default_parameters()
# default params 
pi = np.copy(pdef)  

jAvg = separateFluxes(model,s,pi,tsteps)
