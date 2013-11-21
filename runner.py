## Separated CK/ATPase fluxes
## Params


Vmax_MM_f_idx=7
Vmax_MM_b_idx=6
Vmax_Mi_b_idx=14
Vmax_Mi_f_idx=15
ATPcyt_idx = 0
ADPcyt_idx = 1
Vcyt_idx = 29
Vims_idx = 30

import vanbeek_model_2007 as model
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
class empty:pass

## Functions 

# 'propagates' states forward by dt, similar to iterator.py
def propagator():
  initvals=model.init_values()
  params=model.default_parameters()
  t0=0; tss=10; dtn=0.01; tstepsSS = np.linspace(t0, tss, (tss-t0)/dtn+1)


  # disable CK activity
  vCK = 0
  params[Vmax_MM_f_idx] *= vCK
  params[Vmax_MM_b_idx] *= vCK
  params[Vmax_Mi_f_idx] *= vCK
  params[Vmax_Mi_b_idx] *= vCK


  ## get steady state first 
  statesSS= odeint(model.rhs,initvals,tstepsSS,(params,))
  finalStateSS = statesSS[-1::,]
  # plot steady state
  pSS, = plt.plot(tstepsSS,statesSS[:,ATPcyt_idx])

  ## big loop 
  TF = 50; dT = 5  # validated
  TF = 2; dT = 0.05 # validated
  noiseScale = 1.  # to simulate errors in PDE 

  # init w prev
  t0 = tss
  finalStateSS= np.ndarray.flatten(finalStateSS)
  statesPrev = finalStateSS

  # this is where the PDE 'concentration' values will be stored
  # I'm saving all variables here for simplicity, but really we'd 
  # only be concerned w ATP/ADP
  localVals = np.copy(finalStateSS)

  iters = TF/dT
  for i in np.arange(iters):
    print "iter ", i 
    ## get forward solution
    tf = t0 + dT
    tsteps = np.linspace(t0, tf, (tf-t0)/dtn+1)
    initvalsi = np.ndarray.flatten(statesPrev)
    statesi= odeint(model.rhs,initvalsi,tsteps,(params,))
    # plot next
    pEst,=plt.plot(tsteps,statesi[:,ATPcyt_idx],'k--')

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
    statesF[ATPcyt_idx] = localValsi[ATPcyt_idx] + noise[0]
    statesF[ADPcyt_idx] = localValsi[ADPcyt_idx] + noise[1]
    localVals = localValsi
    pC,=plt.plot(tf,localVals[ATPcyt_idx],'k.') 
    

    # update
    t0 = tf
    statesPrev = statesF
  

    #plt.show(); quit()


  plt.xlim([9.5,12])            
  #pSS.remove()
  pEst.remove()
  pC.remove()
  plt.title("ODE steady state + propagation of 'PDE' states by ODE fluxes")
  plt.legend([pSS,pEst,pC],["ODE Steady state","ODE propagation", "'Simulated' PDE update"],loc=1)
  plt.gcf().savefig("propagated.png") 

#  iterator mode shows that we can 'continue' a cellML simulation by initializing initvals with previous state values
def iterator():
  initvals=model.init_values()
  params=model.default_parameters()
  t0=0; tss=10; dtn=0.01; tstepsSS = np.linspace(t0, tss, (tss-t0)/dtn+1)


  # disable CK activity
  vCK = 0
  params[Vmax_MM_f_idx] *= vCK
  params[Vmax_MM_b_idx] *= vCK
  params[Vmax_Mi_f_idx] *= vCK
  params[Vmax_Mi_b_idx] *= vCK


  ## get steady state first 
  statesSS= odeint(model.rhs,initvals,tstepsSS,(params,)) 
  finalStateSS = statesSS[-1::,]
  # plot steady state
  plt.plot(tstepsSS,statesSS[:,ATPcyt_idx])

  ## big loop 
  TF = 50; dT = 5  # validated
  TF = 1; dT = 0.1 # validated

  # init w prev
  t0 = tss
  statesPrev = finalStateSS

  iters = TF/dT
  for i in np.arange(iters):  
    tf = t0 + dT
    tsteps = np.linspace(t0, tf, (tf-t0)/dtn+1)
    initvalsi = np.ndarray.flatten(statesPrev)    
    statesi= odeint(model.rhs,initvalsi,tsteps,(params,)) 

    # update
    statesPrev = statesi[-1::,]
    t0 = tf

    # plot next
    plt.plot(tsteps,statesi[:,ATPcyt_idx])


  plt.xlim([9,11])  
  plt.gcf().savefig("iterated.png") 

# validate our approach of evolving the ODE system over a time interval
# by fluxes extracted from the ode system solved at intermediate intervals 
def separateFluxesValidation():
  (tExact,statesExact,ts,ATP_stateflux) = separateFluxesTest(mode="stateflux")
  (tExact,statesExact,ts,ATP_sepflux) = separateFluxesTest(mode="nucleotideFlux")
  # jVol, jBound
  (tExact,statesExact,ts,ATP_volboundflux) = separateFluxesTest(mode="volboundFlux")
  (tExact,statesExact,ts,ATP_volboundCKATPaseflux) = separateFluxesTest(mode="volCKATPaseFluxes")


  plt.clf()
  plt.plot(tExact,statesExact[:,ATPcyt_idx],'k-',label="Exact")
  plt.plot(ts,ATP_stateflux,'r-.',label="by state flux")
  plt.plot(ts,ATP_sepflux,'b-.',label="by JAxP")
  plt.plot(ts,ATP_volboundflux,'g-.',label="by Jvol, Jbound")
  plt.plot(ts,ATP_volboundCKATPaseflux,'m-.',label="by JvolATPase,JvolCK, Jbound")
  plt.legend()
  file="separate.png"
  print "Generating ", file 
  plt.gcf().savefig(file)

# compare exact solution for ODE system, versus solution generated
# by flux estimates from ODE system solved on intermediate intervals
# mode: stateFlux - evolve nucleotide concentrations using [c(t1)-c(t0)/(t1-t0)]
# mode: nucleotideFlux - evolve nucleotide concentrations using jATP, jADP
# mode: volboundFlux - evolve nucleotide concentrations using volume (CKMM, ATP) and boundary (diff) fluxes
def separateFluxesTest(mode="stateflux"):
  s=model.init_values()
  pdef=model.default_parameters()
  # default params 
  pi = np.copy(pdef)
  results_d = empty()
  results_d.vCK = 1.0


  # Vcyte
  #print "Overriding Vcyt for now"
  #pi[Vcyt_idx]=1.
  #pi[Vims_idx]=1.
  Vcyt = pi[Vcyt_idx]
  Vims  = pi[Vims_idx]

  # we simulate the system to start from steady state 
  t=0.0; dt=5.; tsteps = np.linspace(t, t+dt, 10000)              
  ss= odeint(model.rhs,s,tsteps,(pi,))
  s = ss[-1,]

   # Test small range  
  (statesi,fluxes) = separateFluxes(model,s,pi,tsteps)
  jATP=fluxes.jATP;
  jADP=fluxes.jADP;
  dcdt_ATP=fluxes.dcdt_ATP;
  dcdt_ADP=fluxes.dcdt_ADP

  assert(np.abs(dcdt_ATP-jATP) < 1.0), "Poor match; choose smaller t interval"
  assert(np.abs(dcdt_ADP-jADP) < 1.0), "Poor match; choose smaller t interval"

  # get exact solution   
  T0 = dt
  localVals = np.copy(s)
  dT = 0.01
  TF=T0+3.0; tExact = np.linspace(T0, TF, (TF+dT)/dT)              
  statesExact = odeint(model.rhs,s,tExact,(pi,))
  #plot(tExact,statesExact[:,ATPcyt_idx])

  # do incrementally 
  statesPrev = s
  dT = 0.01  # time interval for each 'PDE' step 
  dtn = 0.001# time intervals for each ODE step
  t0 = T0
  iters = (TF-T0)/dT
  
  ATPs=[]
  ADPs=[]
  ts = []
  for i in np.arange(iters):
    #print "iter ", i
    ## get forward solution
    tf = t0 + dT
    tsteps = np.linspace(t0, tf, (tf-t0)/dtn+1)
    initvalsi = np.ndarray.flatten(statesPrev)

    ## just using concentration fluxes (from propagator code) 
    #statesi= odeint(model.rhs,initvalsi,tsteps,(pi,))
    # plot next
    #pEst,=plt.plot(tsteps,statesi[:,ATPcyt_idx],'k--')
    # get init/final states 
    #states0    = np.ndarray.flatten(statesi[0,]) # matched 'initvalsi'
    #statesF    = np.ndarray.flatten(statesi[-1::,])

    ## estimate change in local variables based on ode model  
    #plt.plot(tf,statesPrev[ATPcyt_idx],'.')                        
    #dStatedt = (statesF-states0) / dT
    #localValsi = localVals + dStatedt*dT

    (statesi,fluxes) = separateFluxes(model,initvalsi,pi,tsteps)
    jATP=fluxes.jATP;
    jADP=fluxes.jADP;
    dcdt_ATP=fluxes.dcdt_ATP;
    dcdt_ADP=fluxes.dcdt_ADP
    jVolATP=fluxes.jVolATP
    jVolATPaseATP=fluxes.jVolATPaseATP
    jVolCKATP=fluxes.jVolCKATP
    jBoundATP=fluxes.jBoundATP
    jVolADP=fluxes.jVolADP
    jVolATPaseADP=fluxes.jVolATPaseADP
    jVolCKADP=fluxes.jVolCKADP
    jBoundADP=fluxes.jBoundADP
    statesF    = np.ndarray.flatten(statesi[-1::,])     

   
    if(mode=="stateflux"):
      localValsi = localVals[[ATPcyt_idx,ADPcyt_idx]] + np.array([dcdt_ATP,dcdt_ADP])*dT     
    elif(mode=="nucleotideFlux"):
      localValsi = localVals[[ATPcyt_idx,ADPcyt_idx]] + np.array([jATP,jADP])*dT     
    elif(mode=="volboundFlux"):
      jATP=jVolATP + jBoundATP
      jADP=jVolADP + jBoundADP
      localValsi = localVals[[ATPcyt_idx,ADPcyt_idx]] + np.array([jATP,jADP])*dT     
    elif(mode=="volCKATPaseFluxes"): 
      jATP=jVolATPaseATP +jVolCKATP + jBoundATP
      jADP=jVolATPaseADP +jVolCKADP + jBoundADP
      localValsi = localVals[[ATPcyt_idx,ADPcyt_idx]] + np.array([jATP,jADP])*dT     

    else:
      raise RuntimeError("Not understood") 


    # when using PDE, the updated localVals might not exactly match the state values, so we update here
    statesF[ATPcyt_idx] = localValsi[ATPcyt_idx]
    statesF[ADPcyt_idx] = localValsi[ADPcyt_idx] 
    localVals = localValsi
    #pC,=plt.plot(tf,localVals[ATPcyt_idx],'k.')


    # update
    t0 = tf
    statesPrev = statesF
  
    # store 
    ATPs.append(statesF[ATPcyt_idx])
    ADPs.append(statesF[ADPcyt_idx])
    ts.append(tf)

  # plot 
  ATPs= np.asarray(ATPs)
  ts = np.asarray(ts)
  plt.plot(tExact,statesExact[:,ATPcyt_idx],'k-',label="Exact")
  plt.plot(ts,ATPs,'r-.',label="Js")
  plt.legend()
  plt.gcf().savefig(mode+".png")

  return(tExact,statesExact,ts,ATPs)


# This function computes the effective substrate flux 
# over a given time interval. 
# The results contain the concentration flux, the average flux due to the
# individual 'J' terms in the ode mode, as well as the averaged 'J' terms
# themselves
# not sure if this is the correct/efficient way of doing this
# runs iterator, grabs flux values at each time point 
def separateFluxes(model,s,pi,tsteps,warn=False):

  # run model to get states vs time  
  states = odeint(model.rhs,s,tsteps,(pi,))


  # all fluxes acting on ATP/ADPcyt
  # WARNING: these monitors are specific to a given ode file 
  # e.g. pay attention to the order of the monitor calls in ode
  #monitor(J_CKMi) 
  #monitor(J_syn)
  #monitor(J_diff_PCr) 
  #monitor(J_diff_Cr) 
  #monitor(J_diff_Pi) 
  JCKMM_idx=5 #monitor(J_CKMM)
  JdiffATP_idx = 6 #monitor(J_diff_ATP)
  JdiffADP_idx= 7 #monitor(J_diff_ADP)
  Jhyd_idx = 8 # monitor(J_hyd) 
  jSums= np.zeros(9)
  
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
    #print "####################"
    #ddprint "si",si[[0,5]]
    #print "m3",m3
    #print "m3",m3[6:9]

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
  dcdt = dc/dT

  #print "s0", states[0,ATP_idx]
  #print "JCKMM %f" % jAvg[JCKMM_idx]
  #print "Jhyd %f" % jAvg[Jhyd_idx]
  #print "JdiffADP %f" % jAvg[JdiffADP_idx]
  #print "JdiffATP %f" % jAvg[JdiffATP_idx]

  #dy[0] = (J_diff_ATP - J_CKMM - J_hyd)/V_cyt
  #dy[1] = (J_CKMM + J_diff_ADP + J_hyd)/V_cyt
  Vcyt = pi[Vcyt_idx]
  jATP=(jAvg[JdiffATP_idx] - jAvg[JCKMM_idx] - jAvg[Jhyd_idx])/Vcyt
  jADP=(jAvg[JdiffADP_idx] + jAvg[JCKMM_idx] + jAvg[Jhyd_idx])/Vcyt
  #print "ATP J: %f " % jATP
  #print "ADP J: %f " % jADP
  dcdt_ATP =dcdt[ATPcyt_idx]
  dcdt_ADP =dcdt[ADPcyt_idx]
  #print "dATP/dt: %f " % (dcdt_ATP)
  #print "dADP/dt: %f " % (dcdt_ADP)

  fluxes = empty()
  fluxes.jATP = jATP
  fluxes.jADP = jADP
  fluxes.dcdt_ATP = dcdt_ATP
  fluxes.dcdt_ADP = dcdt_ADP
  # jBound (boundary flux due to diffusion across sarcomere/IMS boundary) 
  fluxes.jBoundATP= jAvg[JdiffATP_idx]/Vcyt
  fluxes.jBoundADP = jAvg[JdiffADP_idx]/Vcyt
  # jBound (boundary flux due to diffusion across sarcomere/IMS boundary) 
  fluxes.jVolATP= ( - jAvg[JCKMM_idx] - jAvg[Jhyd_idx])/Vcyt
  fluxes.jVolADP= (   jAvg[JCKMM_idx] + jAvg[Jhyd_idx])/Vcyt

  # for CK/ATPase sepate 
  fluxes.jVolCKATP= ( - jAvg[JCKMM_idx])/Vcyt
  fluxes.jVolCKADP= (   jAvg[JCKMM_idx])/Vcyt
  fluxes.jVolATPaseATP= ( - jAvg[Jhyd_idx])/Vcyt
  fluxes.jVolATPaseADP= (   jAvg[Jhyd_idx])/Vcyt

  if(warn):
    print "Diff, dATP/dt vs Js: %f" % (\
      fluxes.dcdt_ATP-(fluxes.jBoundATP+fluxes.jVolATP))

  return (states,fluxes)

# validation of fig 2a in 
# van Beek, Cell Physiology, vol 293 no 3 
def validation():
  s=model.init_values()
  pdef=model.default_parameters()
  t=0; dt=10; dtn=0.001; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)

  # default params 
  pi = np.copy(pdef)
  states = odeint(model.rhs,s,tsteps,(pi,))

  jhyds = []
  Jhyd_idx = 8
  for i,t in enumerate(tsteps):
    # get current state
    si = states[i,:]
    # extract minitored fluxes 
    jis = model.monitor(si, t, pi)
    jhyds.append(jis[Jhyd_idx])

  jhyds = np.asarray(jhyds)
  plt.figure()
  plt.plot(tsteps,jhyds,'k-')
  plt.ylabel("j(ATP hydrolysis) [uM/s]") 
  plt.xlabel("t [s]") 
  plt.gcf().savefig("fig2vb.png") 
  
  print "WARNING: not a unit test - need to compare w van Beek fig"  



# validation of implementation of ODE model
# compare against fig10, row 1 of 
# van Beek, Cell Physiology, vol 293 no 3 
def validation2(): 
  s=model.init_values()
  pdef=model.default_parameters()
  t=0; dt=10; dtn=0.01; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)

  # default params 
  pi = np.copy(pdef)
  results_d = empty()
  results_d.vCK = 1.0
  results_d.states = odeint(model.rhs,s,tsteps,(pi,))
  
  # no ck
  pi= np.copy(pdef)
  vCK = 1e-9         
  pi[Vmax_MM_f_idx] *= vCK
  pi[Vmax_MM_b_idx] *= vCK
  pi[Vmax_Mi_f_idx] *= vCK
  pi[Vmax_Mi_b_idx] *= vCK
  results_nock = empty()
  results_nock.vCK = vCK
  results_nock.states= odeint(model.rhs,s,tsteps,(pi,))
  
  
  # low ck
  pi= np.copy(pdef)
  vCK = 0.02
  pi[Vmax_MM_f_idx] *= vCK
  pi[Vmax_MM_b_idx] *= vCK
  pi[Vmax_Mi_f_idx] *= vCK
  pi[Vmax_Mi_b_idx] *= vCK
  results_lowck = empty()
  results_lowck.vCK = vCK
  results_lowck.states= odeint(model.rhs,s,tsteps,(pi,))
  
  # 3x
  pi= np.copy(pdef)
  vCK = 3.0       
  pi[Vmax_MM_f_idx] *= vCK
  pi[Vmax_MM_b_idx] *= vCK
  pi[Vmax_Mi_f_idx] *= vCK
  pi[Vmax_Mi_b_idx] *= vCK
  results_3x= empty()
  results_3x.vCK = vCK
  results_3x.states= odeint(model.rhs,s,tsteps,(pi,))
  
  # 1/2
  pi= np.copy(pdef)
  vCK = 1e-9
  pi[Vmax_MM_f_idx] *= vCK
  pi[Vmax_MM_b_idx] *= vCK
  results_MiCKonly= empty()
  results_MiCKonly.vCK = vCK
  results_MiCKonly.states= odeint(model.rhs,s,tsteps,(pi,))
  
  # plot 
  plt.figure()
  plt.title("van beek: ADPcyt vs time") 
  ##plt.plot(tsteps,results_3x.states[:,ADPcyt_idx],'k.', label="vCK=%f"%results_3x.vCK)
  plt.plot(tsteps,results_d.states[:,ADPcyt_idx],'k-',linewidth=2.0,label="vCK=%f"%results_d.vCK)
  plt.plot(tsteps,results_lowck.states[:,ADPcyt_idx],'k.-',label="vCK=%f"%results_lowck.vCK)
  ##plt.plot(tsteps,results_nock.states[:,ADPcyt_idx],'k--', label="vCK=%f"%results_nock.vCK)
  ##plt.plot(tsteps,results_MiCKonly.states[:,ADPcyt_idx],'r--', label="MiCK only")
  plt.legend()
  plt.ylim([0,200])
  plt.xlim([9.2,10])
  plt.ylabel("[ADP] [uM]")
  plt.xlabel("time [s]")
  plt.gcf().savefig("ADP_vb.png") 
  
  plt.plot(tsteps,results_d.states[:,ADPcyt_idx],'k-',linewidth=2.0,label="vCK=%f"%results_d.vCK)
  
  
  #print "Plots of J_syns do not correspond to times (have emnail to Johan) " 
  #J_syns= np.asarray(model.ns["J_syns"])
  #plot(J_syns[:,])
  
  print "WARNING: not a unit test - need to compare w van Beek fig"  
  


import sys
#
# Revisions
#       10.08.10 inception
#

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Script for testing stuff with van beek model 
 
Usage:
"""
  msg+="  %s -validation/-iterator" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-validation"):
      validation()
      validation2()
      separateFluxesValidation()
      quit()
    if(arg=="-propagator"):
      propagator()
      quit()
    if(arg=="-test"):
      #separateFluxesTest(mode="volboundFlux")
      separateFluxesValidation()
      quit()

