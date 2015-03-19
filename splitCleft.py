
import numpy as np
import matplotlib.pylab as plt 

##
## PLOTTING etc 
## 


def plotter(model,tsteps,states,monitored):
    plt.subplot(1,3,1)
    RI =1- np.sum(states[:,[idxR,idxO,idxI]],axis=1)
    plt.plot(tsteps,states[:,idxR], label="R")
    plt.plot(tsteps,states[:,idxO], label="O")
    plt.plot(tsteps,states[:,idxI], label="I")
    plt.plot(tsteps,RI, label="RI")
    plt.legend(loc=0)

    plt.subplot(1,3,2)
    plt.plot(tsteps,states[:,idxO], label="O")
    plt.legend(loc=0)

    plt.subplot(1,3,3)
    plt.plot(tsteps[1:],monitored[:,0])
    
def validatePlots(refModel,refStates,splitModel,splitStates,tsteps,slabel="Cai",plotRoot=None):
  # state plots   
  plt.figure()  
  idxCai=refModel.state_indices(slabel)
  plt.plot(tsteps,refStates[:,idxCai],label="single")
  idxCai=splitModel.state_indices(slabel)
  plt.plot(tsteps,splitStates[:,idxCai],'g.',label="split")
  plt.legend()    
  plt.title("Total %s"%slabel)  
  if plotRoot!=None:   
    fileName="%s_%s.png"%(plotRoot,slabel)    
    plt.gcf().savefig(fileName,dpi=300)  
  
  # junctional Ca  
  plt.figure()
  idxCai=refModel.state_indices("Ca_jct1")
  plt.plot(tsteps,refStates[:,idxCai],label="single/jct")
  idxCai=splitModel.state_indices("Ca_jct1")
  #print idxCai
  plt.plot(tsteps,splitStates[:,idxCai],'g+',label="split/jct1")
  idxCai=splitModel.state_indices("Ca_jct2")
  #print idxCai  
  plt.plot(tsteps,splitStates[:,idxCai],'r.',label="split/jct2")
  plt.legend()    
  plt.title("Junctional Cai") 
  if plotRoot!=None:   
    fileName="%s_junctionalCas.png"%(plotRoot)    
    plt.gcf().savefig(fileName,dpi=300)  
    
    
def ValidateFluxes(ref,model,split,msplit,mlabel="J_Ca_SL_myo",mlabelAlt=None):
    if mlabelAlt==None:
      GetMonitored(ref,model,tsteps,mlabel=mlabel)
    else:    
      GetMonitored(ref,model,tsteps,mlabel=mlabelAlt)
    GetMonitored(split,msplit,tsteps,mlabel=mlabel)
    plt.plot(tsteps,ref.js,'b-')
    plt.plot(tsteps,split.js,'g.')
    plt.title("%s"%mlabel)
    
def GetMonitored(case,model,tsteps,mlabel="J_Ca_SL_myo"):
  case.js = np.zeros(np.shape(tsteps)[0])
  for i, ti in enumerate(tsteps):    
    r=model.monitor(case.states[i,:],ti,case.p)
    case.js[i]=r[model.monitor_indices(mlabel)]

##
## PICKLE
##
import cPickle as pickle

def WritePickle(states,params,fileName='data.pickle'):
  data1={'states':states,'params':params}
  output = open(fileName, 'wb')
  pickle.dump(data1, output)
  output.close()

def ReadPickle(fileName='data.pickle'):
  
  pkl_file = open(fileName, 'rb')
  data1 = pickle.load(pkl_file)
  states = data1['states']
  params = data1['params']  
  pkl_file.close()
  return states, params
    


##
## Models 
## 
class empty:pass


#t=0; dt=4; dtn=dt/2000.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
#t=0; dt=25; dtn=dt/13.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
from scipy.integrate import odeint

## Params 
# for baseline model (note jct1 instead of jct)
import shannon_2004_hack as model
paramList=["h","j","m","Xr","Xs","X_tos","Y_tos","R_tos","X_tof","Y_tof","d","f",
           "fCaB_SL","fCaB_jct1","R","O","I",
           "Ca_TroponinC","Ca_TroponinC_Ca_Mg",
           "Mg_TroponinC_Ca_Mg","Ca_Calmodulin","Ca_Myosin","Mg_Myosin","Ca_SRB",
           "Na_jct1_buf","Na_SL_buf","Na_jct1","Na_SL","Nai","Ca_Calsequestrin",
           "Ca_SLB_SL","Ca_SLB_jct1","Ca_SLHigh_SL","Ca_SLHigh_jct1","Ca_SR","Ca_jct1",
           "Ca_SL","Cai","V"]

# parameters split for 'split' model 
import shannon_splitcleft as msplit
changedList = [
"Fx_Na_jct1",
"Fx_Na_jct2",
"Fx_NaBk_jct1",
"Fx_NaBk_jct2",
"Fx_NaK_jct1",
"Fx_NaK_jct2",
"Fx_Ks_jct1",
"Fx_Ks_jct2",
"Fx_Cl_jct1",
"Fx_Cl_jct2",
"Fx_ICaL_jct1",
"Fx_ICaL_jct2",
"Fx_NCX_jct1",
"Fx_NCX_jct2",
"Fx_SLCaP_jct1",
"Fx_SLCaP_jct2",
"Fx_CaBk_jct1",
"Fx_CaBk_jct2",
"ks1",
"ks2",
"KSRleak1",
"KSRleak2",
"Bmax_SLB_jct1",
"Bmax_SLB_jct2",
"Bmax_SLHigh_jct1",
"Bmax_SLHigh_jct2",
"Bmax_jct1",
"Bmax_jct2"]

# Split model 
# Find parameters changed in split model and define lists specific to each junction 
s=msplit.init_state_values()
p=msplit.init_parameter_values()

paramsJct1only = np.zeros(np.shape(changedList)[0])
for i,statei in enumerate(changedList):
  paramsJct1only[i]=p[msplit.parameter_indices(statei)] 
  #print p[msplit.parameter_indices(statei)]   

paramsJct2only = paramsJct1only*0
paramsJct2only[1::2]=paramsJct1only[::2]
paramsJct2only[::2] =paramsJct1only[1::2]
#print paramsJct1only
#print paramsJct2only

# run parameters 
stim_period = 400
stim_start = 0 
t=0; dt=1000; dtn=dt/1000.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)


##
## Runs
##
refPickle = "ref.pkl"
refSplitPickle = "refsplit.pkl"

## Baseline 
def RunBaseline(pklName="ref.pkl"):
  ref = empty()
  s=model.init_state_values()
  p=model.init_parameter_values()
  p[model.parameter_indices("stim_period")] = stim_period
  p[model.parameter_indices("stim_start")] = stim_start
  ref.states = odeint(model.rhs,s,tsteps,(p,),mxstep=1000,hmax=.03,rtol=1e-12, atol=1e-12)

  WritePickle(ref.states,p,pklName)

  return ref




## Apply parameters for each junction 

def RunBaselineSplit(pklName=refSplitPickle):
  w1 = 0.5 # equal weighting between 1 and 2 
  msplit.w1 = w1
  wParams = paramsJct1only*w1 + paramsJct2only*(1-w1)
  jct12 = empty()

  s=msplit.init_state_values()
  p=msplit.init_parameter_values()
  for i,statei in enumerate(changedList):
    p[msplit.parameter_indices(statei)] = wParams[i]
    #print p[msplit.parameter_indices(statei)] 
    
  p[msplit.parameter_indices("stim_period")] = stim_period
  p[msplit.parameter_indices("stim_start")] = stim_start

  jct12.states = odeint(msplit.rhs,s,tsteps,(p,),mxstep=10000,hmax=.03,rtol=1e-12, atol=1e-12)
  jct12.p = p

  WritePickle(jct12.states,jct12.p,pklName)



# if run
doRun=0
if doRun:
  ref = RunBaseline(pklName=refPickle)
else:
  ref = empty()
  ref.states, ref.p = ReadPickle(refPickle)


doRun=0
if doRun:
  jct12 = RunBaselineSplit(pklName=refSplitPickle)
else:
  jct12 = empty()
  jct12.states, jct12.p = ReadPickle(refSplitPickle)

validatePlots(model,ref.states,msplit,jct12.states,tsteps,plotRoot="bothactive_splitXX")

