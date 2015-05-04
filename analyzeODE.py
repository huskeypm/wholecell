"""
For processing ODE outputs in support of Satin/Despa collaborations 
"""
import cPickle as pickle
import runner 
import matplotlib.pylab as plt
import numpy as np 
from runShannonTest import *
import taufitting as tf
runner.init()

mM_to_uM=1e3

def readPickle(name = "PCa0.75kss0.25.pickle"):
  print "Reading " + name  
  pkl_file = open(name, 'rb')
  data1 = pickle.load(pkl_file)  
  pkl_file.close()
     
  return data1  
    

def analyOut(data1,state="Cai",label=""):
  si = data1['s']
  pi = data1['p']

  np.shape(si)
  temp = si[:,runner.model.state_indices(state)]
  # take last third
  idx = np.int(np.shape(temp)[0]*2/3)
  #print idx
  ts = data1['t']
  ts = ts[idx:]
  temp = temp[(1+idx):]
  plt.plot(ts,temp*mM_to_uM,label=label)
    
  PCa = pi[runner.model.parameter_indices("PCa")]
  ks = pi[runner.model.parameter_indices("ks")]
  minCai = np.min(temp)
  maxCai = np.max(temp)  
    
  return (PCa,ks,minCai,maxCai)
    
#name =  "/tmp/PCa1.00ks1.00.pickle"      
#d = readPickle(name)
#pca,ks,minCai,maxCai = analyOut(d)

def ProcessOneDOutputs(var1Name,names,allVars,state="Cai",xlim=None,ylim=None,offsetMin=False):
  print "WARNING: does not include time steps" 
  for i,name in enumerate(names):               
      print name
      d = readPickle(name+".pickle")
      print np.shape(d['s'])
      s = d['s']
      si = s[:,runner.model.state_indices(state)]         
     
      # recenter each to minimum
      if offsetMin:
        # assume last third is in steady state
        inds  = np.shape(si)[0]
        inds = np.int(inds/3.)
        si = si - np.min(si[-inds:])
          
      plt.plot(si*mM_to_uM,label="%s=%4.2f"%(var1Name,allVars[i]))  

  if offsetMin:
    plt.title("Minima offset to 0. uM")

  plt.legend(loc=2)

  if xlim!=None:
    plt.xlim(xlim)
  if ylim!=None:
    plt.ylim(ylim)



  plt.ylabel("[%s] [uM]" % state)  
  plt.xlabel("timesteps []") # t [ms]") 
  name = state+"transients%s"%(var1Name)
  plt.gcf().savefig(name+".png")

def ProcessTwoDOutputs(allKeys,allVars,state="Cai",ylims=None,stim_period=1000,nameTag=None):
  vars1 = allVars[0]
  var1Name = allKeys[0]
  vars2 = allVars[1]
  var2Name = allKeys[1]

  outsMin = np.zeros([np.shape(vars1)[0],np.shape(vars2)[0]])
  outsMax = np.zeros(outsMin.shape)


  # assuming there exist two iterated var 
  #for i, PCa in enumerate(PCas):
  for i, var1Val in enumerate(vars1):
    plt.figure()  
    plt.title("%s=%3.2f"%(var1Name,var1Val))
    #for j, ks in enumerate(kss):
    for j, var2Val in enumerate(vars2):
        name =namer(var1Name,var1Val,var2Name,var2Val,stim_period=stim_period,tag=nameTag)+".pickle"
        #print name
        try: 
          d = readPickle(name) 
        except: 
          print name + " was not found. Skipping" 
          continue 
        #print np.shape(d['s'])
        dummy,dummy,minCai,maxCai = analyOut(d,state=state,label="%s=%3.2f"%(var2Name,var2Val))
        outsMin[i,j]= minCai 
        outsMax[i,j]= maxCai 
    if ylims!=None:    
      plt.ylim(ylims)  
    plt.ylabel(state+" [uM]")
    plt.xlabel("t [ms]")  
    plt.legend(loc=0,ncol=3)  
    name = state+"transients%s%3.2f"%(var1Name,var1Val)
    plt.gcf().savefig(name.replace(".","p")+".png")
    
  return outsMin,outsMax
    
# labeli - label corresponding to key[i]
def TwoDPlots(allKeys,allVars,outsMin, outsMax,label0="",label1="",state="Cai"):
    xv, yv = np.meshgrid(allVars[1],allVars[0])

    plt.figure()
    plt.subplot(2,2,1)
    plt.title(state+" min")
    plt.pcolormesh(xv,yv,outsMin*mM_to_uM)#,shading='gouraud')
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.ylabel(label0)
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.title(state+" max")
    plt.pcolormesh(xv,yv,outsMax*mM_to_uM)
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.ylabel(label0)
    plt.colorbar()
    plt.xlabel(label1)        
    
    
    plt.subplot(2,2,3)
    plt.title(state+" diff")
    plt.pcolormesh(xv,yv,(outsMax-outsMin)*mM_to_uM)
    plt.ylabel(label0)
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.colorbar()
    plt.xlabel(label1)         
    plt.tight_layout()
    
    name = state+"_extrema.png"
    plt.gcf().savefig(name)        

def PlotFluxes(t,j,idx1,label1="flux1",idx2=None,label2="flux2"):

  fig, ax1 = plt.subplots()
  rects = []
  labels= []
  #print "WARNING: should pull Cm, Vol, F from shannon model"
  #i_to_j = 2e-2 # [A/F] --> [uM/ms]
  labels.append(label1)
  rect1 =ax1.plot(t,j[:,idx1],'k',label=label1)
  rects.append(rect1[0])
  ax1.set_xlabel('time [ms]')
  ax1.set_ylabel(label1)
  
  
  if idx2!= None:
    ax2 = ax1.twinx()
    labels.append(label2)
    rect2 = ax2.plot(t,j[:,idx2],'k--',label=label2)
    rects.append(rect2[0])
    ax2.set_ylabel(label2)
    leg = ax2.legend( (rects) , (labels), loc=2, fancybox=True)
  else: 
    leg = ax1.legend( (rects) , (labels), loc=2, fancybox=True)
  
  
  
  leg.get_frame().set_alpha(1.0) 

idxCai = runner.model.state_indices("Cai")
idxCaSR = runner.model.state_indices("Ca_SR")
idxjRyR = runner.model.monitor_indices("j_rel_SR")

# Compute quantities of interest from transient data
def ProcessTransients(case,pacingInterval,tstart=8000):
    # compute rate of \catwo transient decline 
    tau = tf.GetTau(case,pacingInterval,tstart=tstart,idxCai=idxCai)
    
    
    # get transient amplitudes 
    # grab suitable spot for statistics 
    tsub, caisub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCai)
    delCa = (np.max(caisub) - np.min(caisub))*mM_to_uM
    
    # get SR transient amplitudes 
    tsub, casrsub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCaSR)
    delCaSR = (np.max(casrsub) - np.min(casrsub))*mM_to_uM
    
    # get max RyR
    tsub, jRyRsub= tf.GetInterval(case,pacingInterval,tstart=tstart,getFlux=True, idx=idxjRyR)
    #plt.figure()
    #plt.plot(jRyRsub)
    #jRyR = j[:, runner.model.monitor_indices(idx) ]
    maxRyR = np.max(jRyRsub)*mM_to_uM
    
    print tau, delCa, delCaSR, maxRyR
    #pctChg = GetExtreme(t,cai,subMin=7900,subMax=8500,si=7910)
    #pctChgSR = GetExtreme(t,casr,subMin=7900,subMax=8500,si=7910)
    
    return tau,delCa, delCaSR, maxRyR

def TransientBarPlots(cases,caseNames,resultsA,resultsB=False,tag=""):
    plt.subplot(1,3,1)
    idxs = np.arange(len(cases))
    width = 0.5
    plt.bar(idxs,resultsA.pctChgs,width)
    if resultsB:
      plt.bar(idxs+width,resultsB.pctChgs,width,color="r")
    plt.title("$\Delta Ca_i^{2+}$")
    plt.ylabel('$\Delta$ Ca [uM]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)

    plt.subplot(1,3,2)
    plt.bar(idxs,resultsA.pctChgSRs,width)
    if resultsB:
     plt.bar(idxs+width,resultsB.pctChgSRs,width,color="r")
    plt.title("$\Delta Ca_{SR}^{2+}$")
    plt.ylabel('$\Delta$ Ca [uM]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)


    ax = plt.subplot(1,3,3)
    rects1 = ax.bar(idxs,resultsA.taus,width)
    if resultsB:
      rects2 = ax.bar(idxs+width,resultsB.taus,width,color="r")
    plt.title("Decay constant")
    plt.ylabel('tau [ms]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)
    #plt.legend()
    plt.ylim([0,300])
    #label1="%3.1f Hz"%(1000/A.)
    #label2="%3.1f Hz"%(1000/B.)
    #ax.legend((rects1[0], rects2[0]), (label1,label2),loc=2)

    plt.tight_layout()
    plt.gcf().savefig(tag+"miscdata.png",dpi=300)

    plt.figure()
    plt.subplot(1,3,1)
    plt.bar(idxs,resultsA.maxRyRs,width)
    if resultsB:
      plt.bar(idxs+width,resultsB.maxRyRs,width,color="r")
    plt.title("$RyR Max$")
    plt.ylabel('$jRyR$ [uM/ms]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)


    plt.tight_layout()
    plt.gcf().savefig(tag+"jRyRmiscdata.png",dpi=300)

class empty:pass
mM_to_uM = 1e3
ms_to_s = 1e-3
# Collect all data 
def ProcessAllTransients(cases,caseTags,pacingInterval,tag="",cols=[],name=None):
    
    baseline = readPickle(caseTags[0]+tag+".pickle")
    caseA = readPickle(caseTags[1]+tag+".pickle")
    caseB= readPickle(caseTags[2]+tag+".pickle")
    #pca1p25vmax1p25 = readOut("PCa1.25ks1.00vMax1.25"+tag+".pickle")
    if len(cases)>3:
      caseC = readPickle(caseTags[3]+tag+".pickle")

    
    taus=[]
    pctChgs=[]
    pctChgSRs=[]
    maxRyRs=[]
    for i, case in enumerate(cases):
      var = eval(case)
      tau,pctChg,pctChgSR,maxRyR  = ProcessTransients(var,pacingInterval)
      taus.append(tau)
      pctChgs.append(pctChg)
      pctChgSRs.append(pctChgSR)  
      maxRyRs.append(maxRyR)    

    plt.figure()    
    plt.subplot(1,2,1) 
    for i, case in enumerate(cases): 
      var = eval(case)
      s = var['s']
      t = var['t']
      cai = s[1:,idxCai]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,cai[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("Ca [uM]")          
    #plt.legend(loc=3,ncol=2)
    plt.title("Cai")
    plt.ylim([0,1.2]) 
    
    plt.subplot(1,2,2)    
    for i, case in enumerate(cases):
      var = eval(case)  
      s = var['s']
      t = var['t']
      casr = s[1:,idxCaSR]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,casr[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("Ca [uM]")
    plt.ylim([250,1e3])    
    plt.legend(loc=2,ncol=1)
    plt.title("CaSR")

    plt.tight_layout()
    plt.gcf().savefig(name,dpi=300)

    plt.figure()
    for i, case in enumerate(cases):
      var = eval(case)  
      j = var['j']
      t = var['t']
      jRyR = j[:,idxjRyR]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,jRyR[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("j [uM/ms]")
    plt.legend(loc=0,ncol=1)      
    plt.ylim([0,200])    
    plt.xlim([8.5,8.7])    
    plt.title("jRyR")
    plt.tight_layout()
    plt.gcf().savefig("jRyR"+name,dpi=300)
    
    
    
    results = empty()
    results.taus = taus; 
    results.pctChgs = pctChgs; 
    results.pctChgSRs = pctChgSRs
    results.maxRyRs = maxRyRs
    
    return results



