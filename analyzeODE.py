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

class empty:pass
mM_to_uM = 1e3
ms_to_s = 1e-3

### 
### I/O 
###  
def readPickle(name = "PCa0.75kss0.25.pickle"):
  print "Reading " + name  
  pkl_file = open(name, 'rb')
  data1 = pickle.load(pkl_file)  
  pkl_file.close()

  return data1  

def LoadPickles(caseDict):
  for key,case in caseDict.iteritems():
    print "# ", key
    print "Loading "  , case.name
    case.data = readPickle(case.name)
    

###
### Mostly plotting 
### 
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
    
  pi = data1['p']
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


# trange can set 't' limit
def PlotPickleData(data1,data2=None,idxName="V",ylabel="V (mV)",trange=None,
    case1legend = None, case2legend=None
    ):    
  #  idx1=runner.model.state_indices(idxName)     
  # fluxes
  ms_to_s = 1e-3
  class empty:pass

  def mycont(data): 
    datac = empty()
    datac.t = data['t'] * ms_to_s
    datac.s = data['s']
    datac.s_idx = data['s_idx']
    datac.j = data['j']
    datac.j_idx = data['j_idx']

    if idxName in datac.j_idx:
      datac.v = datac.j
      datac.v_idx = datac.j_idx
    # states 
    elif idxName in datac.s_idx:
      datac.v = datac.s
      datac.v_idx = datac.s_idx
    else:
      print idxName, " not found"
      datac.v =None

    return datac

  datac1 = mycont(data1) 
  if data2!=None:
    datac2 = mycont(data2) 


  fig = plt.figure()

  if trange==None:
    fig.add_subplot(111)

  else:
    trange = np.asarray(trange) 
    plt.subplot(1,2,2)
    if datac1.v !=None:
      idx1 = datac1.v_idx.index(idxName)
      plt.plot(datac1.t,datac1.v[:,idx1],'k',label=case1legend)
    if data2!=None and datac2.v !=None:
      idx2 = datac2.v_idx.index(idxName)
      plt.plot(datac2.t,datac2.v[:,idx2],'r',label=case2legend)
    plt.xlim(trange*ms_to_s)
    plt.subplot(1,2,1)


  if datac1.v !=None:
    idx1 = datac1.v_idx.index(idxName)
    plt.plot(datac1.t,datac1.v[:,idx1],'k',label=case1legend)
  if data2!=None and datac2.v !=None:
    idx2 = datac2.v_idx.index(idxName)
    plt.plot(datac2.t,datac2.v[:,idx2],'r',label=case2legend)
  plt.xlabel('time [s]')
  plt.ylabel(ylabel)
  plt.legend(loc=3)
  plt.tight_layout()
  
  


def PlotFluxes(t,j,idx1=None,idx1Name="i_Ca",label1="flux1",idx2=None,label2=None):      

  raise RuntimeError("No longer using this. Try/revise PlotPickle()")
  if idx1==None:
    idx1=runner.model.monitor_indices(idx1Name)    
  #if idx1==None:
  #  idx1=runner.model.monitor_indices(label1)      

  fig, ax1 = plt.subplots()
  rects = []
  labels= []
  #print "WARNING: should pull Cm, Vol, F from shannon model"
  #i_to_j = 2e-2 # [A/F] --> [uM/ms]
  labels.append(label1)
  #rect1 =ax1.plot(t,j[:,idx1],'k',label=label1)
  rect1 =ax1.plot(t,j[idx1,:],'k',label=label1)
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

    plotAsDelta=False
   
    plt.subplot(1,3,2)
    if plotAsDelta:
      plt.bar(idxs,resultsA.pctChgSRs,width)
      if resultsB:
       plt.bar(idxs+width,resultsB.pctChgSRs,width,color="r")
      plt.title("$\Delta Ca_{SR}^{2+}$")
      plt.ylabel('$\Delta$ Ca [uM]')
      plt.xticks(idxs+0.5*width, (caseNames),rotation=70)
    else:
      plt.bar(idxs,resultsA.minCaSRs,width,color='b')
      plt.bar(idxs,resultsA.pctChgSRs,width,color='g',bottom=resultsA.minCaSRs)
      #print resultsA.maxCaSRs
      #if resultsB:
      # plt.bar(idxs+width,resultsB.pctChgSRs,width,color="r")
      plt.title("$Ca_{SR}^{2+}$")
      plt.ylabel("Min/$\Delta Ca^{2+}$ [uM]")
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

# Collect all data 
def ProcessAllTransients(cases,caseTags,pacingInterval,tag="",\
                         cols=[],name=None,root=""):
    baseline = readPickle(root+caseTags[0]+tag+".pickle")
    caseA = readPickle(root+caseTags[1]+tag+".pickle")
    caseB= readPickle(root+caseTags[2]+tag+".pickle")
    #pca1p25vmax1p25 = readOut("PCa1.25ks1.00vMax1.25"+tag+".pickle")
    if len(cases)>3:
      caseC = readPickle(root+caseTags[3]+tag+".pickle")

    
    taus=[]
    pctChgs=[]
    pctChgSRs=[]
    maxRyRs=[]
    minCaSRs=[]
    maxCaSRs=[]
    for i, case in enumerate(cases):
      var = eval(case)
      tau,pctChg,pctChgSR,minCaSR,maxCaSR,maxRyR  = ProcessTransients(var,pacingInterval)
      print minCaSR, maxCaSR
      taus.append(tau)
      pctChgs.append(pctChg)
      pctChgSRs.append(pctChgSR)  
      maxRyRs.append(maxRyR)    
      minCaSRs.append(minCaSR)    
      maxCaSRs.append(maxCaSR)    

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
    results.maxCaSRs = maxCaSRs
    results.minCaSRs = minCaSRs
    
    return results


# for despa data 
def DespaPlots(caseDict, stateName=None, monitorName=None,doLegend=True,label="",loadPickle=False,\
               xlim=[7000,10000],ylim=None,
               root=None,normalize=False):
    if loadPickle:
      LoadPickles(caseDicts)
#print "python runShannonTest.py "+" ".join(healthyArgs)+" -name "+healthyName+" &"

    if ylim:
      ylimUpdate=False
    else:
      ylim=[1e9,-1e9]
      ylimUpdate=True
    for key, case in caseDict.iteritems():
      t = case.data['t']
      j = case.data['j']
      s = case.data['s']
       
      if stateName!=None:  
        idxCai = runner.model.state_indices(stateName)
        #sca = s[1:,idxCai] - np.min(s[-ind:,idxCai])
        sca = s[1:,idxCai]  
      elif monitorName!=None:
        idxCai = runner.model.monitor_indices(monitorName)
        #sca = s[1:,idxCai] - np.min(s[-ind:,idxCai])
        sca = j[:,idxCai]  

      print "%s Diast/Systolic %f/%f " % \
            (case.tag,np.min(sca[xlim[0]:xlim[1]]),np.max(sca[xlim[0]:xlim[1]]))

      
      if normalize:
        sca = sca - np.min(sca[xlim[0]:xlim[1]])
        sca = sca/np.max(sca[xlim[0]:xlim[1]])
      plt.plot(t,sca,label=case.label)
            
      # determine bounds   
      if ylimUpdate:
        ylim[0] = np.min([ylim[0],np.min(sca)])
        ylim[1] = np.max([ylim[1],np.max(sca)])

    #plt.title("%s transients (offset by diastolic [%s])" %(stateName,stateName))    
   # plt.ylabel("[%s] - min([%s])" %(stateName,stateName))
    if stateName!=None:  
      plt.title("%s transients" %(stateName))    
      plt.ylabel("[%s] [mM])" %(stateName))
      fileName = stateName  
    elif monitorName!=None:
      plt.title("%s " %(monitorName))    
      plt.ylabel(label)
      fileName = monitorName
        
    if doLegend:
      plt.legend(loc=0)    
    plt.xlabel("t [ms]")
    if xlim!=None:
      plt.xlim(xlim)
    plt.ylim(ylim)
    if root==None:
      root=""
  
    if normalize:
      normalize="Norm"
    else:
      normalize=""        
    plt.gcf().savefig(root+"despa_%s%s.png"%(fileName,normalize),dpi=300)

###
### Data processing 
###

## Compute quantities of interest from transient data
# transient decay rates, 
# changes in SR calcium 
# ryr flux rates 
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
    minCaSRsub = (np.min(casrsub))*mM_to_uM
    maxCaSRsub = (np.max(casrsub))*mM_to_uM
    
    # get max RyR
    tsub, jRyRsub= tf.GetInterval(case,pacingInterval,tstart=tstart,getFlux=True, idx=idxjRyR)
    #plt.figure()
    #plt.plot(jRyRsub)
    #jRyR = j[:, runner.model.monitor_indices(idx) ]
    maxRyR = np.max(jRyRsub)*mM_to_uM
    
    print tau, delCa, delCaSR, maxRyR
    #pctChg = GetExtreme(t,cai,subMin=7900,subMax=8500,si=7910)
    #pctChgSR = GetExtreme(t,casr,subMin=7900,subMax=8500,si=7910)
    
    return tau,delCa, delCaSR, minCaSRsub,maxCaSRsub,maxRyR

## Power spectral density stuff 
        
import scipy.fftpack as fftp
def doPSD(sca):
  dm = sca - np.mean(sca,axis=0)
  M = fftp.fft(dm,axis=0)
  psd = np.abs(M)**2
  return psd  


#psd = doPSD(sca)
#plt.plot(psd)

def PlotFrequencies(s1,statei=0):
    # real domain 
    plt.figure()
    plt.pcolormesh(s1.T)
    plt.figure()
    plt.plot(s1[:,statei],label="%d"%statei)
    plt.legend()
    
    # get PSD
 #   dm = s1 - np.mean(s1)
 #   S = fftp.fft(dm,axis=0)
 #   psd2 = np.abs(S*S)
    
    psd2 = doPSD(s1)
    
    plt.figure()
    plt.pcolormesh(psd2.T)
    plt.figure()

    energy = np.sum(psd2,axis=0)
    psdScaled = psd2/energy
#plt.plot(psd2[:,0])
    plt.plot(psdScaled[:,statei],label="%d"%statei)

    return psdScaled

def PSDAnaly(s1,ranger=[2,200],verbose=True):
    ## raw data 
    #plt.figure()
    #plt.plot(s1)
    
    #plt.figure()
    #pcolormesh(s1.T,cmap=cm.gray)
    #plt.colorbar()
    
    ## Demeaned
    dc = np.mean(s1,axis=0)
    dm = s1 - dc
    dm /= np.max(s1,axis=0)   
    if verbose: 
      plt.figure()
      plt.bar(np.arange(np.shape(dc)[0]),dc)
    #plt.figure()
    #pcolormesh(dm.T,cmap=cm.gray)
    #plt.colorbar()
    
    # Power Spectral Density 
    S = fftp.fft(dm,axis=0)
    psd2 = np.abs(S*S)
    # log to make peaks more observable 
    psd2 = np.log(psd2)
    #print np.shape(psd2)
    # grab low freq 
    psd2 = psd2[ranger[0]:ranger[1],:]
    
    if verbose:
      plt.figure()
      from matplotlib import cm 
      plt.pcolormesh(psd2.T,cmap=cm.gray)
      plt.colorbar()
    
    
    #plt.figure()
    #plt.plot(psd2[ranger[0]:ranger[1],0])
    return dc, psd2

#  indSS = 2e3 # collect statistics after this time point [ms] (looking for steady state)
def ProcessDecomp(caseDict, \
                  wanted1="baseline",wanted2="incrleak",
                  indSS=2e3, # [ms]
                  xlimV=[89e3,90e3],
                  root="./",
                  ranked=20
                 ):  
  wanted =[wanted1,wanted2]
  
  # we want to compare only two of the cases, so select these here    
  #wanted =["baseline1Hz","baseline0.25Hz"]# ,"2xincrleak0.25Hz"]
  #wanted =["baseline0.25Hz","2xincrleak0.25Hz"]
  subCaseDict = dict()
  for wantedi in wanted:
    for key, case in caseDict.iteritems():
    # will be out of order  
    #if key in wanted:
      if key == wantedi:
        print "Selecting ", key      
        subCaseDict[ case.tag ]  = case
        #print np.shape(case.data['s']) 

  LoadPickles(subCaseDict)

  # for two cases, perform psd analysis to obtain mean and PSD
  caseComp = []
  maxTimeStep = -1
  for key, case in subCaseDict.iteritems():
      print case.name
      t = case.data['t']
      j = case.data['j']
      s = case.data['s']

      #sub = s[1e4:5e4,0:10]
      #sub = s[1e4:5e4,]
      sub = s[indSS:maxTimeStep,]
      case.dc, case.psd2 = PSDAnaly(sub,verbose=False)
      caseComp.append(case)


  ## display comparison of transients 
  plt.figure()
  label1 = "Cai"
  label2 = "Ca_SR"
  #idx =  module.state_indices( label )
  #plt.plot(tsteps,results[:,idx],label = label)
  fig, ax1 = plt.subplots()
  ax2 = ax1.twinx()
  i=0
  cols = ["r","b"]
  for key, case in subCaseDict.iteritems():
      print case.name
      ti = case.data['t']
      si = case.data['s']
      s_idx = case.data['s_idx']
      idx1 = s_idx.index(label1)
      idx2 = s_idx.index(label2)
      ax1.plot(ti,si[:,idx1],cols[i]+"-",label = case.label)
      ax2.plot(ti,si[:,idx2],cols[i]+"--")
      i+=1

  lastT = ti[-1]
  ax1.set_ylim([0,1.5e-3])
  ax1.set_xlim([indSS,lastT])
  ax2.set_ylim([0.0,0.6])
  ax2.set_xlim([indSS,lastT])
  plt.title("Ca2+ transients") 
  ax1.set_ylabel("%s [mM]"%label1)
  ax2.set_ylabel("%s [mM]"%label2)
  plt.legend()
  plt.gcf().savefig(root+"transients_%s_%s.png"%(wanted1,wanted2),dpi=300)

  ## action potential 
  plt.figure()
  label1 = "V"
  fig, ax1 = plt.subplots()
  i=0
  cols = ["r","b"]
  for key, case in subCaseDict.iteritems():
      print case.name
      ti = case.data['t']
      si = case.data['s']
      s_idx = case.data['s_idx']
      idx1 = s_idx.index(label1)
      ax1.plot(ti,si[:,idx1],cols[i]+"-",label = case.label)
      i+=1

  lastT = ti[-1]
  #ax1.set_ylim([0,1.5e-3])
  ax1.set_xlim(xlimV)
  plt.title("Action potential") 
  ax1.set_ylabel("%s [V]"%label1)
  plt.legend()
  plt.gcf().savefig(root+"AP_%s_%s.png"%(wanted1,wanted2),dpi=300)

  
  
  
  stateChg = (caseComp[1].dc-caseComp[0].dc)/caseComp[0].dc
  sort_index = (np.argsort(np.abs(stateChg)))[::-1]
  #### CUT AND PASTE FROM GOTRANNED CODE (shannon_2004.ode)
  if 1:
      print "SHOULD BE ABLE TO SURMISE THIS DIRECTLY FROM NEW PICKLE FILES"
      state_inds = dict([("h", 0), ("j", 1), ("m", 2), ("Xr", 3), ("Xs", 4),        ("X_tos", 5), ("Y_tos", 6), ("R_tos", 7), ("X_tof", 8), ("Y_tof", 9),        ("d", 10), ("f", 11), ("fCaB_SL", 12), ("fCaB_jct1", 13), ("R", 14),        ("O", 15), ("I", 16), ("Ca_TroponinC", 17), ("Ca_TroponinC_Ca_Mg",        18), ("Mg_TroponinC_Ca_Mg", 19), ("Ca_Calmodulin", 20), ("Ca_Myosin",        21), ("Mg_Myosin", 22), ("Ca_SRB", 23), ("Na_jct1_buf", 24),        ("Na_SL_buf", 25), ("Na_jct1", 26), ("Na_SL", 27), ("Nai", 28),        ("Ca_Calsequestrin", 29), ("Ca_SLB_SL", 30), ("Ca_SLB_jct1", 31),        ("Ca_SLHigh_SL", 32), ("Ca_SLHigh_jct1", 33), ("Ca_SR", 34),        ("Ca_jct1", 35), ("Ca_SL", 36), ("Cai", 37), ("V", 38)])
  
  #### END
  
  # create reverse lookup
  revDict = dict()
  for key, case in state_inds.iteritems():
    #print key, case
    revDict[case]=key  
  idx=1
  
  # grabbing top-twenty modulated states
  bestidx = sort_index[0:ranked]
  beststates = []
  for i,idx in enumerate(bestidx):
          if idx not in revDict:
              raise ValueError("Unknown state: '{0}'".format(idx))
              
          #print revDict[idx],"pct %4.2f"%stateChg[idx],                           "0 %4.1e"% caseComp[0].dc[idx],"1 %4.1e"% caseComp[1].dc[idx]    
          beststates.append(revDict[idx])
          #print beststates[i]
          #indices.append(state_inds[state])
          
  
  
  # Plot comparative data 
  
  # In[50]:
  
  width=0.3        
  plt.figure()
  fig, ax = plt.subplots()
  ind = np.arange(ranked)        
  
  dc0s = caseComp[0].dc
  norm = 1/dc0s[bestidx]
  rects1 = ax.bar(ind, dc0s[bestidx]*norm, width,color='r')
  
  dc1s = caseComp[1].dc
  rects2 = ax.bar(ind+width, dc1s[bestidx]*norm, width,color='b')
  
  ax.set_xticks(ind+width)
  ax.set_xticklabels( beststates,rotation=90 )
  
  
  lb1 =caseComp[0].label 
  lb2= caseComp[1].label
  plt.title("%s vs %s" % (lb1,lb2))
  ax.legend( (rects1[0], rects2[0]), (lb1,lb2),loc=0 )
  ax.set_ylabel("%chg wrt WT")
  plt.tight_layout()
  
  
  #plt.gcf().savefig(root+versionPrefix+"comparative.png",dpi=300)
  plt.gcf().savefig(root+"comparative_%s_%s.png"%(wanted1,wanted2),dpi=300)
  
  
  # In[ ]:


# Compute quantities of interest from transient data
# DUPE # def ProcessTransients(case,pacingInterval,tstart=8000):
# DUPE #     # compute rate of \catwo transient decline 
# DUPE #     tau = tf.GetTau(case,pacingInterval,tstart=tstart,idxCai=idxCai)
# DUPE #     
# DUPE #     
# DUPE #     # get transient amplitudes 
# DUPE #     # grab suitable spot for statistics 
# DUPE #     tsub, caisub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCai)
# DUPE #     delCa = (np.max(caisub) - np.min(caisub))*mM_to_uM
# DUPE #     
# DUPE #     # get SR transient amplitudes 
# DUPE #     tsub, casrsub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCaSR)
# DUPE #     delCaSR = (np.max(casrsub) - np.min(casrsub))*mM_to_uM
# DUPE #     
# DUPE #     # get max RyR
# DUPE #     tsub, jRyRsub= tf.GetInterval(case,pacingInterval,tstart=tstart,getFlux=True, idx=idxjRyR)
# DUPE #     #plt.figure()
# DUPE #     #plt.plot(jRyRsub)
# DUPE #     #jRyR = j[:, runner.model.monitor_indices(idx) ]
# DUPE #     maxRyR = np.max(jRyRsub)*mM_to_uM
# DUPE #     
# DUPE #     print tau, delCa, delCaSR, maxRyR
# DUPE #     #pctChg = GetExtreme(t,cai,subMin=7900,subMax=8500,si=7910)
# DUPE #     #pctChgSR = GetExtreme(t,casr,subMin=7900,subMax=8500,si=7910)
# DUPE #     
# DUPE #     return tau,delCa, delCaSR, maxRyR
# DUPE # 
# DUPE # 




