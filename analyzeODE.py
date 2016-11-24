"""
For processing ODE outputs in support of Satin/Despa collaborations 
"""
import cPickle as pickle
import runner 
import matplotlib.pylab as plt
import numpy as np
import analyzeODE as ao 
from runShannonTest import *
import taufitting as tf
runner.init()

class empty:pass
mM_to_uM = 1e3
ms_to_s = 1e-3

plotlyAuth=None
# For plotly support 
def PlotViaPlotly(casesSubset,state):
  downsample=10
  import plotly
  if plotlyAuth==None:
    with open('/net/share/pmke226/PLOTLY', 'r') as f:
        first_line = f.readline()
    plotlyKey = first_line
    plotly.tools.set_credentials_file(
      username='huskeypm', api_key='3x5rz5d19r') #plotlyKey)
  import plotly.tools as tls
  import plotly.plotly as py

  fig = plt.figure()
  ax = plt.subplot(111)
  case1=casesSubset[0]
  case2=casesSubset[1]
  title ="%s_%s_%s"%(state,case1.label,case2.label)
  ax.set_title("%s: %s,%s"%(state,case1.label,case2.label))
  for i,case in enumerate(casesSubset):
    pkg = GetData(case.data,state)
    ax.plot(pkg.t[::downsample],pkg.valsIdx[::downsample], label=case.label)
  plotly_fig = tls.mpl_to_plotly( fig )

  # Adding custom attributes to legend
  plotly_fig['layout']['showlegend'] = True
#layout = go.Layout(
#    xaxis=dict(
#        range=[2, 5]
#    ),
#    yaxis=dict(
#        range=[2, 5]
#    )
#)
#fig = go.Figure(data=data, layout=layout)

  plot_url = py.iplot(plotly_fig, filename = title)
  print plot_url.resource

### 
### I/O 
###  
def writePickle(name,p,p_idx,s,s_idx,j,j_idx,t):
  # store to pickle
  # using 'asarray' since my 'j' was getting stored as its transpose 
  data1 = {'p':p,'s':s,'t':t,'j':np.asarray(j),\
           'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
  #print "j again: ", len(j) 
  #print "j_idx: ",np.shape(j_idx)
  if ".pkl" not in name:
    name += ".pkl"
  output = open(name, 'wb')
  pickle.dump(data1, output)
  output.close()
  print "SUCCESS! Wrote output to", name

def readPickle(name = "PCa0.75kss0.25.pkl",verbose=True):          
  if verbose: 
    print "Reading " + name  
  pkl_file = open(name, 'rb')
  data1 = pickle.load(pkl_file)  
  pkl_file.close()

  return data1  

def LoadPickles(caseDict,noOverwrite=False,verbose=True):
  for key,case in caseDict.iteritems():
    if verbose:
      print "# ", key
      print "Loading"  , case.name

    if hasattr(case,'data') and noOverwrite==True:
      print "Skipping read, since already populated"
    else: 
      case.data = readPickle(case.name,verbose=verbose)
    
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
    return       

def GetData(data,idxName): 
    datac = empty()
    datac.t = data['t'] * ms_to_s
    datac.s = data['s'] / mM_to_uM
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

    idx = datac.v_idx.index(idxName)
    datac.valsIdx = datac.v[:,idx] 
    return datac

# trange can set 't' limit
def PlotPickleData_OLD(data1,data2=None,idxName="V",ylabel="V (mV)",trange=None,
    case1legend = None, case2legend=None,ylim=False,
    color='r'):    
  #  idx1=runner.model.state_indices(idxName)     
  # fluxes

  datac1 = GetData(data1,idxName) 
  if data2!=None:
    datac2 = GetData(data2,idxName) 

  fig = plt.figure()

  if trange==None:
    fig.add_subplot(111)

  else:
    trange = np.asarray(trange) 
    plt.subplot(1,2,2)
    if datac1.v !=None:
      idx1 = datac1.v_idx.index(idxName)
      plt.plot(datac1.t,datac1.v[:,idx1],color,label=case1legend)
    if data2!=None and datac2.v !=None:
      idx2 = datac2.v_idx.index(idxName)
      plt.plot(datac2.t,datac2.v[:,idx2],'k',label=case2legend)
    if ylim != False:
      plt.ylim(ylim)
    plt.xlim(trange*ms_to_s)
    plt.tight_layout()
    plt.subplot(1,2,1)

  if datac1.v !=None:
    idx1 = datac1.v_idx.index(idxName)
    plt.plot(datac1.t,datac1.v[:,idx1],color,label=case1legend)
  if data2!=None and datac2.v !=None:
    idx2 = datac2.v_idx.index(idxName)
    plt.plot(datac2.t,datac2.v[:,idx2],'k',label=case2legend)
  plt.xlim(0,60)
  plt.xlabel('time [s]',fontsize=14)
  if idxName == "Cai"or"Ca_SR":
      plt.ylabel(ylabel+' [uM]',fontsize=14)
  if idxName == "Nai":
      plt.ylabel(ylabel+' [mM]',fontsize=14)
  if idxName == "V":
      plt.ylabel(ylabel+' [mV]',fontsize=14)
  legend = plt.legend(loc=3)
  legend.get_frame().set_facecolor('white')
  plt.tight_layout()
  
### Below definition made by BDS on 10/13/2016 ###  
def Plot_Pickle_Data(rootOutput,datas,state=None,colors=None,xlabel=None,ylabel=None,
                   Full_image_xlim=None,Zoomed_image_xlim=None,plot_ylim=None,unit_scaler=None,
                   legends=None,time_range=None,legendMover1=1.0,legendMover2=1.0):

    Zoomed_image = plt.subplot(1,2,2)
    Full_image = plt.subplot(1,2,1)
    File_Name_Cases = ""
 
    def GetData(data,idxName):
    	datac = empty()
    	datac.t = data['t'] * ms_to_s
    	datac.s = data['s'] / mM_to_uM
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

    	idx = datac.v_idx.index(idxName)
    	datac.valsIdx = datac.v[:,idx]
    	return datac    

    for i,data in enumerate(datas):
        extracted_data = GetData(data,state)
        
        idx = extracted_data.v_idx.index(state)
        Full_image.plot(extracted_data.t,(extracted_data.v[:,idx]*unit_scaler),colors[i],label=legends[i])
        Zoomed_image.plot(extracted_data.t,(extracted_data.v[:,idx]*unit_scaler),colors[i],label=legends[i])
        File_Name_Cases += "%s_" %legends[i]

        Full_image.legend(loc=3)

    #plt.locator_params(nbins=6)
    Full_image.locator_params(nbins=8)
    Zoomed_image.locator_params(nbins=8)
    Full_image.set_xlim(Full_image_xlim)
    Zoomed_image.set_xlim(Zoomed_image_xlim)
    Full_image.set_ylim(plot_ylim)
    Zoomed_image.set_ylim(plot_ylim)
    
    plt.xlabel(xlabel, weight="bold",fontsize=14)
    plt.ylabel(ylabel, weight="bold",fontsize=14)
    plt.tight_layout()

    art = []
    lgd = plt.legend(bbox_to_anchor=(legendMover1,legendMover2))
    art.append(lgd)

    outFile = rootOutput+"Intracellular_%s_%splots.png"%(state,File_Name_Cases)
    print outFile
    plt.gcf().savefig(outFile,additional_artists=art,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

def GetData(data,idxName):
    datac = empty()
    datac.t = data['t'] * ms_to_s
    datac.s = data['s'] / mM_to_uM
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

    idx = datac.v_idx.index(idxName)
    datac.valsIdx = datac.v[:,idx]
    return datac

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

#idxCai = runner.model.state_indices("Cai")
#idxCaSR = runner.model.state_indices("Ca_SR")
#idxjRyR = runner.model.monitor_indices("j_rel_SR")


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
    
    ## Demeaned (looks for zeros) 
    dc = np.mean(s1,axis=0)
    dm = s1 - dc

    # abs. value 
    daMax = np.max(s1,axis=0)
    daMin = np.min(s1,axis=0)
    amp = daMax-daMin
    eps = 1e-9
    nonZero = np.argwhere(np.abs(daMax)>eps)
    dm[:,nonZero] = dm[:,nonZero]/daMax[nonZero] 
    isZero = np.argwhere(np.abs(daMax)<=eps)
    dm[:,isZero] = eps

    if verbose: 
      plt.figure()
      plt.bar(np.arange(np.shape(dc)[0]),dc)
      plt.legend()
    #plt.figure()
    #pcolormesh(dm.T,cmap=cm.gray)
    #plt.colorbar()
    
    # Power Spectral Density 
    S = fftp.fft(dm,axis=0)
    psd2 = np.abs(S*S)
    # log to make peaks more observable 
    psd2[ np.argwhere( np.abs(psd2) < eps )  ]= eps
    psd2 = np.log(psd2)
    #print np.shape(psd2)
    # grab low freq 
    psd2 = psd2[ranger[0]:ranger[1],:]
    
    if verbose:
      plt.figure()
      from matplotlib import cm 
      plt.pcolormesh(psd2.T,cmap=cm.gray)
      plt.colorbar()
      plt.legend()
    
    
    #plt.figure()
    #plt.plot(psd2[ranger[0]:ranger[1],0])
    return daMin,daMax,amp,dc, psd2

# State Decomposition Analysis (SDA)
# indSS = 2e3 # collect statistics after this time point [ms] (looking for steady state)
def ProcessDecomp(): 
  raise RuntimeError("ProcessDecomp is antiquated. Use StateDecompositionAnalysis()")

def StateDecompositionAnalysis(caseDict, \
                  wanted1="baseline",wanted2="incrleak",
                  indSS=[2e3,-1], # Time aftter which data is used for comparative analysis. -1 signifies getting last time step  [ms]
                  xlim=None,
                  root="./",
                  ranked=20,
                  ignoreList = ["dNa_SL_buf"],
                  mode="states" , # fluxes
                  cols = ["r","b"],
                  doPlot=False,
                  sortby="mean" # max, min, amp
                 ):  
  print "Starting"
  return  
  # we want to compare only two of the cases, so select these here    
  #wanted =["baseline1Hz","baseline0.25Hz"]# ,"2xincrleak0.25Hz"]
  #wanted =["baseline0.25Hz","2xincrleak0.25Hz"]
  subCaseDict = dict()
  wanted =[wanted1,wanted2]
  for idx, wantedi in enumerate(wanted):
    print wantedi
    for key, case in caseDict.iteritems():
    #if key in wanted:
      if key == wantedi:
        print "Selecting ", key      
        case.idx = idx 
        subCaseDict[ case.tag ]  = case
        #print np.shape(case.data['s']) 

  # Load data
  # Should be loaded dum dum 
  #LoadPickles(subCaseDict,noOverwrite=True)   


  # decide on which data to pull
  if mode=="states":
    v_key = 's'
  elif mode=="fluxes":
    v_key = 'j'
  else:
    raise RuntimeError(mode +" not understood") 
  

  ## for two cases, perform psd analysis to obtain mean and PSD
  caseComp = [None,None]
  for key, case in subCaseDict.iteritems():
      print case.name
      t = case.data['t']

      v = case.data[v_key] 

      #sub = s[1e4:5e4,0:10]
      #sub = s[1e4:5e4,]
      sub = v[int(indSS[0]):int(indSS[1]),]
      case.min,case.max, case.amp,case.dc, case.psd2 = PSDAnaly(sub,verbose=False)
      # abs 
      #case.max = np.abs( case.max ) 
      #case.dc  = np.abs( case.dc )   
      caseComp[ case.idx ] = case 


  ## displaycomparison
  # get shortest traj
  lastT = 1e30
  for key, case in subCaseDict.iteritems():
      ti = case.data['t']
      # pick whichever T is smallest: max length of either trajector or the xlim bound
      if xlim!=None:
        lastT = np.min([ti[-1],lastT,xlim[1]])
      else: 
        lastT = np.min([ti[-1],lastT])             
  if xlim==None:
    xlim = [0,lastT]

  def PlotValueComparison(label1,xmin=0):
    plt.figure()
    fig, ax1 = plt.subplots()
    i=0
    xlim=[xmin,xmin+1e4]
    #ylims = [1e90,-1e90]
    for key, case in subCaseDict.iteritems():
      ti = case.data['t']

      vi = case.data[v_key] 
      v_idx = case.data['%s_idx'%v_key]

      idx1 = v_idx.index(label1)

      #ti_vals = ti[xlim[0]:xlim[1]]             
      #vi_vals = vi[xlim[0]:xlim[1],idx1]
      ti_vals = ti; vi_vals=vi[:,idx1]
      ax1.plot(ti_vals,vi_vals,cols[i]+"-",label = case.label)

      i+=1

    #ax1.set_ylim([0,1.5e-3])
    ax1.set_xlim(xlim)
    #plt.title("Action potential") 
    ax1.set_ylabel("%s [unk]"%label1)
    plt.legend()
    plt.gcf().savefig(root+"%s_%s_%s_%s.png"%(mode,label1,wanted1,wanted2),dpi=300)

  ## Quantify (pct error) change in mean value of each state
  # normalize by case 1
  def donorm(subj,ref):  
    eps = 1e-9
    # means 
    nonZero = np.argwhere(np.abs(ref.dc)>eps)
    subj.dcn = np.zeros( np.shape(ref.dc) ) 
    subj.dcn[ nonZero ] = subj.dc[ nonZero ] / ref.dc[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.max)>eps)
    subj.maxn = np.zeros( np.shape(ref.max) ) 
    subj.maxn[ nonZero ] = subj.max[ nonZero ] / ref.max[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.min)>eps)
    subj.minn = np.zeros( np.shape(ref.min) ) 
    subj.minn[ nonZero ] = subj.min[ nonZero ] / ref.min[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.amp)>eps)
    subj.ampn = np.zeros( np.shape(ref.amp) ) 
    subj.ampn[ nonZero ] = subj.amp[ nonZero ] / ref.amp[ nonZero ]

    #maxdiff
    nonZero = np.argwhere(np.abs(ref.amp)>eps)
    subj.maxdiff= np.zeros( np.shape(ref.amp) )
    subj.maxdiff[ nonZero ] = subj.maxdiff[ nonZero ] - ref.maxdiff[ nonZero ]

  donorm(caseComp[0],caseComp[0])
  donorm(caseComp[1],caseComp[0])

  # sort from largest to smallest change 
  if sortby=="mean" or sortby=="dc":
    valueChg = caseComp[1].dcn-caseComp[0].dcn
  elif sortby=="max":
    valueChg = caseComp[1].maxn-caseComp[0].maxn
  elif sortby=="min":
    valueChg = caseComp[1].minn-caseComp[0].minn
  elif sortby=="amp":
    valueChg = caseComp[1].ampn-caseComp[0].ampn
  sort_index = (np.argsort(np.abs(valueChg)))[::-1]


  #Assuming that both pickle files have same states/ode model. 
  v_idx = case.data['%s_idx'%v_key]
  value_inds = {key: idx for (idx, key) in enumerate(v_idx)}
  # create reverse lookup
  value_inds_rev = {idx: key for (idx, key) in enumerate(v_idx)}
  
  # grabbing top-twenty modulated states
  bestidx    = []
  bestvalues = []
  lines = []
  lines.append("name caseno val diff_val\n")

  stored = 0 
  for i,idx in enumerate(sort_index):
          if idx not in value_inds_rev:
              raise ValueError("Unknown state/flux: '{0}'".format(idx))

          if value_inds_rev[idx] in ignoreList:
            #print "Skipping ", value_inds_rev[idx]
            continue
              
          # print raw values
          if sortby=="mean":
            line=value_inds_rev[idx]+" pct %4.2f "%valueChg[idx]+ \
                                      " 0 %4.1e %4.1e/%4.1e"% (caseComp[0].max[idx], caseComp[0].dc[idx],caseComp[0].dcn[idx])+\
                                      " 1 %4.1e %4.1e/%4.1e"% (caseComp[1].max[idx], caseComp[1].dc[idx],caseComp[1].dcn[idx])    
          elif sortby=="max":
            line = "MAX %-20s"%value_inds_rev[idx]+ " 0 %4.1e 1 %4.1e [%3.1f] "%( caseComp[0].max[idx], caseComp[1].max[idx],100*caseComp[1].max[idx]/caseComp[0].max[idx]-100)
          elif sortby=="min":
            line = "MEAN %-20s"%value_inds_rev[idx]+ " 0 %4.1e 1 %4.1e [%3.1f] "%( caseComp[0].min[idx], caseComp[1].min[idx],100*caseComp[1].min[idx]/caseComp[0].min[idx]-100)
          elif sortby=="amp":
            line = "AMP  %-20s"%value_inds_rev[idx]+ " 0 %4.1e = (%4.1e - %4.1e) 1 %4.1e [%3.1f] "%( caseComp[0].amp[idx], 
                                                                                   caseComp[0].max[idx],caseComp[0].min[idx],
                                                                                   caseComp[1].amp[idx],100*caseComp[1].amp[idx]/caseComp[0].amp[idx]-100)
          print line

          # store 
          lines.append(line)
          bestidx.append(idx)
          bestvalues.append(value_inds_rev[idx])
          #print beststates[i]
          #indices.append(state_inds[state])

          stored+=1
          if stored>=ranked and ranked>0:
            break 
          
  ## Plot State Data 
  if mode == "states":
    plotValues = ["V","Cai","Ca_SR"] + bestvalues
  elif mode =="fluxes":
    plotValues = ["i_Cab"] + bestvalues

  #
  if doPlot:
    for label in plotValues:
      PlotValueComparison(label,xmin=indSS[0])
  
  
  ## Plot comparative data 
  width=0.3        
  plt.figure()
  fig, ax = plt.subplots()
  ind = np.arange(stored)        
  
  # store info for later analysis
  caseComp[0].bestidx =  bestidx
  caseComp[0].bestvalues =  bestvalues
  caseComp[1].bestidx = bestidx
  caseComp[1].bestvalues = bestvalues
  if sortby=="mean":
    vals0s = caseComp[0].dcn
    vals1s = caseComp[1].dcn
  elif sortby=="max":
    vals0s = caseComp[0].maxn
    vals1s = caseComp[1].maxn
  elif sortby=="min":
    vals0s = caseComp[0].minn
    vals1s = caseComp[1].minn
  elif sortby=="amp":
    vals0s = caseComp[0].ampn
    vals1s = caseComp[1].ampn
 
  # bar plot 
  rects1 = ax.bar(ind, vals0s[bestidx], width,color='r')
  rects2 = ax.bar(ind+width, vals1s[bestidx], width,color='b')
  
  ax.set_xticks(ind+width)
  ax.set_xticklabels( bestvalues,rotation=90 )
  
  
  lb1 =caseComp[0].label 
  lb2= caseComp[1].label
  plt.title("%s vs %s (%s)" % (lb1,lb2,sortby))
  ax.legend( (rects1[0], rects2[0]), (lb1,lb2),loc=4 )
  ax.set_ylabel("fold chg wrt WT")
  plt.tight_layout()
  
  #plt.gcf().savefig(root+versionPrefix+"comparative.png",dpi=300)
  filePrefix = root+"comparative_%s_%s_%s_%s"%(mode,wanted1,wanted2,sortby)
  plt.gcf().savefig(filePrefix+".png",dpi=300)

  import json
  f = open(filePrefix+'.txt', 'w')
  json.dump(lines, f)
  f.close()

  return caseComp

def PlotMorotti(cases,
                case1Name='rabbit_5',case2Name='mouse_5',
                trange=[2.0e3,2.2e3],
                root ="./" # path for prionting figures 
                ):
    case1=cases[case1Name]
    case2=cases[case2Name]    
    states = ["V", "Cai", "Nai"]
    ctr=0
    for i,state in enumerate(states):
        PlotPickleData(case1.data,data2=case2.data,idxName=state, 
                          ylabel=state,trange=trange,
                          case1legend=case1.caseName,
                          case2legend=case2.caseName)
        plt.tight_layout()
        #title = case1.caseName+"_mouserabbit_comp"+state
        title = "mouse_rabbit_compare_%.2d_"%ctr+state
        plt.gcf().savefig(root+"/"+title+".png")
        ctr+=1

    fluxes = ["i_Na", "i_CaL", 
              "i_kur",# I think this is Morotti's I_K,slow
              "i_ss",
              "i_tof", "i_tos","i_Ks","i_Kr","i_K1",
              "i_NaCa","i_NaK",
              "i_Kp"]

    #I_Ca, IK,Slow, Iss, Ito, IKs, IKr, IK1, INaCa, INaK
    for i,flux in enumerate(fluxes):
        PlotPickleData(case1.data,data2=case2.data,idxName=flux,
                          ylabel="%s [A/F]"%flux,trange=trange,
                            case1legend=case1.caseName,
                          case2legend=case2.caseName)                        
        #title = case1.caseName+"_mouserabbit_comp"+flux
        title = "mouse_rabbit_compare_%.2d_"%ctr+flux
        ctr+=1
        plt.gcf().savefig(root+"/"+title+".png")
                
