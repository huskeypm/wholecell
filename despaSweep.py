# ### Configure single run cases 
  
class empty:pass
import analyzeODE as ao
import numpy as np
import matplotlib.pylab as plt
root = "/net/share/pmke226/data/150609_despa/"
versionPrefix = "150609_"
  
def init():
  pacing = 1.0 # Hz
  stim_period = 1000/pacing # [ms]
  T_sec = 5
  T = T_sec * 1e3 # [ms] 
  
  caseDict = dict()
  
  # baseline 
  case = empty()
  case.tag = "baseline" # tag for filename 
  case.label = "Baseline" # label for simulation when plotted
  case.args = [ # pass in cellml/gotran arguments (eg. in shannon_2004.ode)
  "-stim %d" % stim_period,
  "-T %d" % T
  ]
  case.name = root+"despa_%s_%d.pickle"%(case.tag,stim_period) # filename 
  caseDict[case.tag] = case # 'save' this case to dictionary (referenced by 'tag')
  
  # incr leak
  case = empty()
  case.tag = "incrleak"
  case.label = "HIP (baseline SERCA) "
  case.args = [
  "-stim %d" % stim_period, # stimulation period 
  "-T %d" % T,  # simulation length [ms]
  "-var G_CaBk 5. ", # set variable G_CaBk to 5.0
  "-var G_NaBk 5. ",
  ]
  case.name = root+"despa_%s_%d.pickle"%(case.tag,stim_period)
  caseDict[case.tag] = case
  
  
  # incr leak/descr SERCA
  case = empty()
  case.tag = "incrleakdecrSERCA"
  case.label = "HIP (reduced SERCA)"
  case.args = [
  "-stim %d" % stim_period,
  "-T %d" % T,
  "-var G_CaBk 5. ",
  "-var G_NaBk 5. ",
  "-var V_max_Jpump 0.75 "
  ]
  case.name = root+"despa_%s_%d.pickle"%(case.tag,stim_period)
  caseDict[case.tag] = case
  
  
  # incr leak/descr SERCA
  case = empty()
  case.tag = "decrSERCA"
  case.label = "75% SERCA (UCD)"
  case.args = [
  "-stim %d" % stim_period,
  "-T %d" % T,
  "-var V_max_Jpump 0.75 "
  ]
  case.name = root+"despa_%s_%d.pickle"%(case.tag,stim_period)
  caseDict[case.tag] = case
  
  
  for key,case in caseDict.iteritems():
      print "# ", key
      print "python runShannonTest.py "+" -jit ".join(case.args)+" -name "+case.name+" &"
  #print "python runShannonTest.py "+" ".join(healthyArgs)+" -name "+healthyName+" &"

  return caseDict


def process(caseDict, wanted1="baseline",wanted2="incrleak"):  
  wanted =[wanted1,wanted2]
  
     

      
      
  # we want to compare only two of the cases, so select these here    
  #wanted =["baseline1Hz","baseline0.25Hz"]# ,"2xincrleak0.25Hz"]
  #wanted =["baseline0.25Hz","2xincrleak0.25Hz"]
  indSS = 2e3 # collect statistics after this time point [ms] (looking for steady state)
  stateName = "Ca_SR"
  subCaseDict = dict()
  for key, case in caseDict.iteritems():
    if key in wanted:
        print "Selecting ", key      
        subCaseDict[ case.tag ]  = case
        #print np.shape(case.data['s']) 

  ao.LoadPickles(subCaseDict)

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
      case.dc, case.psd2 = ao.PSDAnaly(sub)
      caseComp.append(case)
  
  
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
  num=20
  bestidx = sort_index[0:num]
  beststates = []
  for i,idx in enumerate(bestidx):
          if idx not in revDict:
              raise ValueError("Unknown state: '{0}'".format(idx))
              
          print revDict[idx],"pct %4.2f"%stateChg[idx],                           "0 %4.1e"% caseComp[0].dc[idx],"1 %4.1e"% caseComp[1].dc[idx]    
          beststates.append(revDict[idx])
          print beststates[i]
          #indices.append(state_inds[state])
          
  
  
  # Plot comparative data 
  
  # In[50]:
  
  width=0.3        
  fig, ax = plt.subplots()
  ind = np.arange(num)        
  
  dc0s = caseComp[0].dc
  norm = 1/dc0s[bestidx]
  rects1 = ax.bar(ind, dc0s[bestidx]*norm, width,color='r')
  
  dc1s = caseComp[1].dc
  rects2 = ax.bar(ind+width, dc1s[bestidx]*norm, width,color='b')
  
  ax.set_xticks(ind+width)
  ax.set_xticklabels( beststates,rotation=90 )
  
  
  lb1 =caseComp[0].label 
  lb2= caseComp[1].label
  plt.title("%s vs %s" % (wanted1,wanted2))
  ax.legend( (rects1[0], rects2[0]), (lb1,lb2),loc=0 )
  ax.set_ylabel("%chg wrt WT")
  
  
  #plt.gcf().savefig(root+versionPrefix+"comparative.png",dpi=300)
  plt.gcf().savefig("comparative_%s_%s.png"%(wanted1,wanted2),dpi=300)
  
  
  # In[ ]:
  
  
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


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -init" % (scriptName)
  msg+="  %s -process case1 case2" % (scriptName)
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
    if(arg=="-init"):
      init()
      exit()
    if(arg=="-process"):
      caseDict=init()
      process(caseDict,sys.argv[i+1],sys.argv[i+2])
      exit()
  





  raise RuntimeError("Arguments not understood")




