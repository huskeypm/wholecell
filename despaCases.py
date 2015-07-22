# ### Configure single run cases 
  
class empty:pass
import analyzeODE as ao
import numpy as np
import matplotlib.pylab as plt
root = "/net/share/pmke226/data/150609_despa/"
versionPrefix = "150609_"
  

def init():
  # params applying to all cases. 
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
      print "python runShannonTest.py -jit "+" ".join(case.args)+" -name "+case.name+" &"
  #print "python runShannonTest.py "+" ".join(healthyArgs)+" -name "+healthyName+" &"

  return caseDict



  
  
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
      ao.ProcessDecomp(caseDict,sys.argv[i+1],sys.argv[i+2])
      exit()
  





  raise RuntimeError("Arguments not understood")




