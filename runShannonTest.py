"""
For performing parameter sweeps and running shannon model with dictionaries of parameter values 
"""

# Revisions 
# Changed default stimulation to 701
import numpy as np
import runner 
runner.init()
idxNCX = runner.model.monitor_indices("i_NaCa")

#def namer(PCa,ks,vMax=None,stim=None):
def namer(var1Name, var1Val, var2Name=None,var2Val=None,stim_period=1000,tag=None):
    #loc = "/u1/huskeypm/srcs/wholecell/"
    loc = "./"
    name =  loc+"run"
    name+=  "_%s%3.2f"%(var1Name,var1Val)
    if var2Name!=None:    
      name+=  "_%s%3.2f"%(var2Name,var2Val)

    name+="_stim%d"%stim_period

    if tag!=None:
      name+= "_"+tag

    return name 



class empty:
    def __init__(self,p,s,ts,js):
      self.p = p
      self.s = s
      self.ts = ts
      self.js = js

def runParams(
  runner=None,
  varDict=None,       
  name = "out",
  stim_period = 1000,
  mxsteps = 1000,
  deltaT = 3000, # duration
  dt = 1. # interval 
  ):

  # rescale params 
  pi = runner.model.p
  if varDict!=None:
    for key, value in varDict.iteritems():
      print "Rescale %s by %f" %(key,value)
      idx = runner.model.parameter_indices(key)
      param=pi[idx]
      pi[idx]=param*value
  #ks=pi[runner.model.parameter_indices("ks")]
  #pi[runner.model.parameter_indices("ks")]=ks*np.float(rks)
  #V_max=pi[runner.model.parameter_indices("V_max_Jpump")]
  #pi[runner.model.parameter_indices("V_max_Jpump")]=V_max*np.float(rVmax)

  runner.model.p = pi
  (p,s,t,j)=runner.runner(dt=deltaT, dtn=dt, stim_period=stim_period,mxsteps=mxsteps)
#dummy = runner.plotting(p,sres,tsres,jsres,case="fast_healthy")
#res = empty(p,sres,tsres,jsres)
  data1 = {'p':p,'s':s,'t':t,'j':j}
  #print np.shape(data1['s'])

  import cPickle as pickle 
  if ".pickle" not in name:
    name += ".pickle"
  output = open(name, 'wb')
  pickle.dump(data1, output)
  output.close()

# Print's command lines for running param sweep
# FixedParm is used to pass in a modified parameter that is not being swept
# over. Could probably be generalized into another varDict
def GenSweptParams(varDict,\
    stim_period=1000,T=10000,fixedParm=None,fixedParmVal=None,nameTag=None):

  # create list of input args (for command line) 
  allArgs=[]
  allVars=[]
  keys=[]
  for key,value in sorted(varDict.items()):
    #print key, value                
    keys.append(key)

    var1vals = [np.float(x) for x in value]               
    var1s = np.linspace(var1vals[0],var1vals[1],
      np.int((var1vals[1]-var1vals[0])/var1vals[2])+1)
    allVars.append(var1s)
    #print var1s
  
    args1 = []
    for i, rvar1 in enumerate(var1s):
      args1.append("-var %s %f"%(key,rvar1))
    allArgs.append(args1)
    #print args1


    

  # cmd and timing 
  cmdpre = "python runShannonTest.py"
  cmdpre+= " -stim %d" % stim_period
  cmdpre+= " -T %d" % T                     


  cmdFixedParms = ""
  if fixedParm!=None:
    cmdFixedParms = " -var %s %f" % (fixedParm,fixedParmVal)
    


  # iter over one var
  names = []
  if len(allArgs)==1:
    for i, arg1 in enumerate(allArgs[0]): 
        var1 = (allVars[0])[i]
        name = namer(keys[0],var1,stim_period=stim_period,tag=nameTag)
        cmd = cmdpre
        cmd+= " "+arg1 
        cmd+= cmdFixedParms
        cmd+= " -name "+name 
        cmd+= " &"
        print cmd
        names.append(name)

  elif len(allArgs)==2:

    for i, arg1 in enumerate(allArgs[0]): 
      for j, arg2 in enumerate(allArgs[1]):
        var1 = (allVars[0])[i]
        var2 = (allVars[1])[j]
        name = namer(keys[0],var1,keys[1],var2,stim_period=stim_period,tag=nameTag)
        cmd = cmdpre
        cmd+= " "+arg1 
        cmd+= " "+arg2 
        cmd+= cmdFixedParms
        cmd+= " -name "+name 
        cmd+= " &"
        print cmd
        names.append(name)

  else:
   raise RuntimeError("Not supported") 

  return names,keys,allVars 



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
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -sweep nameVar startVal endVal incrVal " % (scriptName)
  msg+=" \n or\n "
  msg+="  %s -var nameVar val  " % (scriptName)
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
  pi = runner.model.p
  stim = 1000 # [ms] 
  dt = 1. # [ms] 
  name="out"
  deltaT = 10000 # [ms] 
  sweep = False
  varDict = dict()              
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if("-var" in arg):
      varName =sys.argv[i+1] 
      varVal =sys.argv[i+2] 
      varDict[varName] = np.float(varVal)
    if("-sweep" in arg):
      varName =sys.argv[i+1] 
      varVals =sys.argv[(i+2):(i+5)] 
      varDict[varName] = varVals
      sweep=True
      
    if(arg=="-dt"):
      dt=np.float(sys.argv[i+1])
      
    if(arg=="-T"):
      deltaT=np.float(sys.argv[i+1])
    if(arg=="-name"):
      name=sys.argv[i+1] 
    if(arg=="-stim"):
      stim=sys.argv[i+1] 

  # execute
  if sweep:
    GenSweptParams(varDict)# var1Name,var1Vals)
  else: 
    runParams(runner=runner,varDict=varDict,\
              name=name,deltaT=deltaT,dt=dt, stim_period = stim)
  
