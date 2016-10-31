"""
For performing parameter sweeps and running shannon model with dictionaries of parameter values 
"""

# Revisions 
# Changed default stimulation to 701
import numpy as np
import runner 
import downSamplePickles
runner.init()
idxNCX = runner.model.monitor_indices("i_NaCa")
import analyzeODE as ao
import gotranJIT

class empty:
    def __init__(self,p,s,ts,js):
      self.p = p
      self.s = s
      self.ts = ts
      self.js = js

#def WritePickle(name,p,p_idx,s,s_idx,j,j_idx,t):
#  raise RuntimeError("WARNING: need to antiquate this function and use analyzeODE version")
#  #ao.writePickle(name,p,p_idx,s,s_idx,j,j_idx,t)
#  data1 = {'p':p,'s':s,'t':t,'j':np.asarray(j),\
#           'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
#  #print np.shape(j) 
#  import cPickle as pickle
#  if ".pickle" not in name:
#    name += ".pickle"
#  output = open(name, 'wb')
#  pickle.dump(data1, output)
#  output.close()
#  print "SUCCESS! Wrote output to ", name

#def namer(PCa,ks,vMax=None,stim=None):

def namer(var1Name, var1Val, var2Name=None,var2Val=None,stim_period=1000,tag=None):
    
    raise RuntimeError("This is old, use NamerBetter. If need help see BDS.")

    #loc = "/u1/huskeypm/srcs/wholecell/"
    loc   = "/net/share/bdst227/Despa/Despa_Simulations_Data/"
    name  =  loc+"run"
    name +=  "_%s_%3.2f"%(var1Name,var1Val)
    if var2Name!=None:    
      name += "_%s_%3.2f"%(var2Name,var2Val)

    name += "_stim_%d"%stim_period

    if tag != None:
      name += tag

    return name 

# Returns name of file to be submited
### Made by BDS on 10/30/2016 ###
# Made it so that name of file can handle more cases.
# Does percentage math for you and puts into name. 
# Switches out pesky "." for "p" so computers do not confuse tags.
def NamerBetter(temp, var1Name, var1Val, var2Name=None,var2Val=None,stim_period=1000):

    leak_Base = 7.539e-4
    nka_Base = 5.0
    SERCA_Base = 7.02e-3
    
    name  = "mouse_"
    name += "Temp_%3.2f_" %(temp)
    
    # Leak name value
    if var1Name == "G_CaBk":
        leak_Change = var1Val 
        leak_Percentage = (leak_Change / leak_Base) * 100
        name += "leak%3.2fpct_" %(leak_Percentage)
	print "Be carefull and double check that leak is correct!!!!!!"
    elif var2Name == "G_CaBk":
        leak_Change = var2Val 
        leak_Percentage = (leak_Change / leak_Base) * 100
        name += "leak%3.2fpct_" %(leak_Percentage)
        print "Be carefull and double check that leak is correct!!!!!!"
    else:
        name += "leak100pct_"
        
    # NKA name value
    if var1Name == "I_NaK_max":
        nka_Change = var1Val 
        nka_Percentage = (nka_Change / nka_Base) * 100
        name += "nka%3.2fpct_" %(nka_Percentage)
    elif var2Name == "I_NaK_max":
        nka_Change = var2Val 
        nka_Percentage = (nka_Change / nka_Base) * 100
        name += "nka%3.2fpct_" %(nka_Percentage)
    else:
        name += "nka100pct_"
        
    # SERCA name value
    if var1Name == "V_max_Jpump":
        SERCA_Change = var1Val 
        SERCA_Percentage = (SERCA_Change / SERCA_Base) * 100
        name += "SERCA%3.2fpct_" %(SERCA_Percentage)
    elif var2Name == "V_max_Jpump":
        SERCA_Change = var2Val 
        SERCA_Percentage = (SERCA_Change / SERCA_Base) * 100
        name += "SERCA%3.2fpct_" %(SERCA_Percentage)
    else:
        name += "SERCA100pct_"

    stim_period_Hz = stim_period / 1000
    name += "freq%dHz" %(stim_period_Hz)
    
    name = name.replace(".","p")

    return name


def runParams(
  runner=None,
  varDict=None,       
  name = "out",
  stim_period = 1000,
  mxsteps = 1000,
  deltaT = 3000, # duration
  dt = 1., # interval 
  downsampleRate = 1  
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

  # Print pickle file of all states 
  import cPickle as pickle 
  raise RuntimeError("Dying here before writing pickle, since need to store j_index etc. Why are you using the non JIT code anyway?!!") 
  if ".pickle" not in name:
    name += ".pickle"
  output = open(name, 'wb')
  pickle.dump(data1, output)
  output.close()

# Print's command lines for running param sweep
### Made by BDS on 10/24/2016 ###
# Generalized code that constant variables and variable variables can be passed in same dict. 
# One day try to make the code so can handle more than 2 variable variables.
def GenSweptParamsBetter(
    stateDict, # states that are 'swept' over
    varDict, # variables that are 'swept' over 
    Time = 10000, stim_period = 1000, iters = 3,
    odeName = "shannon_2004_mouse.ode",
    nameTag=".pkl",dt=0.1,downsampleRate=None):

    #Command_line_input_pre  = "nohup python daisychain.py"
    commandLineInputPre  = "python daisychain.py"
    commandLineInputPre += " -dt %f" % dt
    commandLineInputPre += " -jit"
    commandLineInputPre += " -stim %d" % stim_period
    commandLineInputPre += " -T %d" % Time
    commandLineInputPre += " -iters %d" % iters
    if downsampleRate != None:
        commandLineInputPre += " -downsampleRate %d" % downsampleRate
   
    for key,value in sorted(stateDict.items()):           
    	#print keys
        if len(value) == 1:
        	commandLineInputPre += " -state %s %f" % (key, value[0])
	else:
		raise RuntimeError("This code can't sweep over states at this time :(")
 
    # create list of input args (for command line)
    allArgs=[]
    allVars=[]
    Non_Fixed_keys=[]

    for key,value in sorted(varDict.items()):           
        #print keys
        if key == "T":
            temp = value[0]
	    print "Temp: ", temp
	if len(value) == 1:
            commandLineInputPre += " -var %s %f" % (key, value[0])
        else:        
            Non_Fixed_keys.append(key)
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
    
    # iter over one var
    names = []
    if len(allArgs)==1:
        for i, arg1 in enumerate(allArgs[0]):
       	    var1 = (allVars[0])[i]
            name = namer(Non_Fixed_keys[0],var1,stim_period=stim_period,tag=nameTag)
            commandLineInput  = commandLineInputPre
            commandLineInput += " " + str(arg1)
            commandLineInput += " -odeName " + str(odeName)
            commandLineInput += " -name " + name
            commandLineInput += " &"
            print commandLineInput
        
    elif len(allArgs)==2:
        for i, arg1 in enumerate(allArgs[0]):
            for j, arg2 in enumerate(allArgs[1]):
                var1 = (allVars[0])[i]
                var2 = (allVars[1])[j]
                name = NamerBetter(temp,Non_Fixed_keys[0],var1,Non_Fixed_keys[1],var2,stim_period)
                commandLineInput  = commandLineInputPre
                commandLineInput += " " + str(arg1)
                commandLineInput += " " + str(arg2)
                commandLineInput += " -odeName " + str(odeName)
                commandLineInput += " -name " + name
                commandLineInput += " &"
                print commandLineInput
		
    else:
        raise RuntimeError("Not supported")
        
    return names,keys,allVars

# FixedParm is used to pass in a modified parameter that is not being swept
# over. Could probably be generalized into another varDict
def GenSweptParams(
    varDict, # variables that are 'swept' over 
    varDictFixed=None, # variables that are passed in but only assume one value 
    odeName ="shannon_2004_mouse.ode",
    stim_period=1000,T=10000,
    fixedParm=None,fixedParmVal=None, 
    nameTag=None):

  raise RuntimeError("This is old, use GenSweptParamsBetter. If need help see BDS.") 

  if fixedParm!=None:
    raise RuntimeError("Antiquated; use fixedVarDict") 

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
  #cmdpre = "nohup python runShannonTest.py"
  cmdpre = "nohup python daisychain.py"
  cmdpre+= " -jit "
  cmdpre+= " -stim %d" % stim_period
  cmdpre+= " -T %d" % T                     


  cmdFixedParms = ""
  if varDictFixed!=None:
    for key,value in varDictFixed.items():
      cmdFixedParms += " -var %s %f" % (key, value)

  # iter over one var
  names = []
  if len(allArgs)==1:
    for i, arg1 in enumerate(allArgs[0]): 
        var1 = (allVars[0])[i]
        name = namer(keys[0],var1,stim_period=stim_period,tag=nameTag)
        cmd  = cmdpre
        cmd += " "+arg1 
        cmd += cmdFixedParms
        cmd += " -odeName "+odeName 
        cmd += " -name "+name 
        cmd += " &"
        print cmd
        names.append(name)

  elif len(allArgs)==2:

    for i, arg1 in enumerate(allArgs[0]): 
      for j, arg2 in enumerate(allArgs[1]):
        var1 = (allVars[0])[i]
        var2 = (allVars[1])[j]
        name = namer(keys[0],var1,keys[1],var2,stim_period=stim_period,tag=nameTag)
        cmd  = cmdpre
        cmd += " "+arg1 
        cmd += " "+arg2 
        cmd += cmdFixedParms
        cmd += " -odeName "+odeName 
        cmd += " -name "+name 
        cmd += " &"
        print cmd
        names.append(name)

  else:
   raise RuntimeError("Not supported") 

  return names,keys,allVars 

###
### JIT support  
###

# Shamelessly snagged from gotranrun.py 
def GetMonitored(module, ode,tsteps,results,model_params):
    monitored = []
    for expr in sorted(ode.intermediates+ode.state_expressions):
        #print expr.name
        if expr.name not in monitored: # safety
            monitored.append(expr.name)

    monitored_plot = []        
    for expr in sorted(ode.intermediates):
        #print expr.name
        if expr.name not in monitored_plot: # safety
            monitored_plot.append(expr.name)

    monitor_inds = np.array([monitored.index(monitor) \
                             for monitor in monitored_plot], dtype=int)
    monitored_get_values = np.zeros(len(monitored), dtype=np.float_)

    # Allocate memory
    plot_values = np.zeros((len(monitored_plot), len(results)))  
    
    for ind, (time, res) in enumerate(zip(tsteps, results)):
        module.monitor(res, time, model_params, monitored_get_values)
        plot_values[:,ind] = monitored_get_values[ monitor_inds ]    

    plot_values = np.transpose(plot_values) 
    
    return monitored_plot, plot_values

#tstop=20000 # works 
#tstop=180000 # kicks the bucket after 130 s
#dt = 0.1
#stim_period = 1000.  # works (after adjusting odeint params)
#stim_period = 500.
# WARNING: here that parameters are SET, not rescaled, in contrast to runParams function
def runParamsFast(odeName = "shannon_2004.ode",name="out",\
                  varDict=None,stateDict=None,dt=0.1,dtn=2000,stim_period=1000.,mxsteps=None,downsampleRate=1 ):

  params = gotranJIT.init()
  params.tstop = dtn
  params.dt = dt
  #params.downsampleRate = downsampleRate

  # SET Params
  params.parameters = ['stim_period',stim_period]  
  if varDict!=None:
    for key, value in varDict.iteritems():
      params.parameters.extend([key,value])  

  # SET states 
  if stateDict!=None:   
    params.init_conditions=[]
    for key, value in stateDict.iteritems():
      #print key, value 
      params.init_conditions.extend([key,value])  
    
  ## JIT compile module, simulate  
  #ks=25. # default  [1/ms]
  #KSRleak=5.348e-6 # default [1/ms]
  #params.parameters = ['ks',ks,'KSRleak',KSRleak]  
  #print params.parameters    
  results,module,tsteps,model_params,ode=gotranJIT.main(odeName, params)  
  
  # to ensure consistency with old code 
  s = results
  t = tsteps
  p = model_params
  s_idx = [state.name for state in ode.full_states]
  p_idx = [param.name for param in ode.parameters]
  assert(len(s_idx)==len(results[0,:])),"Error in state vector. Ask PKH"
  assert(len(p_idx)==len(model_params)),"Error in parameter vector. Ask PKH"

  # get monitored fluxes   
  j_idx,j = GetMonitored(module, ode,tsteps,results,model_params)  
  print "j", np.shape(j)
  print "ji", len(j_idx)   
  
  if downsampleRate >1:     
      sDs,jDs,tDs = downSamplePickles.downsampleData(s,j,t,downsampleRate)
      print "jds", np.shape(jDs)
      print "name ", name 
      red = name 
      ao.writePickle(red,p,p_idx,sDs,s_idx,jDs,j_idx,tDs)
  else:
      ao.writePickle(name,p,p_idx,s,s_idx,j,j_idx,t)

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
  msg += "  %s -sweep nameVar startVal endVal incrVal " % (scriptName)
  msg += " \n or\n "
  msg += "  %s -var nameVar val  " % (scriptName)
  msg += """
  
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
  downsampleRate = 1  
  sweep = False
  useJIT=True # There shouldn't be a compelling reason to set this to false 
  varDict = dict()              
  odeName = "shannon_2004.ode"
  for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
    if("-odeName" in arg): 
      odeName = sys.argv[i+1]
      # wont work for useJIT=False
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
    #if(arg=="-downsampleRate"):
    if(arg=="-dSr" or arg=="-downsampleRate" ):
      downsampleRate = np.float(sys.argv[i+1])
    if(arg=="-name"):
      name=sys.argv[i+1] 
    if(arg=="-stim"):
      stim=sys.argv[i+1] 

    if(arg=="-jit"):
      useJIT = True

  # execute
  if sweep:
    if useJIT:
      raise RuntimeError("Not yet supported") 
    GenSweptParams(varDict)# var1Name,var1Vals)
  else: 
    if useJIT:
      runParamsFast(varDict=varDict,\
              odeName = odeName,
              name=name,dtn=deltaT,dt=dt, stim_period = stim, downsampleRate = downsampleRate)
    else:
      runParams(runner=runner,varDict=varDict,\
              name=name,deltaT=deltaT,dt=dt, stim_period = stim, downsampleRate = downsampleRate)
  
