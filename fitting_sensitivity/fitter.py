"""
Automated model fitting
"""



#from pathos.multiprocessing import ProcessingPool as Pool

import multiprocessing
from os import getpid
import sys
sys.path.append("../")
import runShannonTest as rs
import numpy as np 


##
## Run parameters
##
generation =1
progenyNumber= 2
odeName = "../shannon_2004.ode"
#odeName = "../simpleRyR.ode"

#
# Launches an ode simulation on a single core 
#
def worker(procnum):
    print 'I am number %d in process %d' % (procnum, getpid())
    name = "gen%dprog%d"%(generation, progenyNumber)
    varDict = dict()
    stateDict = dict()
    varName = "PCa"
    varVal  = 1.
    varDict[varName] = np.float(varVal)
    rs.runParamsFast(odeName=odeName,name=name,
                       varDict=varDict,stateDict=stateDict)
    return procnum,getpid()


#
# Simply launches #jobs of #numCores and collects results
#
def SimpleRun(numCores=3,jobs=5):
    pool = multiprocessing.Pool(numCores = numCores)
    #print pool.map(worker, range(5))
    outputs = dict( pool.map(worker, range(jobs)) )
    
    for key, value in outputs.iteritems():
      print key,value

#
# Launches randomized parameters and collects 
# results 
#
def ParameterSensitivity(
  numCores=5,  # number of cores over which jobs are run
  numRandomDraws=3  # number of random draws for each parameter
  ):
  
  # Define parameters to vary and their default 
  # values  
  # should probably initialize ode and grab 
  # defaults. Lazy here
  parmDict = dict()

  defaultVal = 5.4e-4; stdDev = 1e-5
  parmDict["PCa"] = ( defaultVal, stdDev )

  defaultVal = 5.311e-3; stdDev = 1e-4
  parmDict["V_max_Jpump"] = ( defaultVal, stdDev )
  

  # Create a list of jobs with randomized parameters
  jobList = []
  
  for key,values in parmDict.iteritems():
    print key
    
    ## generate random pertubrations
    # draw from normal distribution
    mu,sigma = values
    randomDraws = np.random.normal(mu,sigma,numRandomDraws)

    # create a list of jobs 
    print randomDraws
    listN = [{key,val} for val in randomDraws]
    jobList+=listN
    

  print jobList



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
  Automated fitting
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
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
    # calls 'SimpleRun' with the next argument following the argument '-validation'
    if(arg=="-simpleRun"):
      #arg1=sys.argv[i+1] 
      #SimpleRun(arg1)
      SimpleRun(numCores=5,jobs=5)
      exit()

    if(arg=="-parameterSensitivity"):
      ParameterSensitivity(numCores=5)
      exit()
  





  raise RuntimeError("Arguments not understood")




