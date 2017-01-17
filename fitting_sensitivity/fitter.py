"""
Automated model fitting

Main idea: 
  User defines M parameters that are randomized
  User defines N outputs that are recorded 
  M 'threaded' jobs are launched on n<=M processors.
  After completion, N outputs from the M jobs are recorded 
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

##
## Inputs desired 
##
## Define parameters to vary and their default 
# values  
# should probably initialize ode and grab 
# defaults. Lazy here
parmDict = dict()

paramName = "PCa"
defaultVal = 5.4e-4; stdDev = 5e-5
parmDict[paramName] = ( defaultVal, stdDev )

paramName = "V_max_Jpump"
defaultVal = 5.311e-3; stdDev = 5e-4 
parmDict[paramName] = ( defaultVal, stdDev )


##
## Outputs desired
##
class outputObj:
    def __init__(self,name,mode):
      self.name = name 
      self.mode = mode 
      self.timeRange = [5e4,10e4]  # NEED TO ADD 
      self.result = None

#outputs = ["Cai","Nai"] 
outputList = { "Cai":outputObj("Cai","max"), 
               "Nai":outputObj("Nai","mean") }
             # decayRate:outputObj("decayRate","tau") 

class empty:pass


##
## Workers - these perform the 'threaded' tasks 
## 

#
# Launches an ode simulation on a single core 
#
def worker(procnum):
    print 'I am number %d in process %d' % (procnum, getpid())
    name = "gen%dprog%d"%(generation, progenyNumber)
    varDict = dict()
    varName = "PCa"
    varVal  = 1.
    varDict[varName] = np.float(varVal)
    rs.runParamsFast(odeName=odeName,name=name,
                       varDict=varDict)
    return procnum,getpid()


#
# Launch, same as worker, but accepts parameter dictionary
#
def workerParams(jobDict):
    jobNum =jobDict['jobNum']
    dtn = jobDict['jobDuration'] # [ms]
    varDict = jobDict['varDict']
    print "Worker bee %d, Job %d "%(getpid(),jobNum)

    verbose = False
    if verbose:
      for key,val in varDict.iteritems() :
        print "  ",key,val

    ## create varDict for runParams

    ## Launch job with parameter set 
    name = None  # don't write a pickle
    returnDict = dict() # return results vector 
    rs.runParamsFast(odeName=odeName,name=name,
                     varDict = varDict,
                     dtn=dtn,
                     returnDict=returnDict)


    ## do output processing 
    data = returnDict['data']
    outputResults = ProcessWorkerOutputs(data,outputList,tag=jobNum)
    if verbose: 
      for key,val in outputResults.iteritems() :
        print "  ",key,val.result


    ## package additional useful information 
    results = empty()
    results.outputResults = outputResults 
    results.pid = getpid()
    results.jobDict = jobDict
    results.jobNum = jobNum 
    
    
    return jobNum,results  



##
##  Data processing
## 

import analyzeODE as ao

# computes desired outputs based on data returned from ode model 
def ProcessWorkerOutputs(data,outputList,tag=99):
  outputResults = {}
  for key,obj in outputList.iteritems():
    dataSub = ao.GetData(data, obj.name)
    #print np.shape(dataSub.valsIdx)
    if key=="Cai": # for debug
      np.savetxt("test%d"%tag,dataSub.valsIdx)
    
    if obj.mode == "max":
        result = np.max(dataSub.valsIdx)
    elif obj.mode =="mean":
        result = np.mean(dataSub.valsIdx)
    else:
        raise RuntimeError("%s is not yet implemented"%output.mode)
    #output.result = result    
    resultObj = copy.copy(obj)
    resultObj.result = result
    #outputResults.append( resultObj ) 
    outputResults[key]=resultObj 

  return outputResults 
    




def test():
    name = None  # don't write a pickle
    returnDict = dict() # return results vector 
    rs.runParamsFast(odeName=odeName,name=name,
                     dtn=10,
                     returnDict=returnDict)

    #print returnDict['results']
    data = returnDict['data']
    print data.keys()
    dataSub = ao.GetData(data,"Cai")
    # probably also want to give this a range 
    output = np.mean(dataSub.valsIdx)
    print output
    #print vals
    #print np.mean(vals)

##
## Master processes that call workers
## 


#
# Simply launches #jobs of #numCores and collects results
#
def SimpleRun(numCores=3,jobs=5):
    pool = multiprocessing.Pool(processes= numCores)
    #print pool.map(worker, range(5))
    outputs = dict( pool.map(worker, range(jobs)) )
    
    for key, value in outputs.iteritems():
      print key,value

#
# Launches randomized parameters and collects 
# results 
#
import copy 
def ParameterSensitivity(
  numCores=5,  # number of cores over which jobs are run
  numRandomDraws=3,  # number of random draws for each parameter
  jobDuration = 2000 # job run time, [ms]
  ):

  
  
  ## Create 'master' varDict list 
  defaultVarDict= dict()
  numParams = 0
  for parameter,values in parmDict.iteritems():
    defaultVarDict[parameter]=values[0]  # default value
    print "Defaults: ", parameter, values[0]
    numParams+=1

  ## determine core count 
  numJobs = numRandomDraws*numParams
  numCores = np.min( [numCores, numJobs])
  print "Using %d cores for %d jobs"%(numCores,numJobs)

  ## Create a list of jobs with randomized parameters
  jobList = []
  ctr=0 
  for parameter,values in parmDict.iteritems():
    print parameter, " random draws:"
    
    ## generate random pertubrations
    # draw from normal distribution
    mu,sigma = values
    randomDraws = np.random.normal(mu,sigma,numRandomDraws)

    # create a list of jobs 
    print randomDraws
    #listN = [{parameter:val,'jobNum':i} for i,val in enumerate(randomDraws)]
    #jobList+=listN
    for val in randomDraws:
      varDict = copy.copy(defaultVarDict)
      varDict[parameter] = val 

      jobDict =  {'varDict':varDict,'jobNum':ctr,'jobDuration':jobDuration} 
      jobList.append( jobDict )
      ctr+=1
    
  #print jobList

  ## Run jobs
  pool = multiprocessing.Pool(processes= numCores) 
  jobOutputs = dict( pool.map(workerParams, jobList ) )
  #for key, results in outputs.iteritems():
  #  print key

  ## push data into a pandas object for later analysis
  myDataFrame = PandaData(jobOutputs,csvFile="example.csv")

  return myDataFrame

##
## I/O
## 

# Stores job information into a dict that can be used with pandas 
def StoreJob(job1):
    pandasDict = dict()
    tag = "%d_%d"%(job1.jobNum,job1.pid)
    pandasDict['jobID']=tag

    # pull out inputs
    varDict = job1.jobDict['varDict']
    for param,value in varDict.iteritems():
        #print param, value
        pandasDict[param] = value

    # pull out its results vector
    outputResults = job1.outputResults
    for output,result in outputResults.iteritems():
        #print output, result.result
        pandasDict[output] = result.result

    return pandasDict

import pandas as pd
#
# stores all data into a pandas object, which simplifies analyses later 
#
def PandaData(jobOutputs,csvFile="example.csv"):
  masterDict = dict()
 
  # get dictionary for each job and append it to a 'master' dictionary
  for workerNum, jobObj in jobOutputs.iteritems():
    jobDict = StoreJob(job1= jobObj)
    jobID = jobDict['jobID']
    masterDict[jobID]=jobDict

    
  # store data in pandas dataframe 
  df = pd.DataFrame(masterDict)    
  df.to_csv(csvFile)
  return df 



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

    if(arg=="-test"): 
      test()
      exit()
  





  raise RuntimeError("Arguments not understood")




