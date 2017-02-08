import sys
sys.path.append("./fitting_sensitivity/") 

import multiprocessing
from os import getpid
import runShannonTest as rs
import numpy as np
import analyzeODE as ao
import copy
import pandas as pd
import taufitting as tf
import matplotlib.pylab as plt
import fitter


class outputObj:
    #def __init__(self,name,mode):
    def __init__(self,name,mode,timeRange,truthValue):
      self.name = name
      self.mode = mode
      self.timeRange = timeRange #[5e4,10e4]  # NEED TO ADD 
      self.truthValue = truthValue
      self.result = None

#outputs = ["Cai","Nai"]
#outputListDefault = { "Nai":outputObj("Nai","mean"),
#                      "Cai":outputObj("Cai","max")} 
outputListDefault = { "Nai":outputObj("Nai","mean",[5e4,10e4],12.0e-3),
                      "Cai":outputObj("Cai","amp",[5e4,10e4],10000) }
             # decayRate:outputObj("decayRate","tau") 

class empty:pass

def workerParams(jobDict):
    #print "poop"
    odeName = "shannon_2004_mouse.ode"
    jobNum =jobDict['jobNum']
    dtn = jobDict['jobDuration'] # [ms]
    varDict = jobDict['varDict']
    print "Worker bee %d, Job %d "%(getpid(),jobNum)

    #print "varDict: ", varDict

    outputList = jobDict['outputList']
    #print "outputList: ", outputList
    #print "outputListDefault: ", outputListDefault
    if outputList == None:
        outputList = outputListDefault
        print "No outputList given, using outputListDefault."

    verbose = False
    if verbose:
      for key,val in varDict.iteritems() :
        print "  ",key,val

    ## create varDict for runParams
    #print "before runParamsFast"
    ## Launch job with parameter set 
    name = None  # don't write a pickle
    returnDict = dict() # return results vector 
    rs.runParamsFast(odeName=odeName, 
                     name=name,
                     varDict = varDict,
                     dt=0.1,
                     dtn=dtn,
                     stim_period=1000.0,
                     returnDict=returnDict)

    #print "after runParamsFast"

    ## do output processing 
    data = returnDict['data']
    #print "DATA: ", data
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

def ProcessWorkerOutputs(data,outputList,tag=99):
  outputResults = {}
  #print "made it to ProcessWorkerOutputs"
  #print "outputList: ", outputList
  for key,obj in outputList.iteritems():
    #print "key: ", key, "obj: ", obj
    #print "outputList: ", outputList
    #print "in the for loop"
    #print "obj.timeRange: ", obj.timeRange
    dataSub = ao.GetData(data, obj.name)

    #print "dataSub: ", dataSub
    #print "dataSub.valsIdx: ", dataSub.valsIdx
    result = ao.ProcessDataArray(dataSub,obj.mode,obj.timeRange,key=key)

    #output.result = result    
    resultObj = copy.copy(obj)
    resultObj.result = result
    #outputResults.append( resultObj ) 
    outputResults[key]=resultObj

  return outputResults

def PandaData(jobOutputs,csvFile="example.csv"):
  raise RuntimeError("Not using") 
  masterDict = dict()

  # get dictionary for each job and append it to a 'master' dictionary
  for workerNum, jobObj in jobOutputs.iteritems():
    jobDict = StoreJob(job1= jobObj)
    jobID = jobDict['jobID']
    masterDict[jobID]=jobDict


  # store data in pandas dataframe 
  df = pd.DataFrame(masterDict)
  df = df.T
  df.to_csv(csvFile)
  return df

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

# Genetic algorithm that randomizes the provided parameters (1 for now), selects the solution that minimizes the error, and repeats this process for a given number of iterations 
def fittingAlgorithm(
  myVariedParam, # Supports a single param currently = "Bmax_SL", 
  numCores=5,  # number of cores over which jobs are run
  numRandomDraws=3,  # number of random draws for each parameter
  jobDuration = 2000, # job run time, [ms]
  paramVarDict = None,
  outputList = None,
  truthValues = None,
  sigmaScaleRate = 1., # rate at which sigma is reduced by iteration 
  numIters = 10):
    
  trialParamVarDict = copy.copy( paramVarDict ) 
  
  iters = 0
  allErrors = []
  errorsGood_array = []
  flag = True
  randomDrawAllIters = []
  bestDrawAllIters = []
  while flag:

      ## Create 'master' varDict list 
      iters += 1

      numParams = 0

      defaultVarDict = dict()

      if trialParamVarDict != None:
          parmDict = trialParamVarDict
          
      print "parmDict: " , parmDict

      for parameter,values in parmDict.iteritems():
          defaultVarDict[parameter] = values[0]  # default value
          print "Inputs: ", parameter, values[0]
          numParams+=1

      ## determine core count 
      numJobs = numRandomDraws*numParams
      numCores = np.min( [numCores, numJobs])
      print "Using %d cores for %d jobs"%(numCores,numJobs)

      #print "outputList: ", outputList
      ## Create a list of jobs with randomized parameters
      jobList = []
      ctr=0
      for parameter,values in parmDict.iteritems():

          ## generate random pertubrations
          # draw from normal distribution
          mu,sigma = values
          #print "sigma: ", sigma
          print "Should probably rescale sigma by the tolerated error vs current error" 
          #rescaledSigma = sigma/(sigmaScaleRate * iters)
          rescaledSigma = sigma*np.exp(-sigmaScaleRate * (iters-1)) 
          print "rescaledSigma: ", rescaledSigma, " rate ", sigmaScaleRate
          #rescaledSigma = sigma
          distro = "lognormal"
          if distro=="normal":
            randomDraws = np.random.normal(mu,rescaledSigma,numRandomDraws)
          if distro=="lognormal":
            unif = np.random.normal(0,rescaledSigma,numRandomDraws)
            randomDraws = np.exp(unif) * mu 
            
  
          randomDraws = np.sort(randomDraws)

          # create a list of jobs 
          print parameter, " random draws:"
          print randomDraws
          randomDrawAllIters.append(randomDraws)
          #listN = [{parameter:val,'jobNum':i} for i,val in enumerate(randomDraws)]
          #jobList+=listN
          #print "JobList: ", jobList
          for val in randomDraws:
              varDict = copy.copy(defaultVarDict)
              varDict[parameter] = val

              jobDict =  {'varDict':varDict,'jobNum':ctr,'jobDuration':jobDuration, 'outputList':outputList}
              jobList.append( jobDict )
              ctr+=1
              #print "JobList2: ", jobList

      #print jobList

      ## Run jobs
      if numCores > 1:
          print "Multi-threading"
          pool = multiprocessing.Pool(processes = numCores)
          jobOutputs = dict( pool.map(workerParams, jobList))#, outputList ) )
      else:
          print "Restricting to one job only/assuming results are all that's needed"
          jobNum, results = workerParams(jobList[0])
      
      # Shouldn't have to write csv for these 
      myDataFrame = fitter.PandaData(jobOutputs,csvFile=None) # "example.csv")
      
      #allErrors.append([])
      #errorsGood_array.append([])
      
      jobFitnesses = np.ones( len(myDataFrame.index) )*-1
      jobNums      = np.ones( len(myDataFrame.index),dtype=int )*-1
      for i in range(len(myDataFrame.index)):
          #jobOutputs_copy = jobOutputs.copy()
          #slicedJobOutputs = jobOutputs_copy[slicer[]]
          #allErrors.append([myDataFrame.index[i]])
          #errorsGood_array.append([myDataFrame.index[i]])
          #print myDataFrame.index[i]
          #print myDataFrame.loc[myDataFrame.index[i],'jobNum']

          # score 'fitnesss' based on the squared error wrt each output parameter 
          fitness = 0. 
          for key,obj in outputList.iteritems():
              #print "outputList: ", key
              result = myDataFrame.loc[myDataFrame.index[i],key] 
              error = (result - obj.truthValue) ** 2
              #allErrors[iters-1].append(error)
              

              #if error <= (obj.truthValue * 0.001):
                  #errorsGood_array[iters-1].append(True)
              #else:
                  #errorsGood_array[iters-1].append(False)
              fitness += error
              
          # compute sqrt
          jobFitnesses[i] =  np.sqrt(fitness) 

          # These lines are intended to correct for a discrepancy between the pandas numbering and the job list
          # It works, but its confusing 
          jobNums[i] = myDataFrame.loc[myDataFrame.index[i],'jobNum']  
          myDataFrame.loc[myDataFrame.index[i],'fitness'] = jobFitnesses[i] 
      #
      # Summarize results
      #
      print "myDataFrame: "
      print myDataFrame
          
      print "jobFitnesses: ", jobFitnesses
      # find best job
      pandasIndex = np.argmin( jobFitnesses )
      jobIndex = jobNums[ pandasIndex ]           
      print "jobIndex: ", jobIndex
      #print "jobFitnes: " , jobFitnesses[jobIndex]
      # grab the job 'object' corresponding to that index
      bestJob = jobList[ jobIndex ]
      #print "bestJob: ", bestJob
      # get its input params/values
      bestVarDict = bestJob[ 'varDict' ]
      print "bestVarDict: " , bestVarDict
      
      variedParamVal = bestVarDict[ myVariedParam ]
      bestDrawAllIters.append(variedParamVal)
      # update 'trialParamDict' with new values, [0] represents mean value of paramater
      trialParamVarDict[ myVariedParam ][0]  = variedParamVal 
      # [1] to represent updating stdDev value
      #trialParamVarDict[ myVariedParam ][1]  = variedStdDevVal
      
      #print allErrors
      
      #if errorsGood_array[iters-1].count(False) == 0:
          #errorsGood = True
      #else:
          #errorsGood = False
          #print "Error is not good, need to run another iteration."
      
      #iters += 1
      
      print "iters: ", iters
      print ""
      print "######" 
      print ""
      
      if iters >= numIters: # or errorsGood:
          flag = False
      
  #return results

    #for key, results in outputs.iteritems():
    #  print key

  ## push data into a pandas object for later analysis
  #myDataFrame = PandaData(jobOutputs,csvFile="example.csv")

  #return myDataFrame
  return randomDrawAllIters, bestDrawAllIters    

# Here we try to optimize the sodium buffer to get the correct free Na concentration 
ms_to_s = 1e-3


def validation():
  # define job length and period during which data will be analyzed (assume sys. reaches steady state) 
  jobDuration = 4e3 # [ms] simulation length 
  timeRange = [1.0,jobDuration*ms_to_s] # [s] range for data (It's because of the way GetData rescales the time series)

  ## Define parameter, its mean starting value and the starting std dev 
  # Bmax_SL
  myVariedParam="Bmax_SL"
  paramDict = dict()  
  paramDict[myVariedParam] = [10.0, 1.0]
  
  ## Define the observables and the truth value
  outputList = {"Nai":outputObj("Nai","mean",timeRange,12.0e-3)}
  
  # Run 
  trial(paramDict=paramDict,outputList=outputList)


def test1():
  # define job length and period during which data will be analyzed (assume sys. reaches steady state) 
  jobDuration = 20e3 # [ms] simulation length 
  timeRange = [1.0,jobDuration*ms_to_s] # [s] range for data (It's because of the way GetData rescales the time series)

  ## Define parameter, its mean starting value and the starting std dev 
  # Bmax_SL
  myVariedParam="I_NaK_max"
  paramDict = dict()  
  truthVal = 5.0 
  #paramDict[myVariedParam] = [2*truthVal, 1.0]
  paramDict[myVariedParam] = [2*truthVal, 0.2] # for log normal
  sigmaScaleRate = .15
  
  ## Define the observables and the truth value
  outputList = {"Nai":outputObj("Nai","mean",timeRange,12e-3)}
  
  # Run 
  numRandomDraws = 10 
  numCores = np.min([numRandomDraws,30]) 
  numIters = 20 
  trial(paramDict=paramDict,outputList=outputList,numRandomDraws=numRandomDraws,numCores=numCores,numIters=numIters,sigmaScaleRate=sigmaScaleRate)

  


def trial(
  paramDict,
  outputList,
  numCores = 2, # maximum number of processors used at a time
  numRandomDraws = 2,# number of random draws for parameters list in 'parmDict' (parmDict should probably be passed in)
  jobDuration = 4e3, # [ms] simulation length 
  numIters=2,
  sigmaScaleRate = 1.
  ):

  # get varied parameter (should only be one for now) 
  keys = [key for key in paramDict.iterkeys()]
  variedParam = keys[0]


  
  ## do fitting and get back debugging details 
  allDraws,bestDraws = fittingAlgorithm(
    variedParam,numCores, numRandomDraws, jobDuration, paramDict, outputList,numIters=numIters, sigmaScaleRate=sigmaScaleRate)

  PlotDebuggingData(allDraws,bestDraws,numIters,numRandomDraws,title="Varied param %s"%variedParam) 
  
def PlotDebuggingData(allDraws,bestDraws,numIters,numRandomDraws,title=None):
  # put into array form 
  allDraws = np.asarray(allDraws)
  bestDraws = np.asarray(bestDraws)
  
  # create a matrix of random draws versus iteration
  vals= np.ndarray.flatten(allDraws)
  iters = np.repeat([np.arange(numIters)],numRandomDraws)
  scatteredData= np.asarray(zip(iters,vals))
  
  plt.scatter(scatteredData[:,0], scatteredData[:,1],label="draws")
  plt.plot(np.arange(numIters), bestDraws, label="best")
  plt.legend()
  if title!= None:
    plt.title(title) 
  plt.gcf().savefig("mytest.png")


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
def doit(fileIn):
  1


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
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
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      validation()
      quit()       
    if(arg=="-test1"):
      test1()
      quit()       
  





  raise RuntimeError("Arguments not understood")




