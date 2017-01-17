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
name = "gen%dprog%d"%(generation, progenyNumber)
varDict = dict()
stateDict = dict()
odeName = "../shannon_2004.ode"
#odeName = "../simpleRyR.ode"


def worker(procnum):
    print 'I am number %d in process %d' % (procnum, getpid())
    name = "gen%dprog%d"%(generation, progenyNumber)
    varName = "PCa"
    varVal  = 1.
    varDict[varName] = np.float(varVal)
    rs.runParamsFast(odeName=odeName,name=name,
                       varDict=varDict,stateDict=stateDict)
    return procnum,getpid()


def doit(processes=3,jobs=5):
    pool = multiprocessing.Pool(processes = processes)
    #print pool.map(worker, range(5))
    outputs = dict( pool.map(worker, range(jobs)) )
    
    for key, value in outputs.iteritems():
      print key,value


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
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-run"):
      #arg1=sys.argv[i+1] 
      #doit(arg1)
      doit(processes=5,jobs=5)

      exit()
  





  raise RuntimeError("Arguments not understood")




