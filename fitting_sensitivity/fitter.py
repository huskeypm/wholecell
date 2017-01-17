"""
Automated model fitting
"""

import sys
sys.path.append("../")
import runShannonTest as rs
odeName = "ryr.ode"
generation =1
progenyNumber= 2 
name = "gen%dprog%d"%(generation, progenyNumber)
varDict = dict()
stateDict = dict()
odeName = "../shannon_2004.ode"

def doit():
  rs.runParamsFast(odeName=odeName,name=name,
                       varDict=varDict,stateDict=stateDict)




import multiprocessing

def worker(progenyNumber):
    """thread worker function"""
    print 'Worker:', progenyNumber
    name = "gen%dprog%d"%(generation, progenyNumber)
    rs.runParamsFast(odeName=odeName,name=name,
                       varDict=varDict,stateDict=stateDict)
    return

if __name__ == '__main__':
    jobs = []
    for progenyNumber in range(5):
        name = "gen%dprog%d"%(generation, progenyNumber)
        p = multiprocessing.Process(target=worker, args=(progenyNumber,))
        jobs.append(p)
        p.start()





