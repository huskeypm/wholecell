import sachse
import analyze
import numpy as np

hdfName = "i.h5" # where data will be stored
params = sachse.Params() # where parameters can be modified
params.T = 10.
mode="2D_SSL"; ssl = True;
#mode ="2D_noSSL"; ssl = False;

concsCaClefts=[]
concsCas=[]
concsCaSSLs=[]
#tsi=[]

import copy
  


# range of diffusion constants we want to consider
iters = 4
vs = np.linspace(-2,1,iters)
#print vs
Ds = 10.**vs


def runit(arg="test"): 
  print arg 
  for i,Di in enumerate(Ds):
    parms = copy.deepcopy(params)  
    hdfName="job_D_%3.1f.h5"%Di    
  
    # apply diff const   
    parms.D_CleftCyto = Di
    parms.D_CleftSSL = Di 
    parms.D_SSLCyto = Di 
  
    sachse.tsolve(mode=mode,hdfName=hdfName,params=parms) 

def readit():
  concsCaClefts=[]
  concsCas=[]
  concsCaSSLs=[]
  tsi=0
  for i,Di in enumerate(Ds):
    hdfName="job_D_%3.1f.h5"%Di    
    print "I will read this ", hdfName 
    tsi,concsCaClefti,concsCaSSLi,concsCai=analyze.ReadHdf(hdfFile=hdfName,ssl=ssl)
    concsCaClefts.append(concsCaClefti)
    concsCas.append(concsCai)
    concsCaSSLs.append(concsCaSSLi)
  concsCaClefts = np.array(concsCaClefts)
  concsCas = np.array(concsCas)
  concsCaSSLs = np.array(concsCaSSLs)

  np.savetxt('tsi.txt',tsi)
  np.savetxt('concsCaClefts.txt',concsCaClefts)
  np.savetxt('concsCas.txt',concsCas)
  np.savetxt('concsCaSSLs.txt',concsCaSSLs)
      
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
    # myscript -runMPI <arg>
    if(arg=="-runMPI"):
      arg = sys.argv[i+1] 
      runit(arg=arg)   
      quit()
    if(arg=="-readSingle"):
      readit()   
      quit()
  





  raise RuntimeError("Arguments not understood")




