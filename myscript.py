import sachse
import analyze
import numpy as np

hdfName = "i.h5" # where data will be stored
params = sachse.Params() # where parameters can be modified
params.T = 1000.
params.dt = 2.5
mode="2D_noSSL_simple"; ssl = False;
#mode ="2D_noSSL"; ssl = False;

concsCaClefts=[]
concsCas=[]
concsCaSSLs=[]
#tsi=[]

import copy
  


# range of diffusion constants we want to consider
iters = 10 
iters2 = 5 
vs = np.linspace(-2,1,iters)
#print vs
Ds = 10.**vs
phis=np.linspace(0.1,1.0,iters) 
Kds=np.linspace(-6,-4,iters) 
Kds=10**Kds
buff=10.**-6
#Kds = [0.000001]
#phis = [1.0]

def runit(arg="test"): #,Kdi=1.,phij=1.): 

 #phij = [0.6]
 #Kdi = [0.000001]
# Kdi=10**Kdi
 for i,Di in enumerate(Ds):
  if 0:  
    #parms = copy.deepcopy(params)  
    parms = sachse.Params()
    parms.T = 1000 
    parms.dt = 2.5 
    hdfName="job_D_%3.1f.h5"%Di    
  
    # apply diff const   
    parms.D_CleftCyto = Di
    parms.D_CleftSSL = Di 
    parms.D_SSLCyto = Di 
    
    sachse.tsolve(mode=mode,hdfName=hdfName,params=parms) 
  
#  for i,Kdi in enumerate(Kds):
#    for j,phij in enumerate(phis):
#      print "i ", i," j ",j
#      print "Kdi ",Kdi, " phij ",phij
#      params = sachse.Params()
#      print "params.DCa ",params.DCa
#      Des1=params.DCa/(1+buff/Kdi)
#      Des2=2*phij/(3-phij)

#      Des=Des1*Des2
#      params.DCa_SSL=Des
#      print "Des1 ", Des1, " Des2 ", Des2, " DES ", Des
#      #mode = "2D_noSSL" 
#      hdfName = "Des_DCa_SSL_Kd%3.1f_%3.1f.h5"%(i,phij)
#      print "hdfName  ",hdfName
#      sachse.tsolve(mode=mode,hdfName=hdfName,params=params)


for i,Kdi in enumerate(Kds):
  for j,phij in enumerate(phis):
    if i > 7: 
     print "i ", i," j ",j
     print "Kdi ",Kdi, " phij ",phij
     params = sachse.Params()
     Des1=params.DCa/(1+buff/Kdi)
     Des2=2*phij/(3-phij)
#    print "params.DCa ", params.DCa

     Des=Des1*Des2
     params.D_CleftCyto=Des
     params.D_CleftSSL=Des
     params.D_SSLCyto=Des
     print "Des1 ", Des1, " Des2 ", Des2, " DES ", Des
  #mode = "2D_noSSL" 
     fileName = "Des_otherDs_noSSL_%3.7f_%3.1f"%(Kdi,phij)
     print "hdfName  ",fileName
     tag = "2D_noSSL_simple"
#     sachse.tsolve(mode=mode,hdfName=hdfName,params=params)
     sachse.tsolve(debug=False,params=params,pvdName = "noSSL.pvd", hdfName=fileName+".h5",\
     mode=tag,reactions="simple",buffers=True)

def readit():
  concsCaClefts=[]
  concsCas=[]
  concsCaSSLs=[]
  tsi=0
 # for i,Di in enumerate(Ds):
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
#  Kdi = 1.
#  phij = 1.
  for i,arg in enumerate(sys.argv):
    # myscript -runMPI <arg>
#    if(arg=="-Kdi"):
#      Kdi = np.float(sys.argv[i+1])
#    if(arg=="-phij"):
#      phij = np.float(sys.argv[i+2])
    if(arg=="-runMPI"):
#      runit(Kdi=Kdi)
#      runit(phij=phij)
      quit()
    if(arg=="-readSingle"):
      readit()   
      quit()
  





  raise RuntimeError("Arguments not understood")




