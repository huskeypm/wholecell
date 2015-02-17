import runner 
import numpy as np
runner.init()


class empty:
    def __init__(self,p,s,ts,js):
      self.p = p
      self.s = s
      self.ts = ts
      self.js = js

def doit(\
  name = "out",
  rPCa=1,
  rks=1,
  stim_period = 441,
  mxsteps = 1000,
  dt = 3000):

  # rescale params 
  pi = runner.model.p
  PCa=pi[runner.model.param_indices("PCa")]
  pi[runner.model.param_indices("PCa")]=PCa*np.float(rPCa)
  ks=pi[runner.model.param_indices("ks")]
  pi[runner.model.param_indices("ks")]=ks*np.float(rks)

  runner.model.p = pi
  (p,s,t,j)=runner.runner(dt=dt, stim_period=stim_period,mxsteps=mxsteps)
#dummy = runner.plotting(p,sres,tsres,jsres,case="fast_healthy")
#res = empty(p,sres,tsres,jsres)
  data1 = {'p':p,'s':s,'t':t,'j':j}
  print np.shape(data1['s'])

  import cPickle as pickle 
  output = open(name+'.pickle', 'wb')
  pickle.dump(data1, output)
  output.close()


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
  pi = runner.model.p
  rPCa= 1.0
  rks= 1.0
  name="out"
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-rpca"):
      rPCa=sys.argv[i+1] 
    if(arg=="-rks"):
      rks=sys.argv[i+1] 
    if(arg=="-name"):
      name=sys.argv[i+1] 

  doit(rPCa=rPCa,rks=rks,name=name,dt=10000)
  