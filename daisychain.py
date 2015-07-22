import analyzeODE as ao
# Grabs stuff from previous run 
#prevOut = "run_G_CaBk1.00_G_NaBk1.00_stim2000_1.pickle"
#prevNum=1 # could probably grab this from pickleName
def InitializeNextInSequence(prevOut,prevNum):
  # Determine new pickleoutName 
  nextNum=prevNum+1
  nextOut = prevOut.replace("_%d.pickle"%prevNum,"_%d.pickle"%nextNum)

  # Load in prev data 
  data = ao.readPickle(prevOut) 
  si = data['s']
  s_idx = data['s_idx']
  p = data['p']
  p_idx = data['p_idx']


  # create new dict with states/values
  stateDict = dict()
  v = zip(s_idx,si[-1,:]) # grabs  state values from last iteration 
  for i,g in enumerate(v):
    #print g[0]    
    stateDict[g[0]]=g[1]  

  # create new dict with params/values   
  paramDict = dict()
  v = zip(p_idx,p) # grabs  state values from last iteration 
  for i,g in enumerate(v):
    #print g[0]    
    paramDict[g[0]]=g[1]  
    
  return nextOut,nextNum,stateDict,paramDict   
#s_idx

import runShannonTest as rs
def daisychain(\
    odeName = "shannon_2004.ode",
    dt=0.1,
    dtn=10e3, # elapsed time [ms]
    iters = 3,   
    stim_period=1000.,
    mxsteps=None,
    outBase = "run_G_CaBk1.00_G_NaBk1.00_stim2000"):
  
  
  
  ### sample param dictionary, can add specific parameters here
  paramDict = dict() 
  paramDict["stim_period"] = stim_period
  
  # stateDict
  stateDict = None # use defaults for first iter 
  #if 1:
  #  rs.runParamsFast(odeName=odeName,name=nextName,
  #                 paramDict=None,stateDict=stateDict,dt=dt,dtn=dtn,stim_period=stim_period,mxsteps=mxsteps)
  
  # subsequent runs     
  for i in range(iters):
      # initialize
      if i==0:
        # first run 
        nextNum = 1
        nextName = outBase+"_%d.pickle"%nextNum
      # second run         
      else:    
        prevName = nextName
        prevNum = nextNum
        nextName,nextNum,stateDict,paramDict = InitializeNextInSequence(\
          prevName,prevNum)
      
      
      # hack
      #stateDict["V"]=50 works 
      rs.runParamsFast(odeName=odeName,name=nextName,
                       varDict=paramDict,stateDict=stateDict,dt=dt,dtn=dtn,\
                       stim_period=stim_period,mxsteps=mxsteps)

  # create list of pickle names 
  daiters = range(iters) 
  pickleNames = [ outBase+"_%d.pickle"%(i+1) for i in daiters ]

  return pickleNames



def concatenateTrajs(pickleNames):
  import numpy as np 
  
  allsi = []
  allt = []
  tprev=0
  for i,pickleName in enumerate(pickleNames):
    data = ao.readPickle(pickleName) 
    si = data['s']
    s_idx = data['s_idx']
    t = data['t']
  
    allsi.append(si)  # probably should pre-allocate 
    allt.append(t+tprev)    
  
    # update offset  
    tprev += t[-1]  
    print "t ",tprev  

# test 
  #v = [np.arange(5),np.arange(5)+4, np.arange(5)+9]
  #v= np.array(v)
  #print np.ndarray.flatten(v)
  
  # concatenate times 
  ts = np.array(allt)
  ts = np.array(ts)
  # not general....
  ts = np.ndarray.flatten(ts)
  #print np.shape(ts)
  
  
  # concatenate state values 
  sis = np.array(allsi)
  #print np.shape(sis)
  
  nSteps = np.prod([np.shape(sis)[0],np.shape(sis)[1]])
  nStates = np.shape(sis)[2]
  #print nSteps
  
  allsisi = np.zeros([nSteps,nStates])
  
  # there's a smarter way to hash this out....
  for i in range(nStates):
      sisi = np.ndarray.flatten(sis[:,:,i])    
      allsisi[:,i] = sisi

  return ts, allsisi
  
      
  
  

