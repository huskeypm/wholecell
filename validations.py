from sachse import * 

# Verifying that we get the correct free Calcium if we merege 
# the SSL/Cyto into a single domain 
# TODO understand why there is still some numerical error. Did this at 5 am 
# Total conservation is not exact; might have to do with interpoaltion of the diract function
# used for defining the SSL region 
def validationMergingSSLCyto():
  params = Params()
  params.T = 1000   
  params.dt = 500

  idxCa = 0
  idxFluo= 2 # s.b. generalized 
  idxCaCleft = 4 # s.b. generalized 

  reactions = "ryrOnlySwitch"
  reactions = "ryrOnly"
  reactions = None     

  params.D_SSLCyto = 1e-1 # Can't go much faster than this   
  params.D_CleftSSL= 1e-1
  params.D_CleftCyto = 0.#1e-1
  params.dist = 0.01

  # kill Fluo
  params.Ftot = 0.
  params.cInits[idxFluo] = 0. 

  buffers = True  
  if buffers==False:
    params.Btot = 0.
    params.Ftot = 0.
 
  # separate SSL/Cyto compartments 
  case = "new"
  if 1: 
    mode = "2D_SSL"
    threeComps = tsolve(mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
      reactions = reactions,buffers=buffers) 
 
  # merged compartments 
  print "###\n###\n####\n"
  mode = "2D_noSSL"
  twoComps = tsolve(pvdName="test.pvd",mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
    reactions = reactions,buffers=buffers) 
  msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
  # NOTE: more generous with this error 
  assert(abs(threeComps[idxCa] - twoComps[idxCa]) < 1e-2), msg
  print "PASSED (but worth a couple check)" + msg 

# In[19]:

validationMergingSSLCyto()
