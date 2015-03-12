"""
Validation routines for sachse.py
"""
from sachse import * 


def validationRyRSERCA():
  params = Params()
  params.T = 3      
  params.dt = 1   

  # this expression is taken from simple.py
  #self.iryr = Expression(dirac+"a*exp(-(t-to)/tau)",
  #                           v = 0.05,\
  #                           a=params.ryrAmp,\
  #                           to=params.ryrOffset,\
  #                           tau=params.ryrTau,\
  #                           t=0)
  
  si=0; fi=25; st = 50;
  dt = (fi-si)/float(st)
  ts = np.linspace(0,25,50)
  i = params.ryrAmp*np.exp(-t/tau)
  j = wholeCellI_to_wholeCellJ(wholeCellI)
  dCa = 
  iTot = np.cumsum(j)
  
  
  


  quit()

  idxCa = 0
  idxBuff = 1
  idxFluo= 2 # s.b. generalized 
  idxCaCleft = 4 # s.b. generalized 

  reactions = "ryrOnly"
  #params.ryrAmp = 40.; eps = 5e-2
  #reactions = None; eps = 1e-2

  params.D_SSLCyto = 1e2 # Can't go much faster than this   
  params.D_CleftSSL= 1e2
  params.D_CleftCyto = 1e2 
  params.dist = 0.01

  # kill Fluo
  params.Ftot = 0.
  params.cInits[idxFluo] = 0. 

  buffers = False 
  if buffers==False:
    params.Btot = 0.
    params.cInits[idxBuff] = 0. 
    params.Ftot = 0.
 
  # separate SSL/Cyto compartments 
  # probably not necesarru, since we test below 
  case = "tot"
  #if 1: 
  #  mode = "2D_SSL"
  #  threeComps = tsolve(mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
  #    reactions = reactions,buffers=buffers) 
 
  # merged compartments 
  print "###\n###\n####\n"
  mode = "2D_noSSL"
  twoComps = tsolve(pvdName="test.pvd",mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
    reactions = reactions,buffers=buffers) 
  msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
  # NOTE: more generous with this error 
  assert(abs(threeComps[idxCa] - twoComps[idxCa]) < eps), msg

  reactions = "ryrOnly"

# Verifying that we get the correct free Calcium if we merege 
# the SSL/Cyto into a single domain 
# TODO understand why there is still some numerical error. Did this at 5 am 
# Total conservation is not exact; might have to do with interpoaltion of the diract function
# used for defining the SSL region 
def validationMergingSSLCyto():
  params = Params()
  params.dt = 1.0 

  idxCa = 0
  idxBuff = 1
  idxFluo= 2 # s.b. generalized 
  idxCaCleft = 4 # s.b. generalized 

  reactions = "ryrOnlySwitch"
  params.ryrAmp = 20. ; eps = 5e-2
  #reactions = "ryrOnly"
  #reactions = None; eps = 1e-2

  params.D_SSLCyto = 1e2 # Can't go much faster than this   
  params.D_CleftSSL= 1e2
  params.D_CleftCyto = 1e2 
  params.dist = 0.01

  # kill Fluo
  params.Ftot = 0.
  params.cInits[idxFluo] = 0. 

  buffers = True  
  if buffers==False:
    params.Btot = 0.
    params.cInits[idxBuff] = 0. 
    params.Ftot = 0.
 
  # separate SSL/Cyto compartments 
  case = "new"
  if 1: 
    params.T = 25     
    mode = "2D_SSL"
    threeComps = tsolve(mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
      reactions = reactions,buffers=buffers) 
 
  # merged compartments 
  print "###\n###\n####\n"
  mode = "2D_noSSL"
  params.T = 76     
  params.dt = 0.5
  twoComps = tsolve(pvdName="test.pvd",mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
    reactions = reactions,buffers=buffers) 
  msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
  # NOTE: more generous with this error 
  assert(abs(threeComps[idxCa] - twoComps[idxCa]) < eps), msg
  print "PASSED (but worth a couple check)" + msg 

# In[19]:
def validation_Conversions():
  def dotest(cCaFinal,cCaInit,name=""):
    refConcChange = params.T*params.jTest # [uM]
    concChange = cCaFinal - cCaInit
    msg = "%f != %f " % ( concChange,refConcChange)
    assert( abs(concChange - refConcChange) < 1e-5 ), msg
    print "Passed flux test " + name

  idxCa = 0
  idxCaCleft = 4 # NEED TO NOT HARDCODE

  ## validate wholecell->cyto flux test
  if 1:
   params = Params()
   params.T = 10
   concsFinal = tsolve(doAssert="fluxTest_jVolCyto",\
                      mode="2D_noSSL",
                      params=params,reactions=None,buffers=False, 
                      existsCleft=False,existsSSL=False) 
   dotest(concsFinal[idxCa],params.cCaInit,name="jVolCyto")

  ## Validate wholecell flux test 
  if 1: 
   params = Params()
   params.T = 10
   params.cCaInit = 0.1
   concsFinal = tsolve(doAssert="fluxTest_jSurfCyto",\
                      mode="2D_noSSL",
                       params=params,reactions=None,buffers=False,
		       existsCleft=False,existsSSL=False) 
   dotest(concsFinal[idxCa],params.cCaInit,name="jSurfCyto") 


  ## Validate A/F conversion 
  # apply flux to cytosol domain 
  if 1: 
   iFlux = 1 # [A/F]
   params = Params()
   params.jTest = wholeCellI_to_wholeCellJ(iFlux) 
   params.jTest*=4  # to make about 0.1 uM/ms
   params.T = 10
   concsFinal = tsolve(doAssert="fluxTest_jSurfCyto",\
		       mode="2D_noSSL",
                       params=params,reactions=None,buffers=False,
                       existsCleft=False,existsSSL=False) 
   dotest(concsFinal[idxCa],params.cCaInit,name="fluxTest_jSurfCytoA/F") 

  ## Test scaling of fluxes for cleft 
  # goal here is to find a 10 ms flux that raises the -entire-
  # cell concentration by 1.0 uM
  # 1) We first pump up the cleft for 10 ms, 
  # 2) We use the elevated CaCleft as an init conc, then show the total conc is correct in all compartments 
  ## 1) 
  params = Params()
  params.T = 10
  params.dt = 1 # ms 
  params.D_SSLCyto = 0. 
  params.D_CleftSSL = 0. 
  params.D_CleftCyto = 0. 
  chgConc = 0.5 # [uM] 
  params.jTest = chgConc/(params.dt*params.T) # [uM/ms]
  concsFinal= tsolve(doAssert="fluxTest_jCleft",\
                       params=params,reactions=None,buffers=False)  
  
  caCleftInit = concsFinal[idxCaCleft]

  ## 2) 
  params.T = 100
  params.dt = 10
  params.cInits[idxCaCleft] = caCleftInit
  #params.jTest = 0.
  params.D_SSLCyto = 1
  params.D_CleftSSL = 1
  params.D_CleftCyto = 1
  concsFinal= tsolve(doAssert="conservation",\
                       params=params,reactions=None,buffers=False)  
 
  # check that the change in conc. we anticipated is reflected
  # in cytosol
  assert(abs((params.cCaInit + chgConc) - concsFinal[idxCa]) < 1e-4)
  print "Passed cleft rescalign test" 

# show that in the presence of rapid diffusion, the two and three compartment 
# models are equivalent
def validationRapidDiffusion():
  # Check that compartments equal one another quickly given brief flux 
  if 1: 
    params = Params()
    # NOTE: two cases won't exactly agree, since I havent decided on the best
    # way to merge the 'deleted' ssl volume into the cytosol  
    reactions = "ryrOnlySwitch"
    params.ryrAmp = 1.
    # NOTE: Couldn't avoid having the cleft present a much higher concentration in 
    # the -ssl mode relative to +ssl, since D could not be increased above 1e3
    params.T = 75
    params.dt = 1
    params.D_SSLCyto = 1e2 # Can't go much faster than this   
    params.D_CleftSSL= 1e2  
    params.D_CleftCyto = 1e2  
    params.dist = 0.15
    params.Btot = 0.
    params.Ftot = 0.
    params.cInits[1]=0.
    params.cInits[2]=0.

    idxCa = 0
    mode = "2D_SSL"
    threeComps = tsolve(mode=mode,params=params,hdfName=mode+"_rapid.h5",\
    		reactions = reactions,buffers=False)

    mode = "2D_noSSL"
    twoComps = tsolve(mode=mode,params=params,hdfName=mode+"_rapid.h5",
    			reactions = reactions,buffers=False)
    msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
    assert(abs(threeComps[idxCa] - twoComps[idxCa]) < 1e-4), msg
    print "Passes compartment compare"
  # Check that compartments equal one another quickly given brief flux 
  if 1: 
    params = Params()
    # NOTE: two cases won't exactly agree, since I havent decided on the best
    # way to merge the 'deleted' ssl volume into the cytosol  
    reactions = "ryrOnlySwitch"
    params.ryrAmp = 1.
    # NOTE: Couldn't avoid having the cleft present a much higher concentration in 
    # the -ssl mode relative to +ssl, since D could not be increased above 1e3
    params.T = 150
    params.dt = 2.5
    red = 1e-2
    params.D_SSLCyto = red # Can't go much faster than this   
    params.D_CleftSSL= red  
    params.D_CleftCyto = red  
    params.dist = 0.15
    params.Btot = 0.
    params.Ftot = 0.
    params.cInits[1]=0.
    params.cInits[2]=0.

    idxCa = 0
    mode = "2D_SSL"
    threeComps = tsolve(mode=mode,params=params,hdfName=mode+"_rapid.h5",\
    		reactions = reactions,buffers=False)

    mode = "2D_noSSL"
    twoComps = tsolve(mode=mode,params=params,hdfName=mode+"_rapid.h5",
    			reactions = reactions,buffers=False)
    msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
    assert(abs(threeComps[idxCa] - twoComps[idxCa]) < 1e-4), msg
    print "Passes compartment compare"


def validation():
  #quit()
  validationMergingSSLCyto()
  #raise RuntimeError("NOT FINISHED VALID") 

  ## validate fast/slow diffusion
  validationRapidDiffusion()

  ## flux conversions 
  validation_Conversions()


  ## conservation 
  params = Params()
  idxCa=0
  params.cInits[idxCa]=0.4
  tsolve(doAssert="conservation",reactions=None,buffers=False,params=params)
  print "Passed conservation test"

  ## Reactions
  params = Params()
  params.dt =10
  params.T = 100 
  params.D_SSLCyto = 0
  params.D_CleftSSL = 0 
  params.D_CleftCyto = 0
  # add Ca normally buffered by Buff
  freeCa = 0.1 
  buffedCa = 6.36363636364
  idxCa = 0
  idxCaBuff = 1
  idxCaFluo = 2
  params.cInits[idxCa] = freeCa + buffedCa
  # kill Fluo
  params.Ftot = 0.
  params.cInits[idxCaFluo] = 0. 
  # reset Buff 
  params.cInits[idxCaBuff] = 0. 
  concsFinal=tsolve(doAssert="conservation",params=params, reactions=None,buffers = True)
  eps = 1e-4
  assert(abs(concsFinal[idxCa] - freeCa) < eps) 
  assert(abs(concsFinal[idxCaBuff] - buffedCa) < eps) 
  print "Passed buffering test"


  print "ALL ROUTINES PASSED!"
  

validation()
