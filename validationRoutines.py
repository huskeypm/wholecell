"""
Validation routines for sachse.py

NOTE: successful unit checking is don through validation.py
"""
from sachse import * 
def myassert(ref,actual,eps):
    msg = ref , "!=" , actual 
    assert( np.abs(ref-actual) < eps), msg
    print "PASS!"

def validationSERCA():
  # just spot checks on fluxes
  
  import simple 
  mpi=True
  mpi=False
  params = Params()
  if mpi==False:
    jSERCA  = simple.SERCAExpression()
    mesh = UnitIntervalMesh(200)
    V = FunctionSpace(mesh,"CG",1)
    f = Function(V)

    # at cai 
    #mystr = "vPump * pow(cai,m) / (KmPump + pow(cai,m))"
    def expr(cai,Vmax=params.sercaVmax,Kmf=params.sercaKmf,H=params.sercaH, casrTerm=params.sercaCaSRTerm):
       return -Vmax*(-casrTerm +pow(cai/Kmf, H))/(1 + casrTerm + pow(cai/Kmf, H))

    cai = 0.4
 
    jSERCA.cai = cai # uM
    f.interpolate(jSERCA)
    evaled = np.asarray(f.vector())[0]
    #assert( (expr(cai) - evaled ) < 1e-4 ), "SERCA failed"
    myassert( expr(cai) , evaled , 1e-1 )                  


  # verify that constant flux ok rthrough SERCA pexression  
  # m=0 
  # Km = 0 
  # Vma = 0.1 
  #########
  idxCa = 0
  params.T = 10.        
  params.dt = 1.0 # ms
  params.ryrOffset = 500
  params.sercaVmax = 0.1 # uM/ms
  params.sercaKmf = 0.000 # uM
  params.sercaH = 0.  # Hill coeff 
  params.sercaCaSRTerm = 0. # cancel out backflux 
  # for these params, jSERCA = -0.1 [uM/ms]
  jSERCA = params.sercaVmax / 2. # (since Expression is Vmax/(km^0+ca^0)=V/2.
  dCa = params.T  * jSERCA  
  

  params.cInits[idxCa] = 2.0
  refFinal = params.cInits[idxCa] + dCa
  params.D_SSLCyto = 0.
  params.D_CleftSSL = 0.
  params.D_CleftCyto = 0.
  reactions = "simple" 
  concsFinal= tsolve(mode="2D_noSSL",
                     params=params,reactions=reactions,buffers=False, 
                     existsCleft=False,existsSSL=False) 

  assert( (concsFinal[idxCa] - refFinal) < 1e-4 ),  "SERCA test failed" 
  print "PASS SERCA test"

# Consistency check 
def validationTestGeom(): 
  params = Params()
  idxCaCleft = 4
  idxCaSSL   = 3
  idxCa = 0
  idxCaBuff = 1
  idxCaFluo = 2

  params.ryrOffset = 10.

  params.D_SSLCyto = 100.
  params.D_CleftSSL = 1e-1
  #params.D_CleftCyto = 0.

  withBuffers = True
  if withBuffers:
    refFinal = 0.07233 # PKH 150501
  else: 
    params.sercaVmax = 5e-3
    refFinal = 0.0832763346023#PKH 150501

  params.T = 100 
  #params.T = 20  
  params.dt = 0.25 # ms
  
  testNum = 1
  if testNum==1: 
    # rapid diffusion / unstable  
    params.D_SSLCyto = 1e3
    params.D_CleftSSL = 1e2
    reactions = "ryrOnly"
  elif testNum==2: 
    # slow diffusion 
    params.D_SSLCyto = 10.
    params.D_CleftSSL = 1e-3
    reactions = "simple"

  else: 
    reactions = "simple"
  tag = "2D_SSL"
  concsFinal= tsolve(mode=tag,hdfName="validation"+tag+".h5",
                       params=params,reactions=reactions,buffers=withBuffers)
  
  myassert( concsFinal[idxCa] , refFinal ,  1e-4 )
  print "PASS SERCA test"



# Test that addition of Ca2+ via RyR is approximately offset by SERCA 
# uptake 
#def validationStepwise(withBuffers=False): 
def validationStepwise(withBuffers=False,ryrDebug=False): 
  params = Params()
  idxCaCleft = 4
  idxCaSSL   = 3
  idxCa = 0
  idxCaBuff = 1
  idxCaFluo = 2

  # ryr stuff 
  si=0; fi=300; st = fi
  dt = 1.0 
  ts = np.linspace(si,fi,(fi-si)/dt+1)  
  ts = ts[1:]
  params.ryrOffset = 0.

  # First option tests for constant RyR flux 
  # Second uses usual exponential decay 
  if ryrDebug:
    reactions="ryrOnlySwitch"
    params.ryrAmp = 0.1
    params.ryrTau = 1e9
    params.ryrTerminate = 300.
    i_s = params.ryrAmp*np.ones(np.shape(ts))
    eps = 1e-2  
  else:
    reactions="ryrOnly"
    i_s = params.ryrAmp*np.exp(-ts/params.ryrTau)
    eps = 1  # generous b.c. of use of dirac function and 




  # For simple numerical estimate of total added Ca2+ due to RyR
  params.T = fi - si
  js = wholeCellI_to_wholeCellJ(i_s)
  dt = ts[1]-ts[0]
  dCas = js*dt
  totDCa = np.cumsum(dCas)
  finaldCa  = totDCa[-1]
  finaldCa+= params.cInits[idxCa] # add in init conc 
  print "Tot expected ", finaldCa
  #quit()

  step1=False 
  step2=False 
  step3Test=False 
  step3 = True  
  #########
  ## 1) release RyR Ca2+ into partitioned cleft 
  if step1:
    params.dt = dt # ms
  
    params.D_SSLCyto = 0.
    params.D_CleftSSL = 0.
    params.D_CleftCyto = 0.
  
    concsFinal= tsolve(
                         params=params,reactions=reactions,buffers=withBuffers)
    caCleftInit = concsFinal[idxCaCleft]
  else: 
    caCleftInit = 1.24556843e+04 # PKH 
  

  if step2:
    ## 2) Open cleft and release Ca2+ into cyto 
    params.T = 40  
    params.dt = 10 # smaller steps will usually fail for nobuffer case 
    params.cInits[idxCaCleft] = caCleftInit 
    #params.jTest = 0.
    params.D_SSLCyto = 1
    params.D_CleftSSL = 1
    params.D_CleftCyto = 1
    concsFinal= tsolve(doAssert="conservation",\
                         params=params,reactions=None,buffers=withBuffers)
    finalConc = concsFinal[idxCa] 
    print "concsFinal idx: ",concsFinal[idxCa] #CES
    if withBuffers==False:
      myassert(finaldCa,concsFinal[idxCa],eps)
      print "Tot expected ", finaldCa

  else:
      if withBuffers:
        finalConc = 0.72071# PKH est 
      else:
        finalConc = 41.3290105182
      
      concsFinal= np.array([finalConc,0,0,finalConc,finalConc])

  # verify that constant flux ok rthrough SERCA pexression  
  #########
  ## 3) suck up cyto Ca via SERCA  
  print "CAITLIN - add a flag for me" 

#  ### Constant rate, but WITHOUT using RyR/SERCA expressions  
#  if step3:
#    params = Params()
#    params.T = 90 # ms 
#    params.dt = 10.0 # ms
#    params.ryrOffset = 500
#    # for these params, jSERCA = -0.2 [uM/ms]
#    jSERCA = -0.2  # [uM/ms]
#    params.sercaVmax = jSERCA
#    dCa = params.T  * jSERCA
#  
#  #  params.cInits[idxCa] = 20.1 #CES
#    params.cInits[idxCa] = concsFinal[idxCa] #CES
#    refFinal = params.cInits[idxCa] + dCa
#    params.D_SSLCyto = 0.
#    params.D_CleftSSL = 0.
#    params.D_CleftCyto = 0.
#
#    reactions = "caitlinSERCA"
#    concsFinal= tsolve(mode="2D_noSSL",
#                     params=params,reactions=reactions,buffers=withBuffers,
#                     existsCleft=False,existsSSL=False)
#
#    assert( (concsFinal[idxCa] - refFinal) < 1e-4 ),  "SERCA test failed"
 
  ### Constant SERCA rate, but through using expressions in Simple
  # Equivalent to a constant flux for serca, since we remove cai dependence 
  if step3Test:
    params = Params()
    idxCa = 0
    params.T = 100
    params.dt = 10.0 # ms
    params.ryrOffset = 1500 # delay firing of Ryr until long after the sim
    params.sercaVmax = 0.4 # uM/ms
    params.sercaKmf = 0.000 # uM
    params.sercaH = 0.  # Hill coeff 
    params.sercaCaSRTerm = 0. # cancel out backflux 
    # for these params, jSERCA = -0.1 [uM/ms]
    jSERCA = -params.sercaVmax / 2. # (since Expression is Vmax/(km^0+ca^0)=V/2.
    dCa = params.T  * jSERCA
  
  
    params.cInits[idxCa] = finalConc                  
    refFinal = params.cInits[idxCa] + dCa
    params.D_SSLCyto = 0.
    params.D_CleftSSL = 0.
    params.D_CleftCyto = 0.
    reactions = "simple"
    concsFinal= tsolve(mode="2D_noSSL",
                       params=params,reactions=reactions,buffers=False,
                       existsCleft=False,existsSSL=False)
  
    jSERCAAvg = (finalConc - concsFinal[idxCa])/params.T
    print "Avg jSERCA ", jSERCAAvg

    if withBuffers==False:
      #assert( (concsFinal[idxCa] - refFinal) < 1e-4 ),  "SERCA test failed"
      myassert(concsFinal[idxCa],refFinal, 1e-4)

  ### Uses usual, Ca-dependent SERCA model 
  # Uses a serca model that will be nonzero until the cellular conc is 0. 
  # therefore, we stop this t some time T s.t. conc is reasonable.
  # Note too that dt=1.; otherwise SERCA overshoots cai=0
  params = Params()
  if withBuffers:
    refFinal = 0.07233 # PKH 150501
  else: 
    params.sercaVmax = 5e-1
    refFinal = 0.07235294124753 #PKH 150501

  if step3:
    idxCa = 0
    params.T = 200
    params.dt = 1.0 # ms
    params.ryrOffset = 1500
  
    params.cInits[idxCa] = finalConc 
    #refFinal = params.cInits[idxCa] + dCa
    params.D_SSLCyto = 0.
    params.D_CleftSSL = 0.
    params.D_CleftCyto = 0.
    reactions = "simple"
    concsFinal= tsolve(mode="2D_noSSL",
                       params=params,reactions=reactions,buffers=withBuffers,
                       existsCleft=False,existsSSL=False)
  
    myassert( concsFinal[idxCa] , refFinal ,  1e-4 )
    print "PASS SERCA test"

  print "PASSED STEPWISE with buffer=",withBuffers

def validationRyR(reactions="ryrOnly"):

  # this expression is taken from simple.py
  #self.iryr = Expression(dirac+"a*exp(-(t-to)/tau)",
  #                           v = 0.05,\
  #                           a=params.ryrAmp,\
  #                           to=params.ryrOffset,\
  #                           tau=params.ryrTau,\
  #                           t=0)

  params = Params()
  idxCaCleft = 4
  idxCa = 0

  if reactions == "ryrOnlySwitch":
    si=0; fi=10; st = fi
    ts = np.linspace(si,fi,st+1)
    ts = ts[1:]
    #params.ryrAmp = 1.
    iFlux = params.ryrAmp # A/F
    i_s = np.ones(np.shape(ts)[0]) * iFlux # A/F
    eps = 1e-2
  if reactions == "ryrOnly":
    si=0; fi=300; st = fi
    ts = np.linspace(si,fi,st+1)
    ts = ts[1:]
    params.ryrOffset = 0.
    i_s = params.ryrAmp*np.exp(-ts/params.ryrTau)
    eps = 1e-1   # generous b.c. of use of dirac function and 

  params.T = fi - si
  js = wholeCellI_to_wholeCellJ(i_s)
  dt = ts[1]-ts[0]
  dCas = js*dt
  totDCa = np.cumsum(dCas)
  finaldCa  = totDCa[-1]
  #print "Tot expected ", finaldCa

  #########
  params.dt = dt # ms

  params.D_SSLCyto = 0.
  params.D_CleftSSL = 0.
  params.D_CleftCyto = 0.

  ## test using 'jTest'
  # checked that jtest/ryrOnlySwitch agree PKH
  #PKHparams.jTest = wholeCellI_to_wholeCellJ( iFlux ) # 1 A/F --> um/ms 
  #PKHconcsFinal= tsolve(doAssert="fluxTest_jCleft",\
  #PKH                     params=params,reactions=None,buffers=False)
  #PKHcaCleftInit = concsFinal[idxCaCleft]
  #quit()

  concsFinal= tsolve(
                       params=params,reactions=reactions,buffers=False)
  caCleftInit = concsFinal[idxCaCleft]


  ## 2) 
  params.T = 30  
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
  tot = params.cCaInit + finaldCa
  refFin = concsFinal[idxCa]
  msg = "%f/%f > %f " %(tot,refFin, eps) 
  assert(abs((tot-refFin) / refFin) < eps),msg

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
    threeComps = tsolve(mode=mode,params=params,hdfName="sacshe2TT_%s_%s.h5"%(mode,case),
      reactions = reactions,buffers=buffers) 
 
  # merged compartments 
  #print "###\n###\n####\n"
  mode = "2D_noSSL"
  params.T = 76     
  params.dt = 0.5
  twoComps = tsolve(pvdName="test.pvd",mode=mode,params=params,hdfName="sacshe2TT_%s_%s.h5"%(mode,case),
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

# Master validation function 
def validation(test=1):
  test = int(test)
  
  if test==1:
    validationRyR()
  #raise RuntimeError("NOT FINISHED VALID") 
  if test==13: 
    validationCaitlinSERCA()

  if test==14: 
    validationCaitlinSERCA(withBuffers=True)

  #quit()
  if test==2:
    validationMergingSSLCyto()

  ## validate fast/slow diffusion
  if test==3:
    validationRapidDiffusion()

  ## flux conversions 
  if test==4:
    validation_Conversions()

  if test==5:
    validationsMisc()

  if test==6:
    validationSERCA()


def validationsMisc():
  ## conservation 
  params = Params()
  idxCa=0
  params.cInits[idxCa]=0.4
  tsolve(doAssert="conservation",reactions=None,buffers=False,params=params)
  print "Passed conservation test"

  ## Reactions
def validationBuffering():
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


  

if __name__ == "__main__":
  import sys
  validation(sys.argv[1])
