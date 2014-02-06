# #
# Template for interfacing gotran models 
# Example combining ode model with PDE simulation 

import sys
sys.path.append("/home/huskeypm/localTemp/wholecell")
import runner
import sympy
import numpy as np
import matplotlib.pylab as plt
import numpy as np
from dolfin import * 

namespace = np.__dict__.copy()
TF = 50; dT = 5  
TF = 2; dT = 0.05 
ms_to_s = 1e-3
s_to_ms = 1e3 
mM_to_uM = 1e3 
Cai_idx = 0
CaTnC_idx = 1

reducedDiff=False
debug=False
#ckMode = "mitoonly"

##
## Standard units: 
## [uM], [ms], [um]
## 


##
## ODE model 
##

import vanbeek_model_2007 as model
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
class empty:pass

## contains all ODE model stuff 
# t0 [ms] for compatibility
# tF [ms] for compatibility
class ODEModel():
  def __init__(self,t0=0,tF=10,vCK=0,ckMode="mitoonly"):
    # Look these up w param_indices func
    self.Caicyt_idx = 0
    self.CaTnCcyt_idx = 1
 
    self.dtn=0.01; # [s] for finding steady state
    self.incr = 11;# for propagating 
    self.t0 = t0 * ms_to_s 
    self.tss = tF * ms_to_s 

    #if(debug):
    #  self.t0=0; self.tss=0.1;  # [s] 

  # sets up ODE model and runs for 10 sec (for steady state) 
  def setup(self):
    ## NOTE: standard time units are [s] here 

    # initialize 
    initvals=model.init_values()
    params=model.default_parameters()
    self.params = params 

    tstepsSS = np.linspace(self.t0, self.tss, (self.tss-self.t0)/self.dtn+1)
  
  
  
    ## get steady state first 
    statesSS= odeint(model.rhs,\
      initvals,tstepsSS,(params,))
    finalStateSS = statesSS[-1::,]

    # init w prev
    self.t0 = self.tss
    finalStateSS= np.ndarray.flatten(finalStateSS)
    self.statesPrev = finalStateSS
    #print finalStateSS

    # return relevant subset 
    statesSS_subset = statesSS[ :,[self.Caicyt_idx, self.CaTnCcyt_idx ]]
    return(tstepsSS,statesSS_subset)

  # advance ODE model by dT from t0 for a given state distribution
  # (statesPrev) 
  # dTpde: time step in PDE loop
  def propagateStates(self,t0_ms,tf_ms,dTpde_ms=1.,incr=-1):
    t0 = t0_ms * ms_to_s
    tf = tf_ms * ms_to_s
    
    if(incr<0):
      incr =self.incr

    #incr =  np.int(np.round((tf-t0)/np.float(self.dtn)))+1
    tsteps = np.linspace(t0, tf, incr)
    initvalsi = np.ndarray.flatten(self.statesPrev)
    print "ODE Init", initvalsi[ [self.Caicyt_idx, self.CaTnCcyt_idx] ]
    print "t0 %f tf %f " % (t0,tf)

    #if(mode=="stateFlux"):
    #if(1):
    #  statesi= odeint(model.rhs,initvalsi,tsteps,(self.params,))
      # plot next
      #pEst,=plt.plot(tsteps,statesi[:,self.Caicyt_idx],'k--')

      # get init/final states 
    #  states0    = np.ndarray.flatten(statesi[0,]) # matched 'initvalsi'
    #  statesF    = np.ndarray.flatten(statesi[-1::,])

      # now use ode's estimate for flux to update concentration manually 
    #  dStatedt = (statesF-states0) / dTpde_ms # [uM/ms]
    #  print "PHASE OUT THIS MODE SINCE ALREADY DONE IN SEP. FLUX"

    #elif(mode=="separate"): 
    if(1): 
      tsteps = np.linspace(t0, tf, incr)
      (statesi,fluxes) = runner.separateFluxes(\
        model,initvalsi,self.params,tsteps,warn=True)
      statesF    = np.ndarray.flatten(statesi[-1::,])
      
      # store and rescale uM/s as uM/ms
      self.dcdt_ATP = fluxes.dcdt_ATP * (1/s_to_ms) 
      self.dcdt_ADP = fluxes.dcdt_ADP * (1/s_to_ms) 
      dStatedt = np.array([self.dcdt_ATP,self.dcdt_ADP]) 
      self.jVolATPase = np.array([fluxes.jVolATPaseATP,fluxes.jVolATPaseADP]) * (1/s_to_ms)
      self.jVolCK = np.array([fluxes.jVolCKATP,fluxes.jVolCKADP]) * (1/s_to_ms)
      self.jBoundary = np.array([fluxes.jBoundATP,fluxes.jBoundADP]) * (1/s_to_ms)

      self.statesPrev = statesF

     
    # return relevant subset 
    statesF_subset = statesF[[self.Caicyt_idx, self.CaTnCcyt_idx ]]
    dStatedt_subset = dStatedt[[self.Caicyt_idx, self.CaTnCcyt_idx ]]
    statesTraj_subset = statesi[:,[self.Caicyt_idx, self.CaTnCcyt_idx ]]
    tsteps_ms = tsteps*s_to_ms



    return (statesF_subset,dStatedt_subset,tsteps_ms,statesTraj_subset)

  # store values into steady state 
  def updateStates(self,cATP,cADP):
    old = np.copy(self.statesPrev)
    self.statesPrev[self.Caicyt_idx] = cATP
    self.statesPrev[self.CaTnCcyt_idx] = cADP

    print "ODE/PDE diff: [%f,%f]" % (self.statesPrev[self.Caicyt_idx]-old[self.Caicyt_idx],
                                     self.statesPrev[self.CaTnCcyt_idx]-old[self.CaTnCcyt_idx])

  def getStates(self):
    return (self.statesPrev[self.Caicyt_idx],self.statesPrev[self.CaTnCcyt_idx] )
    
  
  
##
## PLOTTING
##

def PlotSlice(problem,results,filename,mode="ATP",anisotropic=False):
  1 
  # copy from atpSarcomereVanBeek if needed 

def GlobalConc(problem,results):
  # check values
  results.uAvg= np.zeros(2)
  mesh = problem.mesh
  u = results.u
  t = results.t
  for i,ele in enumerate(split(u)):
    tot = assemble(ele*dx,mesh=mesh)
    vol = assemble(Constant(1.)*dx,mesh=mesh)
    results.uAvg[i]=tot/vol
    print "T %f Conc(%d) %f " % (t,i,results.uAvg[i])
  results.c.append(results.uAvg)

  c1 = results.uAvg[Cai_idx]#get last appended value 
  cb1 = results.uAvg[CaTnC_idx]

  return (c1,cb1)

def Report(problem,results,anisotropic=False):

  u = results.u
  t = results.t
  file= results.file
  #file << (u.split()[0], t)
  file << (u,t)
  

  # extract values at certain regions
  ## copy from atpSarcomereVanBeek if needed 
  


##
## PDE stuff
##

immarker=3

class cparams:
  plot=True
  pvd = False
  
  meshdebug = True  
  tag = ""
  # 
  atpasedistromode="uniform"
  ckdistromode="uniform"
  #atpasedistromode="nonuniform"

  # Chose value from smallest spacing (maximal contraction) 
  #16.43 &   [0.57  0.58  0.82] & 2.4 & 0.45 \\ 
  #  (Dz,Dx,Dy)
  Disotropic = np.ones(3)
  Danisotropic = np.array([0.57, 0.58, 0.82])

  # ATP 
  Datp=0.145 # [um^2/ms]
  #usevanbeek tATP   = 4.88e3 # [uM] Michailova

  # Adp 
  Dadp=0.145 # [um^2/ms]
  #c0adp=0.067e3 # [uM]
  #usevanbeek tADP = 0.113e3 # [uM] Michailova

  # for half-sacomere mesh 
  xMax= 1.   # longitudinal (along sarco)
  yMax = 1.2 # circumferential 
  #zMax = 100.  # transverse (along TT) 
  zMax = 7.5    # estiamte assuming cyldrical bundles are 1.2 um in radius
                 # and h=2 um
                 # Vcyl = pi*r*r*h = xMax*yMax*zMax

  xMiddle = xMax/2.


# NOTE: I am leaving this as a two-component PDE, to make it easier to generalize later
# However, I am only concerned w what's happening in the first component
#stateFlux = 0. # [uM/ms]
def Jboundary(stateFlux=0,vol_sa_ratio=1.):
  v=vol_sa_ratio * stateFlux
  return v

# left hand side of box has Ca2+, Ca2+/Buff present 
class InitialConditions(Expression):
  def eval(self, values, x):
    values[0] = self.c0i
    values[1] = self.c1i
  def value_shape(self):
    return (2,)

# Class for interfacing with the Newton solver
class MyEqn(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity)
        self.reset_sparsity = False

##class Extracellular(SubDomain):
#class IMSpace(SubDomain):
#  def inside(self,x,on_boundary):
#    return x[0] < DOLFIN_EPS  and on_boundary # or x[0] > 1.0 - DOLFIN_EPS
#    #return on_boundary
class IMSpace(SubDomain): 
  def inside(self,x,on_boundary):
    result = x[1]>(self.params.yMax-DOLFIN_EPS) and on_boundary
    #print x,result 
    return result


##
## Functions
## 


## propogates a concentration variable in a PDE according to an ODE model 
# mode="totalFlux"    # fluxes are computed based on concentration change 
# mode="separateFlux" # fluxes are separated into boundary and volume fluxes
#
# NOTE: I think I will leave this as a mixed species PDE, so that I can generalize to Ca + mobile buffers
#
def runPDE(\
  params,
  mode="totalFlux", 
  anisotropic=True,
  duration=1e3,		# duration [ms] 
  asserts=False): 

  ## params 
  t0ode = 0
  tFode= 5e3 # [ms]
  t0pde=tFode
  if(debug):
    duration=500
  tFpde=t0pde + duration
  dt = 10. # [ms]

  
  ##
  ## ODE
  ##
  ## get exact solution 
  odeModel = ODEModel(t0=t0ode,tF=tFode,ckMode=params.ckMode)
  odeModel.setup() 
  #print odeModel.getStates()
  # here we propagate model just so that we can compare exact and PDE solutions 
  # but we'll need to 'reset' the statesPrev variable
  statesPrev = np.copy(odeModel.statesPrev)
  (dummy,dummy,tstepsExact,statesTrajExact)= odeModel.propagateStates(\
    t0pde,tFpde,incr=1e4)
  odeModel.statesPrev = statesPrev
  #c0i = odeModel.statesPrev[self.Caicyt_idx]
  #c1i = odeModel.statesPrev[self.CaTnCcyt_idx]
  c0i,cb0i = odeModel.getStates()
  #print c0i,cb0i
  #print odeModel.statesPrev
  
  ##
  ## PDE 
  ##

  ## setup PDE problem 
  if(params.meshdebug):
   # mesh = UnitCube(2,2,2)
    mesh = UnitCube(10,12,10)
  else:
    mesh = UnitCube(25,25,10)
    #mesh = UnitCube(25,25,10)
  # rescale mesh to size of sarcomere 
  mesh.coordinates()[:]= np.array([params.xMax,params.yMax,params.zMax])* mesh.coordinates()

  V = FunctionSpace(mesh,"CG",1)
  ME = V *V

  # Trial and Test functions 
  du = TrialFunction(ME)
  q,v  = TestFunctions(ME)

  # Define function
  u = Function(ME) # current soln
  u0 = Function(ME) # prev soln

  # split mixed functions
  c,cb = split(u)
  c0,cb0 = split(u0)

  # mappings
  indcAtp, indcAdp= set(), set()
  dm_cAtp, dm_cAdp = ME.sub(0).dofmap(), ME.sub(1).dofmap()
  for cell in cells(mesh):
      indcAtp.update(dm_cAtp.cell_dofs(cell.index()))
      indcAdp.update(dm_cAdp.cell_dofs(cell.index()))
  indcAtp = np.fromiter(indcAtp, dtype=np.uintc)
  indcAdp = np.fromiter(indcAdp, dtype=np.uintc)

  cAtp_values = u.vector()[indcAtp]
  cAdp_values = u.vector()[indcAdp]

  # Init conts
  init_cond = InitialConditions()
  init_cond.c0i= c0i
  init_cond.c1i= cb0i
  u.interpolate(init_cond)
  u0.interpolate(init_cond)


  ## expr
  # JOHAN - doesn't work?
  #t = 1
  #print "Todo: Why does this not wokr/"
  #jsyn_expr = Jsyn()
  #jsyn_expr_str = str(jsyn_expr)
  #namespace["t"] = t
  #print eval(jsyn_expr_str,namespace)
  # replacement 
  jboundaryExpr = Expression("j",j=0)
  jvolATPaseExpr = Expression("j",j=0)
  jvolCKExpr = Expression("j",j=0)

  ## Decide on diff const 
  Dii = params.Disotropic

  # homog of myfilaments has 'z' as the long axis, but in half-sarcomere, 'x' is the long axis 
  Dij = Constant([Dii[2],Dii[0],Dii[1]])
  Dij = diag(Dij)
  Dij_ATP= Constant(params.Datp) * Dij
  Dij_ADP = Constant(params.Dadp) * Dij


  ## ODE weak form 
  RHS1 = -inner(Dij_ATP*grad(c),grad(q))*dx
  RHS2 = -inner(Dij_ADP*grad(cb),grad(v))*dx

  # adjust source term depending on how ATPase is distributed
  # This is a work-around (see below)
  # ATPase 
  if(params.atpasedistromode=="uniform"):
    RHS1 +=  jvolATPaseExpr*q*dx
    RHS2 += -jvolATPaseExpr*v*dx
  elif(params.atpasedistromode=="nonuniform"):
    # NOTE: the nonuniform results (even for scale=1) do 
    # not validate exactly. Don't understand why. 130706 
    # Hence, we don't validate this part for unit testing 
    # Compute ATPase density 
    gDensity = ATPaseDensityInit(ME,mesh,params)
    scale = Function(V)
    scale.vector()[:]=gDensity
    #scale.vector()[:]=1.0      
    RHS1 +=  scale*jvolATPaseExpr*q*dx
    RHS2 += -scale*jvolATPaseExpr*v*dx
  # CK 
  if(params.ckdistromode=="uniform"):
    RHS1 +=  jvolCKExpr*q*dx
    RHS2 += -jvolCKExpr*v*dx
  elif(params.ckdistromode=="nonuniform"):
    # NOTE: the nonuniform results (even for scale=1) do 
    # not validate exactly. Don't understand why. 130706 
    # Hence, we don't validate this part for unit testing 
    # Compute ATPase density 
    gDensity = CKDensityInit(ME,mesh,params)
    scale = Function(V)
    scale.vector()[:]=gDensity
    #scale.vector()[:]=1.0      
    RHS1 +=  scale*jvolCKExpr*q*dx
    RHS2 += -scale*jvolCKExpr*v*dx


  ## boundary
  subdomains = MeshFunction("uint",mesh,2)
  boundary = IMSpace()
  boundary.params = params 
  boundary.mark(subdomains,immarker)
  ds = Measure("ds")[subdomains]

  # vol/sa ratio to rescale surf. fluxes
  vol = assemble(Constant(1.)*dx,mesh=mesh)
  sa  = assemble(Constant(1.)*ds(immarker),mesh=mesh)
  vol_sa_ratio = vol/sa
  #print "vol_sa_ratio  [um]",  vol_sa_ratio

  # define flux terms 
  print "WARNING: not working w eval func"
  jL1 =  jboundaryExpr
  jL2 =  -jboundaryExpr # this is not generally correct, since ADP/ATP fluxes may varyn
  RHS1  += jL1*q*ds(immarker)
  RHS2  += jL2*v*ds(immarker)


  # Add in time dependence 
  # (dc/dt) = RHS  --> c1-c0 - dt * RHS = 0
  L1 = c*q*dx - c0*q*dx - dt * RHS1
  L2 = cb*v*dx - cb0*v*dx - dt * RHS2
  L = L1 + L2

  # Compute directional derivative about u in the direction of du (Jacobian)
  # (for Newton iterations) 
  a = derivative(L, u, du)


  # Create nonlinear problem and Newton solver
  # not really needed, but will keep here to generalize later 
  problem = MyEqn(a, L)
  solver = NewtonSolver("lu")
  solver.parameters["convergence_criterion"] = "incremental"
  solver.parameters["relative_tolerance"] = 1e-6

  # Output file
  file = File("output.pvd", "compressed")

  # Step in time
  t = 0.0
  results = empty()
#  problem = empty()

  # Initializations for results 
  # See atpSarc...py
  #results.x_midsarc_aband = np.array([0.50,0.50,0.5])
  results.ts = []
  results.adpGradients = []
  #results.jATP=[]
  #results.jANT=[]
  results.c = []
  problem.mesh=mesh
  problem.u = u
  problem.params = params
#
  results.ME = ME
  results.u = u
  results.file = file


  # initialize 
  c0F = c0i

  ## Time stepping
  jboundarys=[]
  jvolATPases=[]
  jvolCKs=[]
  vals=[[],[]]
  cs=[]
  cbs=[]
  cFs=[]
  ts = np.linspace(t0pde+dt,tFpde,(tFpde-t0pde)/dt)
  #print t0, tFp
  #print ts
  for i,tf in enumerate(ts):
    print "###########################################################################"
    ## update ODE 
    # predict 'forward' solution from ode model '
    #tf = t0+dt
    t0 = tf - dt
    (statesF,stateFlux,dummy,dummy) = odeModel.propagateStates(t0,tf,dt)
    print "ODE states: ", statesF
    print "ODE statesFlux: ", stateFlux

  
    ## update PDE 
    # define flux 
    if(mode=="totalFlux"): 
      atpFlux = stateFlux[Cai_idx] # s.b. same as adp flux
      jboundaryExpr.j =Jboundary(atpFlux,vol_sa_ratio=vol_sa_ratio)

    elif(mode=="separateFlux"): 
      # It appears that ATP/ADP are equal and opposite, so we can make this assumption here 
      # apply boundary fluxes to surface term 
      jboundaryExpr.j =Jboundary(odeModel.jBoundary[Cai_idx],vol_sa_ratio=vol_sa_ratio)
      jvolATPaseExpr.j =odeModel.jVolATPase[Cai_idx]
      jvolCKExpr.j =odeModel.jVolCK[Cai_idx]

      # tests to make sure implemented correctly 
      #print "WARNING: hack " 
      #print stateFlux
      #print odeModel.dcdt_ATP
      #print odeModel.jBoundary[Cai_idx]
      #print odeModel.jVolATPase[Cai_idx]
      ## try all as boundary (works) 
      #atpFlux = odeModel.jBoundary[Cai_idx] + odeModel.jVolATPase[Cai_idx]
      #jboundaryExpr.j =Jboundary(atpFlux,vol_sa_ratio=vol_sa_ratio)
      #jvolATPaseExpr.j = 0.
      ## try all as volume 
      #jboundaryExpr.j =0.
      #jvolATPaseExpr.j = atpFlux
      
      

    # get integrated fluxes 
    totjs = assemble(jboundaryExpr*ds(immarker),mesh=mesh)
    jboundary = totjs/sa * 1/vol_sa_ratio
    totjs = assemble(jvolATPaseExpr*dx,mesh=mesh)
    jvolATPase = totjs/vol
    totjs = assemble(jvolCKExpr*dx,mesh=mesh)
    jvolCK = totjs/vol


    # solve system
    u0.vector()[:] = u.vector()
    solver.solve(problem, u.vector())
    if(params.pvd):
      file << (u,tf)

    # report/store values
    # not very general, but oh well
    ## Collect results 
    results.t = tf
    (c1,cb1)=GlobalConc(problem,results)
    #print "REMOVE THIS SECTION???"
    #pdeVals = []
    #for i,ele in enumerate(split(u)):
    #  tot = assemble(ele*dx,mesh=mesh)
    #  conc = tot/vol
    #  print "Conc(t=%3.2f,%d) %f " % (tf,i,conc)
    #  pdeVals.append(conc)
    #  vals[i].append(conc)

    #print "rename c1,cb1"
    #c1 = pdeVals[Cai_idx]#get last appended value 
    #cb1 = pdeVals[CaTnC_idx]
    #ts.append(t)
    jboundarys.append(jboundary)
    jvolATPases.append(jvolATPase)
    jvolCKs.append(jvolCK)
    results.jdiffs= jboundarys
    results.jATPases = jvolATPases
    results.jCKs = jvolCKs         
    cFs.append(statesF[Cai_idx]) # just want 'ATP' from odeModel
    #cs.append(c1)
    #cbs.append(cb1)

    ## Package results 
    Report(problem,results,anisotropic=anisotropic)


    ## store values in ODE  
    #ts.append(tf)
    t0 = tf
    odeModel.updateStates(c1,cb1)


  ## END LOOP 


#  print "WARNING: merge into plotting"
  # convert lists into arrays
  cFs = np.asarray(cFs)
  cs = np.asarray(results.c)[:,Cai_idx]
  cbs = np.asarray(results.c)[:,CaTnC_idx]
#
#  ## Plot 
  fig=plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.set_title("[ATP]/[ADP] vs time") 
  ax1.plot(tstepsExact,statesTrajExact[:,Cai_idx],label="[ATP] Exact")
  ax1.plot(ts,cs,'k.',label="[ATP] (PDE) ")
  ax1.plot(ts,cFs,label="[ATP]F (forward ODE solutions)")
  ax1.set_ylabel("[ATP] [uM]") 
  ax1.set_xlabel("t [ms]") 
  ax1.legend(loc=0)
  ax2 = ax1.twinx()
  ax2.set_ylabel("[ADP] [uM]") 
  ax2.plot(ts,cbs,'r-.',label="cb (PDE) ")
  plt.gcf().savefig(mode+"conc.png")
#
#
#  plt.figure()
#  jboundarys = np.asarray(jboundarys)
#  jvolATPases = np.asarray(jvolATPases)
#  #plt.plot(ts,jsyns,label="jhyds")
#  plt.plot(ts,jboundarys,label="jboundarys")
#  plt.plot(ts,jvolATPases,label="jvolATPases")
#  plt.legend(loc=0)
#  plt.gcf().savefig(mode+"js.png")

  ext=statesTrajExact[-1,Cai_idx]
  est = c1  # c1 is updated each frame, so this is the most recent
  m = "mode %s: Exact %f != Est %f "%(mode,ext,est)
  if(asserts):
    ext=statesTrajExact[-1,Cai_idx]
    #est=cs[-1]
    est = c1  # c1 is updated each frame, so this is the most recent
    m = "mode %s: Exact %f != Est %f "%(mode,ext,est)
    assert(np.abs(ext-est)<1.0), m
    print "### Mode %s passes!" % mode
  elif (np.abs(ext-est)>1.0):
    print "WARNING: ", m 

  return problem,results 
  

###
### Testing plotting etc 
###
  
# plots isotropic/anisotropic results
def plotting(params,resultsiso,resultsaniso):
  1
  # see atpSarc...py

def testODE():
  t0 = 0.
  tFp = 10e3 # [ms] 
  incr = 100*tFp
  CaTnC_idx = 1

  # vCK = 1
  odeModel = ODEModel(vCK=1.0)
  odeModel.setup()
  (dummy,dummy,tstepsExact,statesTrajExact)= odeModel.propagateStates(t0,tFp,incr=incr) 
  plt.plot(tstepsExact,statesTrajExact[:,CaTnC_idx],label="vCK = 1")       

  # vCK = 1
  odeModel = ODEModel(vCK=0.02)
  odeModel.setup()
  (dummy,dummy,tstepsExact,statesTrajExact)= odeModel.propagateStates(t0,tFp,incr=incr) 
  plt.plot(tstepsExact,statesTrajExact[:,CaTnC_idx],label="vCK = 0")       
  plt.legend(loc=0)
  plt.xlim([9.2e3,10.0e3]) 


  print "WARNING: not a unit test"
  plt.gcf().savefig("testode.png") 


def testODE2():
  t0 = 0.
  tFp = 10e3 # [ms] 
  dt=10
  Cai_idx = 0

  # vCK = 1
  odeModel = ODEModel(vCK=1.0)
  odeModel.setup()
  
  ts = np.linspace(t0+dt,tFp,(tFp-t0)/dt)
  
  jConcs=[]
  jSums=[]
  for i,tf in enumerate(ts): 
    t0 = tf - dt 
    (statesF,stateFlux,dummy,dummy) = odeModel.propagateStates(t0,tf,dt)
    jSum = odeModel.jBoundary + odeModel.jVolATPase + odeModel.jVolCK
    jConcs.append(stateFlux[Cai_idx])
    jSums.append(jSum[Cai_idx])

  plt.plot(jConcs,label="jConc") 
  plt.plot(jSums,label="j individual") 
  plt.legend()
  print "WARNING: not a unit test"
  plt.gcf().savefig("testODE2.png") 

def loop(p,duration=1e3,asserts=True,case="noCK"):
  if(debug):
    duration=10
    

  p.tag = case+"iso"
  problemi,riso = runPDE(params=p,mode="separateFlux",anisotropic=False,asserts=asserts,duration=duration)
 



def Test1(do="all"):
    cparams.plot = False
    pi = cparams()
    loop(pi,duration=1e3,case="fullCK")  

  


import sys
#
# Revisions
#       10.08.10 inception
#

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  params = cparams()

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]

  for i,arg in enumerate(sys.argv):
    if(arg=="-debug"):
      debug = True
      Test1()
    if(arg=="-validation"):
      raise RuntimeError("not yet supported") 
