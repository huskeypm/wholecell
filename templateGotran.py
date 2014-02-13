##
## Generic template for pde models using gotran ode models 
##
## Need to edit all 'TEMPLATE' areas 
##
## Standard units: 
## [mM], [ms], [um]
## -- if units are changed, make sure changes are reflected EVERYWHERE in code 
## 

## TODO buffering by TnC in cytosol is not correct in this model!!!
print'WARNING: buffering by TnC in cytosol is not correct in this model!!!'
# I think I should return Ca + CaTnC into cyto 
## TODO separate buffering in vol from SR flux 

## TODO can't update Ca concentrations, for some reason, without breaking ODE  
## TODO use real geomrtry 
## TODO modify Kmf to agree w guess


import sys
sys.path.append("/home/huskeypm/localTemp/wholecell")
import runner
import sympy
from dolfin import * 

# for cmd line operation 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# set var, namespaces
import shannon_2004 as model
# TODO Check
model = runner.model 
param_indices = runner.model.param_indices
state_indices = runner.model.state_indices
namespace = np.__dict__.copy()
TF = 50; dT = 5  
TF = 2; dT = 0.05 
ms_to_s = 1e-3
s_to_ms = 1e3 
mM_to_uM = 1e3 

## TEMPLATE
# define states
stateNames=("Cai","Ca_TroponinC") 
nStates = len(stateNames) 
Cai_pde_idx=0
CaTnC_pde_idx=1

debugLevel=0

# <codecell>


##
## ODE model 
##

from scipy.integrate import odeint
class empty:pass

## contains all ODE model stuff 
# t0 [ms] for compatibility
# tF [ms] for compatibility
# Assumes shannon model for now (I do have hacked one around) 
class ODEModel():
  # NOTE: except for stimPeriod, pass all other params through 'params' 
  def __init__(self,t0=0,tF=2,stimPeriod=300,params=None):
    self.Cai_ode_idx =  state_indices( stateNames[0] ) 
    self.CaTnC_ode_idx =  state_indices( stateNames[1] ) 

    # Only keep this one 
    self.stim_period_idx = param_indices("stim_period")
    self.stimPeriod = stimPeriod # [ms] 
    self.params = params 

    ## TODO remove
    self.kon_tnc_idx = param_indices("kon_TroponinC")  
    self.koff_tnc_idx = param_indices("koff_TroponinC") 
    self.bmax_tnc_idx = param_indices("Bmax_TroponinC")


    self.Cai_pde_idx =  0
    self.CaTnC_pde_idx = 1 
    
    # if PDE/ODE are at different time units 
    self.time_PDE_to_ODE = 1.    # or .... =  ms_to_s
    self.time_ODE_to_PDE = 1/self.time_PDE_to_ODE
 
    ## TODO verify time units 
    self.dtn=1.00; # [ode time units] for finding steady state
    self.incr = 11;# for propagating 
    self.t0 = t0 * self.time_PDE_to_ODE
    self.tss = tF * self.time_PDE_to_ODE

    #if(debug):
    #  self.t0=0; self.tss=10.0;  # [ode time units] 

  # sets up ODE model and runs for 10 sec (for steady state) 
  def setup(self):
    ## NOTE: standard time units are [ode time units] here 

    # initialize 
    initvals=model.init_values()
    if(self.params==None):
      self.params=model.default_parameters()
      #self.params = params 

    #print self.params[param_indices("Kmf")]
    #quit()

    # params 
    self.params[self.stim_period_idx] = self.stimPeriod 

    print "Equilibrating for %f [ms]" % (self.tss - self.t0)

    tstepsSS = np.linspace(self.t0, self.tss, (self.tss-self.t0)/self.dtn+1)
    #print tstepsSS
  
    ## get steady state first 
    ## TODO CHECK THIS STEADY STATE W PLOITS 
    print "WARHING: CHECK ME" 
    statesSS= odeint(model.rhs,\
      initvals,tstepsSS,(self.params,))
    finalStateSS = statesSS[-1::,]

    # init w prev
    self.t0 = self.tss
    finalStateSS= np.ndarray.flatten(finalStateSS)
    self.statesPrev = finalStateSS
    #print finalStateSS

    # return relevant subset 
    statesSS_subset = statesSS[ :,[self.Cai_ode_idx, self.CaTnC_ode_idx ]]

    plt.plot(tstepsSS,statesSS_subset[:,0])
    plt.gcf().savefig("equil.png") 

    return(tstepsSS,statesSS_subset)


  # advance ODE model by dT from t0 for a given state distribution
  # (statesPrev) 
  # dTpde: time step in PDE loop
  def propagateStates(self,t0_ms,tf_ms,dTpde_ms=1.,incr=-1):
    t0 = t0_ms * self.time_PDE_to_ODE
    tf = tf_ms * self.time_PDE_to_ODE

    
    if(incr<0):
      incr =self.incr

    #incr =  np.int(np.round((tf-t0)/np.float(self.dtn)))+1
    tsteps = np.linspace(t0, tf, incr)
    initvalsi = np.ndarray.flatten(self.statesPrev)
    print "ODE Init [mM]", initvalsi[ [self.Cai_ode_idx, self.CaTnC_ode_idx] ]
    print "t0 %f tf %f [ms]" % (t0,tf)

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
      # s.b. generalized more 
      self.dcdt = np.ndarray.flatten(fluxes.dcdt) * (1/self.time_ODE_to_PDE)

      ## BEGIN TEMPLATE 
      # define volumetric fluxes (occur within PDE region) 
      self.jVol= fluxes.jVol * (1/self.time_ODE_to_PDE)
      self.jSR = fluxes.jSR

      # define all boundary fluxes (occur on PDE boundary)  
      self.jBoundary = fluxes.jBoundary * (1/self.time_ODE_to_PDE)
      ## END TEMPLATE 

      self.statesPrev = statesF

     
    # return relevant subset 
    statesF_subset = statesF[[self.Cai_ode_idx, self.CaTnC_ode_idx ]]
    dStatedt_subset = self.dcdt[[self.Cai_ode_idx, self.CaTnC_ode_idx ]] 
    statesTraj_subset = statesi[:,[self.Cai_ode_idx, self.CaTnC_ode_idx ]]
    tsteps_ms = tsteps*self.time_ODE_to_PDE



    return (statesF_subset,dStatedt_subset,tsteps_ms,statesTraj_subset)

  # store values into steady state 
  def updateStates(self,cCai,cCaTnC):
    old = np.copy(self.statesPrev)
    self.statesPrev[self.Cai_ode_idx] = cCai
    self.statesPrev[self.CaTnC_ode_idx] = cCaTnC

    print "ODE/PDE diff: [%f,%f]" % (self.statesPrev[self.Cai_ode_idx]-old[self.Cai_ode_idx],
                                     self.statesPrev[self.CaTnC_ode_idx]-old[self.CaTnC_ode_idx])

  def getStates(self):
    return (self.statesPrev[self.Cai_ode_idx],self.statesPrev[self.CaTnC_ode_idx] )

# <codecell>

##
## PLOTTING
##


## TEMPLATE adjust to your liking 
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
    print "T %f [ms] Conc(%d) %f [mM]" % (t,i,results.uAvg[i])
    print "T %f [ms] Conc(%d) %f [uM]" % (t,i,results.uAvg[i]*mM_to_uM)

    assert(results.uAvg[i]>0), "Conc is negative!!"

  results.c.append(results.uAvg)

  c1 = results.uAvg[Cai_pde_idx]#get last appended value 
  cb1 = results.uAvg[CaTnC_pde_idx]

  return (c1,cb1)

def Report(problem,results,anisotropic=False):

  u = results.u
  t = results.t
  file= results.file
  #file << (u.split()[0], t)
  file << (u,t)
  

  # extract values at certain regions
  ## copy from atpSarcomereVanBeek if needed 
  

# <codecell>


##
## PDE stuff
##

sslMarker=3

## contains parameters for simulation 
## TEMPLATE 
class cparams:
  plot=True
  pvd = False
 
  # copy over ODE models parameters
  odeParams = model.default_parameters() 
  
  meshdebug = True  
  tag = ""
  # 
  ## TODO rename 
  Caiasedistromode="uniform"
  ckdistromode="uniform"
  #Caiasedistromode="nonuniform"

  # Chose value from smallest spacing (maximal contraction) 
  #16.43 &   [0.57  0.58  0.82] & 2.4 & 0.45 \\ 
  #  (Dz,Dx,Dy)
  Disotropic = np.ones(3)
  Danisotropic = np.array([0.57, 0.58, 0.82])

  # Cai 
  DCai=0.290 # [um^2/ms]

  # CaTnC 
  DCaTnC=0.000 # [um^2/ms]

  # for half-sacomere mesh 
  xMax= 2.   # longitudinal (along sarco)
  yMax = 2.0 # circumferential 
  #zMax = 100.  # transverse (along TT) 
  zMax = 2.0    # estiamte assuming cyldrical bundles are 1.2 um in radius
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

## This represents location of dyadic cleft 
class SSL(SubDomain): 
  def inside(self,x,on_boundary):
    result = x[1]>(self.params.yMax-DOLFIN_EPS) and on_boundary
    #print x,result 
    return result

# <codecell>

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

  ##
  ## Time params 
  ## 
  # ODE steps for equilibration
  t0ode = 0
  tFode= 3e3 # [ms] # for 300 ms stim 
  if(debugLevel>0):
    tFode= 1e3 # [ms]
  
  # pde steps for run
  t0pde=tFode
  tFpde=t0pde + duration
  dt = 10. # [ms]
  if(debugLevel==100):
    dt = 1.  # [ms]

  ##
  ## ODE
  ##
  ## Initialize model/get steadystate 
  odeModel = ODEModel(t0=t0ode,tF=tFode)
  odeModel = ODEModel(t0=t0ode,tF=tFode,params=params.odeParams) 
  odeModel.setup() 
  #print odeModel.getStates()
  # here we propagate model just so that we can compare exact and PDE solutions 
  # but we'll need to 'reset' the statesPrev variable
  statesPrev = np.copy(odeModel.statesPrev)
  (dummy,dummy,tstepsExact,statesTrajExact)= odeModel.propagateStates(\
    t0pde,tFpde,incr=1e4)
  odeModel.statesPrev = statesPrev
  # get initial concentrations from ODE model 
  c0i,cb0i = odeModel.getStates()

  if(debugLevel==100):
    c0i = 1; cb0i = 1; # [mM] 
  
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
  indcCai, indcCaTnC= set(), set()
  dm_cCai, dm_cCaTnC = ME.sub(0).dofmap(), ME.sub(1).dofmap()
  for cell in cells(mesh):
      indcCai.update(dm_cCai.cell_dofs(cell.index()))
      indcCaTnC.update(dm_cCaTnC.cell_dofs(cell.index()))
  indcCai = np.fromiter(indcCai, dtype=np.uintc)
  indcCaTnC = np.fromiter(indcCaTnC, dtype=np.uintc)

  cCai_values = u.vector()[indcCai]
  cCaTnC_values = u.vector()[indcCaTnC]

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
  jboundaryExpr = empty()
  jvolExpr = empty()
  jboundaryExpr.Cai = Expression("j",j=0)
  jboundaryExpr.CaTnC = Expression("j",j=0)
  jvolExpr.Cai = Expression("j",j=0)
  jvolExpr.CaTnC = Expression("j",j=0)

  ## Decide on diff const 
  Dii = params.Disotropic

  # homog of myfilaments has 'z' as the long axis, but in half-sarcomere, 'x' is the long axis 
  Dij = Constant([Dii[2],Dii[0],Dii[1]])
  Dij = diag(Dij)
  Dij_Cai= Constant(params.DCai) * Dij
  Dij_CaTnC = Constant(params.DCaTnC)*Dij # 0.0


  ## PDE weak form 
  # diffusion term  
  RHS1 = -inner(Dij_Cai*grad(c),grad(q))*dx
  # TNC has Dij=0. 
  RHS2 = -inner(Dij_CaTnC*grad(cb),grad(v))*dx

  # buffers (TODO double check, since I don't think this is correct) 
  # based on coupledReactionDiffusion.py in my scripts mercurial repo 
  # Reaction: b + c --kp--> cb,  
  #           b + c <-km--- cb
  # TODO add real parameters 
  reaction = False
  if reaction: 
    kp = odeModel.params[odeModel.kon_tnc_idx]
    km = odeModel.params[odeModel.koff_tnc_idx]
    bT = odeModel.params[odeModel.bmax_tnc_idx]
    R = np.array([
      [-kp,km],     # s
      [kp,-km]      # p
    ])
    RHS1 += (R[0,0]*(bT-cb)*c*q + R[0,1]*cb*q)*dx
    RHS2 += (R[1,0]*(bT-cb)*c*v + R[1,1]*cb*v)*dx

  # volumetric flux
  # This is a work-around (see below)
  # Caiase 
  ## TODO rename as SERCAdistro
  if(params.Caiasedistromode=="uniform"):
    RHS1 +=  jvolExpr.Cai*q*dx
    # no vol flux on tnC
    RHS2 += -jvolExpr.CaTnC*v*dx
  else:
    raise RuntimeError("ouch") 
  
    
  ## boundary
  subdomains = MeshFunction("uint",mesh,2)
  boundary = SSL()
  boundary.params = params 
  boundary.mark(subdomains,sslMarker)
  ds = Measure("ds")[subdomains]


  # vol/sa ratio to rescale surf. fluxes
  vol = assemble(Constant(1.)*dx,mesh=mesh)
  sa  = assemble(Constant(1.)*ds(sslMarker),mesh=mesh)
  vol_sa_ratio = vol/sa
  print "vol_sa_ratio  [um]",  vol_sa_ratio

  # boundary flux terms 
  print "WARNING: not working w eval func"
  jL1 =  jboundaryExpr.Cai
  jL2 =  jboundaryExpr.CaTnC
  RHS1  += jL1*q*ds(sslMarker)
  RHS2  += jL2*v*ds(sslMarker)


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
  results.t=t  
#  problem = empty()

  # Initializations for results 
  # See CaiSarc...py
  #results.x_midsarc_aband = np.array([0.50,0.50,0.5])
  results.ts = []
  results.c = []
  problem.mesh=mesh
  problem.u = u
  problem.params = params

  results.ME = ME
  results.u = u
  results.file = file


  # initialize 
  c0F = c0i

  ## Time stepping
  stateFluxes = [] # TODO remove
  jboundarys=[]
  jvols=[]
  vals=[[],[]]
  jSRs = []; 
  cs=[]
  cbs=[]
  cFs=[]
  ts = np.linspace(t0pde+dt,tFpde,(tFpde-t0pde)/dt)
  #print "ts", ts
  
  (c1,cb1)=GlobalConc(problem,results)
  #print t0, tFp
  #print ts
  for i,tf in enumerate(ts):
    print "###########################################################################"
    ## update ODE 
    # predict 'forward' solution from ode model '
    #tf = t0+dt
    t0 = tf - dt
    print "t0 ", t0, "tf ", tf

    (statesF,stateFlux,dummy,dummy) = odeModel.propagateStates(t0,tf,dt)
    print "ODE statesF: [mM] ", statesF
    print "ODE statesF: [uM] ", statesF[Cai_pde_idx] * mM_to_uM 
    print "ODE states Flux [mM/ms]: ", stateFlux
    print "ODE states Flux [uM/ms]: ", stateFlux * mM_to_uM
    ## TODO remove me? 
    stateFluxes.append(stateFlux)

  
    ## update PDE 
    # define flux 
    if(mode=="totalFlux"): 
      CaiFlux = stateFlux[Cai_pde_idx] 
      jboundaryExpr.Cai.j =Jboundary(CaiFlux,vol_sa_ratio=vol_sa_ratio)
      #jboundaryExpr.j *= dt # I think this is already reflected in PDE weak form 

    # separate fluxes is useful when volumetric fluxes are heterogeneously distributed
    elif(mode=="separateFlux"): 
      # apply fluxes 
      # TODO need to generalize 
      jboundaryExpr.Cai.j =Jboundary(odeModel.jBoundary[odeModel.Cai_ode_idx],vol_sa_ratio=vol_sa_ratio)
      jvolExpr.Cai.j =odeModel.jVol[odeModel.Cai_ode_idx]
      #print jboundaryExpr.Cai.j / vol_sa_ratio 
      #print jvolExpr.Cai.j 
      if(debugLevel==100):
        jboundaryExpr.Cai.j =Jboundary(0.1,vol_sa_ratio=vol_sa_ratio)
        jvolExpr.Cai.j =0.
        #jboundaryExpr.j = 0.0; jvolExpr.j = 0.1
      #jboundaryExpr.j *= dt # I think this is already reflected in PDE weak form 
      #jvolExpr.j *= dt # I think this is already reflected in PDE weak form 

      #print "sdfdsf %f [uM/ms] " % ((jboundaryExpr.Cai.j / vol_sa_ratio + jvolExpr.Cai.j)*mM_to_uM)
      #print "sdfdsf %e [M/ms] " % ((jboundaryExpr.Cai.j / vol_sa_ratio + jvolExpr.Cai.j)) 
      print "WARNING: jSR and dCa are lumped tother into jVold"
      jSRs.append(odeModel.jSR)


      # tests to make sure implemented correctly 
      #print "WARNING: hack " 
      #print stateFlux
      #print odeModel.dcdt_Cai
      #print odeModel.jBoundary[Cai_idx]
      #print odeModel.jVolCaiase[Cai_idx]
      ## try all as boundary (works) 
      #CaiFlux = odeModel.jBoundary[Cai_idx] + odeModel.jVolCaiase[Cai_idx]
      #jboundaryExpr.j =Jboundary(CaiFlux,vol_sa_ratio=vol_sa_ratio)
      #jvolCaiaseExpr.j = 0.
      ## try all as volume 
      #jboundaryExpr.j =0.
      #jvolCaiaseExpr.j = CaiFlux
      
      

    # get integrated fluxes 
    totjs = assemble(jboundaryExpr.Cai*ds(sslMarker),mesh=mesh)
    jboundary = totjs/sa * 1/vol_sa_ratio
    print "WARNING: I probably need to constrain jvolExpr s.t. the assembled value is equal to the ODE value"
    totjs = assemble(jvolExpr.Cai*dx,mesh=mesh)
    jvol = totjs/vol
    print "jboundary [mM/ms]", jboundary
    print "jboundary [uM/ms]", jboundary*mM_to_uM
    print "jvol [mM/ms] ", jvol
    print "jvol [uM/ms] ", jvol *mM_to_uM*mM_to_uM

    print "##\n## WARNING: we are not updating CaTnC!!!\n##"


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

    #print "rename c1,cb1"
    #c1 = pdeVals[Cai_idx]#get last appended value 
    #cb1 = pdeVals[CaTnC_idx]
    #ts.append(t)
    jboundarys.append(jboundary)
    jvols.append(jvol)
    results.jdiffs= jboundarys
    cFs.append(statesF[Cai_pde_idx]) # just want 'Cai' from odeModel
    #cs.append(c1)
    #cbs.append(cb1)

    ## Package results 
    Report(problem,results,anisotropic=anisotropic)


    ## store values in ODE  
    #ts.append(tf)
    t0 = tf
    #if(debugLevel> 0):
    if(1): 
      print "WARNING: skipping conc updatee"
      1
    else:
      # TODO I CAN ONLY DO THIS ONCE TnC is appropriately updated
      odeModel.updateStates(c1,cb1)


  ## END LOOP 


#  print "WARNING: merge into plotting"
  # convert lists into arrays
  results.ts = ts
  results.cFs = np.asarray(cFs)
  results.cs = np.asarray(results.c)[:,Cai_pde_idx]
  results.cbs = np.asarray(results.c)[:,CaTnC_pde_idx]
  results.jSRs = np.asarray(jSRs)
  
#
#  ## Plot 
  if(0): 
    fig=plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("[Cai]/[CaTnC] vs time") 
    ax1.plot(tstepsExact,statesTrajExact[:,Cai_pde_idx],label="[Cai] Exact")
    ax1.plot(ts,results.cs,'k.',label="[Cai] (PDE) ")
    ax1.plot(ts,results.cFs,label="[Cai]F (forward ODE solutions)")
    ax1.set_ylabel("[Cai] [uM]") 
    ax1.set_xlabel("t [ms]") 
    ax1.legend(loc=0)
    ax2 = ax1.twinx()
    ax2.set_ylabel("[CaTnC] [uM]") 
    ax2.plot(ts,results.cbs,'r-.',label="cb (PDE) ")
    plt.gcf().savefig(mode+"conc.png")
#
#
#  plt.figure()
#  jboundarys = np.asarray(jboundarys)
#  jvolCaiases = np.asarray(jvolCaiases)
#  #plt.plot(ts,jsyns,label="jhyds")
#  plt.plot(ts,jboundarys,label="jboundarys")
#  plt.plot(ts,jvolCaiases,label="jvolCaiases")
#  plt.legend(loc=0)
#  plt.gcf().savefig(mode+"js.png")

  ext=statesTrajExact[-1,Cai_pde_idx]
  est = c1  # c1 is updated each frame, so this is the most recent
  m = "mode %s: Exact %f != Est %f "%(mode,ext,est)
  if(asserts):
    ext=statesTrajExact[-1,Cai_pde_idx]
    #est=cs[-1]
    est = c1  # c1 is updated each frame, so this is the most recent
    m = "mode %s: Exact %f != Est %f "%(mode,ext,est)
    assert(np.abs(ext-est)<1.0), m
    print "### Mode %s passes!" % mode
  elif (np.abs(ext-est)>1.0):
    print "WARNING: ", m 

  return problem,results 
  

    

# <codecell>

###
### Testing plotting etc 
###
  
# plots isotropic/anisotropic results
def plotting(params,resultsiso,resultsaniso):
  1
  # see CaiSarc...py

# Simple run to plot what Concntrations s.b. 
def TestODE():
  t0 = 0.
  tFp = 1e3 # [ms] 
  incr = 1e3
  odeModel = ODEModel()                
  odeModel.setup()
  (dummy,dummy,tstepsExact,statesTrajExact)= odeModel.propagateStates(t0,tFp,incr=incr) 
  plt.plot(tstepsExact,statesTrajExact[:,odeModel.Cai_pde_idx],label="Cai")              
  plt.legend(loc=0)
  #plt.xlim([9.2e3,10.0e3]) 


  print "WARNING: not a unit test"
  plt.gcf().savefig("testode.png") 


# This needs to be tightly integrated with 'separate.py' to make sure bookkeeping is 
# correct 
def TestODE2():
  t0 = 0.
  tFp = 1e3 # [ms] 
  dt=25 

  odeModel = ODEModel()
  odeModel.setup()
  
  ts = np.linspace(t0+dt,tFp,(tFp-t0)/dt)
  
  jConcs=[]
  jSums=[]
  for i,tf in enumerate(ts): 
    t0 = tf - dt 
    (statesF,stateFlux,dummy,dummy) = odeModel.propagateStates(t0,tf,dt)
    jSum = odeModel.jBoundary + odeModel.jVol
    #print odeModel.jBoundary[odeModel.Cai_ode_idx]
    #print odeModel.jVol[odeModel.Cai_ode_idx]
    jConcs.append(stateFlux[odeModel.Cai_pde_idx])
    jSums.append(jSum[odeModel.Cai_ode_idx])
  jConcs = np.asarray(jConcs)
  jSums  = np.asarray(jSums)


  plt.plot(jConcs,'b-',label="Flux based on change in concentration") 
  plt.plot(jSums,'b.',label="Sum of individual fluxes") 
  plt.legend(loc=0)
  plt.gcf().savefig("testODE2.png") 

  assert( meanSquaredError(jConcs,jSums) < 1e10), "jConc/jSums don't agree"

def meanSquaredError(a,b):
  return np.mean((a-b)**2)

def loop(p,duration=1e3,asserts=True,case="default",mode="separateFlux"):
  if(debugLevel>0):
    duration=1e3
    
  p.tag = case
  problemi,resultsi= runPDE(params=p,mode=mode,anisotropic=False,\
                         asserts=asserts,duration=duration)


  # pickle 
  import cPickle as pickle
  resultsi.pickleName = case+".pkl"
  data1 = {'results':resultsi}
  data1 = {'ts':resultsi.ts,'cs':resultsi.cs,'cbs':resultsi.cbs,'jSRs':resultsi.jSRs}
  output = open(resultsi.pickleName,    'wb')
  pickle.dump(data1, output)
  output.close()
 
  
  return problemi, resultsi


# Basic test + consistency 
def Test1(do="all"):
    cparams.plot = False
    pi = cparams()

    # Kmf = 0.246e-3 [mM]
    # stim_period = 300 [ms]  
    #print pi.odeParams[param_indices("Kmf")]

    #problemi, resultsi = loop(pi,case="test1",mode="totalFlux")  
    problemi, resultsi = loop(pi,case="test1",mode="separateFlux") 

    # check on first/last [Ca] entry 
    i=0
    #print resultsi.cs[i] - 0.000299
    assert(np.abs(resultsi.cs[i] - 0.000299) < 1e-6), "FAIL: " % resultsi.cs[i]
    i=-1
    assert(np.abs(resultsi.cs[i] - 0.000182) < 1e-6), "FAIL: " % resultsi.cs[i]
    #print resultsi.cs[-1]
    print '## PASSED Test1' 
    
    return (problemi, resultsi)

  

# <codecell>

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
      #debugLevel=100
      debugLevel=10
      Test1()
    if(arg=="-test1"):
      1
    if(arg=="-validation"):
      #TestODE()
      #TestODE2()
      Test1()
      

