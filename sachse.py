#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
class empty:pass 

idxCa = 0 # PDE 
idxBuff = 1
idxFluo = 2
nComp = 1 # compartments
nSpec = 3
nDOF= nComp*nSpec
eps = 0.05

## Units
nm_to_um = 1.e-3
ttRad = 0.25 # [um]


class LeftTT(SubDomain):
  def inside(self,x,on_boundary):
    #edge = (np.abs(x[0]- self.mmin[0]) < DOLFIN_EPS) 
    # Define TT loc 
    yMid = 0.5*(self.mmax[1]+self.mmin[1])
    centroid = np.array([self.mmin[0],yMid])

    # check if point is nearby 
    d = np.linalg.norm(centroid-x[0:2])
    isTT = (d < (ttRad+eps))

    #print x,centroid,d,isTT
    #print x[0], edge, on_boundary
    return on_boundary and isTT

class RightTT(SubDomain):
  def inside(self,x,on_boundary):
    #edge = (np.abs(x[0]- self.mmax[0]) < DOLFIN_EPS) 
    yMid = 0.5*(self.mmax[1]+self.mmin[1])
    centroid = np.array([self.mmax[0],yMid])

    # check if point is nearby 
    d = np.linalg.norm(centroid-x[0:2])
    isTT = (d < (ttRad+eps))

    #print x,centroid,d,isTT
    #print x[0], edge, on_boundary
    return on_boundary and isTT

class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[2]- self.mmax[2]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class Params():
  paraview = True  
  verbose=False 

  # time steps 
  T = 500
  T = 1  

  dt = 1.0   # [ms] 



  # diffusion params 
  #DAb   = Dbulk# [um^2/ms] Diff const within PDE (domain 1) 
  DCa = 1.0  # [um^2/ms]
  DBuff = Constant(0.) # unverified
  DFluo = 1.0 
  Ds = [DCa,DBuff,DFluo]

  # init concs 
  cCainit =1.0  # [uM]
  cBuffinit =1.0  # [uM]
  cFluoinit =1.0  # [uM]
  cInits = [cCainit,cBuffinit,cFluoinit]

  # buffer (mostly TnC) 
  alphabuff=1. # kon buff  [1/uMms] 
  Btot = 1. # [Buff] [uM] 
  betabuff = 1. # koff buff [1/ms]
 
  # fluo
  alphafluo=1. # kon fluo  [1/uMms] 
  Ftot = 1. # [Fluo] [uM] 
  betafluo = 1. # koff fluo [1/ms]

class InitialConditions(Expression):

  def eval(self, values, x):
    # edge  
    #if (x[0] < 0.5 and np.linalg.norm(x[1:2]) < 0.5):
    # corner 
    #oif (np.linalg.norm(x -np.zeros(3) ) < 0.5):
    if 1:
      for i,j in enumerate(self.params.cInits):
            #print i 
            values[i] = j
#      values[0] = self.params.cCainit
#      values[1] = self.params.cBuffinit
#      values[2] = self.params.cFluoinit
  def value_shape(self):
    return (nDOF,)

 # TODO put in real info 
def UpdateBeardVals(nC,bM,mgatpi_uM,fadpi_uM ,fmgi_uM ):
  

  ### [mM] ###
  print "mgatpi %f [uM] adpi %f [uM] mgi %f [uM] " % \
   (mgatpi_uM,fadpi_uM,fmgi_uM)
  cmgatpi = mgatpi_uM * uM_to_mM
  cfadpi = fadpi_uM * uM_to_mM
  cfmgi = fmgi_uM * uM_to_mM
  
  # get nucleotide concentreations at vertex i 
  # convert MgATP into free ATP (used by Beard)
  #catpi = nC.mATP_to_fATP( cmgatpi )                
  ctatpi = nC.mATP_to_tATP( cmgatpi,cfmgi)           
  cfatpi = ctatpi - cmgatpi 
  ctadpi = nC.fADP_to_tADP( cfadpi,cfmgi)           
  cmgadpi = ctadpi - cfadpi 

  ## Compute nucleotide flux using beard model" 
  # Current doesn't use [Mg]
  jtatpi,jtadpi,dummy = bM.jNucleotide(ctatpi*mM_to_M,ctadpi*mM_to_M) # ,cmgi)
  # multiply by '-1', since negative flux represents ATP diffusion 
  # into cytosol 
  jtatpi*=-1
  jtadpi*=-1
  #print "%f %f %f" % (ctatpi,ctadpi,jtatpi)

  #jtatpi = jtadpi = 0. ; print "KILLING JATP flux"
  # jNucleotide returns mM/s, 1 mM/s = 1000 uM/1000 ms = 1 [uM/ms]
  print "[fatp] %f [mgatp] %f [tatp] %f [mM] jtATP %f [mM/s]\n"\
    % ( cfatpi,cmgatpi,ctatpi,jtatpi)
  print "[fadp] %f [mgadp] %f [tadp] %f [mM] jtATP %f [mM/s]\n"\
    % ( cfadpi,cmgadpi,ctadpi,jtadpi)


  # j = dN/dt --> dN = j*dt
  # j [mM/s] <==> j[uM/ms]
  # dt [ms] * j [uM/ms] = delatp [uM] --> delatp [mM]
  deltatpi = params.dt * jtatpi * uM_to_mM
  deltadpi = params.dt * jtadpi * uM_to_mM


  # first will assume that total nucleotide is added (substracted) 
  # from the free nucleotide; 
  # When we add del_FreeATP to the FreeATP concentration, we push 
  # the Mg+fATP<-->MgATP from equilibrium. Therefore, we recompute
  # the system at equlibrium  
  cfatpn = cfatpi + deltatpi
  cfadpn = cfadpi + deltadpi
  #print "fatpn %f (before equil) " % cfatpn


def PrintSlice(mesh,u,dims=3):    
#    mesh = results.mesh
    #dims = np.max(mesh.coordinates(),axis=0) - np.min(mesh.coordinates(),axis=0)
    mmin = np.min(mesh.coordinates(),axis=0)
    mmax = np.max(mesh.coordinates(),axis=0)
    
    #u = results.u_n.split()[0]
    #u = results.u_n
    up = project(u,FunctionSpace(mesh,"CG",1))
    res = 100
    #(gx,gy,gz) = np.mgrid[0:dims[0]:(res*1j),
    #                      dims[1]/2.:dims[1]/2.:1j,
    #                      0:dims[2]:(res*1j)]
    if dims==3:
      (gx,gy,gz) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                            mmin[1]:mmax[1]:(res*1j),  
                                  0:0:1j]
      img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))
      print np.shape(img0)
      img0 = np.reshape(img0[:,:,0],[res,res])
    else: 
      raise RuntimeError("Not supported") 

    # test
    if 0: 
      plt.pcolormesh(img0)
      plt.colorbar()
      plt.gcf().savefig("test.png")
    return img0


# Nonlinear equation 
class MyEquation(NonlinearProblem):
    def __init__(self, a, L,bcs):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bcs = bcs
    def F(self, b, x):
        assemble(self.L, tensor=b)
        for bc in self.bcs:
          bc.apply(b,x)
    def J(self, A, x):
        assemble(self.a, tensor=A)
        for bc in self.bcs:
          bc.apply(A)



# In[18]:

# TODO verify units 
def tsolve(fileName="sarcomere.xml",\
           outName="output.pvd",\
           mode="sachse",\
           params = Params(),\
	   debug=False): 
 
  if debug:
    params.T = 10

  # Create mesh 
  if debug: 
    mesh = UnitCubeMesh(16,16,16)  
    mesh = UnitCubeMesh(5,5,5)
  else: 
    mesh = Mesh(fileName)   # not working with xml files I"ve been generatinf 
  dim =mesh.ufl_cell().geometric_dimension()

  # function spaces
  V = FunctionSpace(mesh, "Lagrange", 1)
  Vs = []
  for i in range(nDOF):
        #print i 
        Vs.append(V)
  ME = MixedFunctionSpace(Vs) # Ca, Buffered Ca, Fluo Ca

  # Define trial and test functions
  du    = TrialFunction(ME)
  vCa,vBuff,vFluo  = TestFunctions(ME)
  vs = [vCa,vBuff,vFluo]

  # Define functions
  u_n   = Function(ME)  # current solution
  u_0  = Function(ME)  # solution from previous converged step

  # split mixed functions
  #print "WARNING: not sure how to generalize"       
  cCa_n,cBuff_n, cFluo_n = split(u_n)
  cns = [cCa_n,cBuff_n, cFluo_n]
  cCa_0,cBuff_0, cFluo_0 = split(u_0)
  c0s = [cCa_0,cBuff_0, cFluo_0]

  ## mark boundaries 
  # problem specific          
  subdomains = MeshFunction("size_t",mesh,dim-1)
  subdomains.set_all(0)
  boundary = LeftTT()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  boundary = RightTT()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  rMarker = 3
  boundary.mark(subdomains,rMarker)
  if dim>2:
    boundary = OuterSarcolemma()
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = 4
    boundary.mark(subdomains,slMarker)

  ## decide on source of probability density  
  init_cond = InitialConditions()
  init_cond.params = params
  u_n.interpolate(init_cond)
  u_0.interpolate(init_cond)

  bcs=[]
  if 1: 
    bc = DirichletBC(ME.sub(0),Constant(1.),subdomains,lMarker)
    bcs.append(bc)
    # doesn't seem to work with ME.sub
    #import view 
    #bc = DirichletBC(V,Constant(1.),subdomains,lMarker)
    #bcs.append(bc)
    ##view.PrintBoundary(mesh,bcs)
    #marked1 = Function(V)
    #bc.apply(marked1.vector())
    #File("marked.pvd") << marked1
    #quit()
    #view.PrintSubdomains(mesh,subdomains)
   
  # define integrator 
  ds = Measure("ds")[subdomains]

  
  ## weak form 
  # params 

  # three D diff const 
  #Dii  = Constant((problem.d_eff[0],problem.d_eff[1],problem.d_eff[2]))
  #Dij = diag(Dii)  # for now, but could be anisotropic
  #a =inner(Dij * grad(u), grad(v))*dx

  ### RHS 
  ## diffusion 
  Ds  = params.Ds     
  RHSis=[]
  for i in range(nDOF):
    #RHSi= -inner(Ds[i]*grad(cns[i]), grad(vs[i]))*dx
    RHSi= -inner(Ds[i]*grad(u_n[i]), grad(vs[i]))*dx
    RHSis.append(RHSi)

  ## buffreactions
  if mode=="sachse":
    #- buffer 
    kp = params.alphabuff 
    btot = params.Btot 
    km = params.betabuff   # probably describe this as a state
    Jbuff = kp*cCa_n*(btot - cBuff_n) - km*cBuff_n   # probably describe this as a state
    RHSis[idxCa] += -Jbuff*vCa*dx
    RHSis[idxBuff] += +Jbuff*vBuff*dx
  
    #- fluo 
    kp = params.alphafluo 
    ftot = params.Ftot 
    km = params.betafluo   # probably describe this as a state
    Jfluo = kp*cCa_n*(ftot - cFluo_n) - km*cFluo_n   # probably describe this as a state
    RHSis[idxCa] += -Jfluo*vCa*dx
    RHSis[idxFluo] += +Jfluo*vFluo*dx
  
    #Jfluo = alphafluo * u * (Ftot - FC) - betafluo*FC   # probably describe this as a state

    ## boundary reactions 
    # RyR
    # TODO add NCX, etc 
    delta0 = 1.  # TBD
    iryr = Expression("exp(-t)",t=0)
    delta = 1.
    F = 1.
    V = 1.
    jRyR = delta0 * iryr / (2*F*delta*V)
    #jRyR = delta0 / (2*F*delta*V)
    #print "jRyR ", jRyR
    # Add left/right TTs
    RHSis[idxCa] += jRyR*vCa*ds(lMarker,domain=mesh)
    RHSis[idxCa] += jRyR*vCa*ds(rMarker,domain=mesh)

    # NCX
    r = Constant(0.001)
    jSarco= r
    RHSis[idxCa] += jSarco*cns[idxCa]*vs[idxCa]*ds(slMarker,domain=mesh)
  

    # SERCA (here we assume SERCA is one of the cyto fluxes)  
    r = Constant(0.001)
    JCyto = r*cns[idxCa]*vs[idxCa]*dx
    RHSis[idxCa] += JCyto

  # eerything communicated through ode model 
  #elif mode=="satin" or mode=="despa":
  elif mode=="ode":
    # TODO 
    jCa = Expression("r",r=0)
    jCa.r = 0.1
    RHSis[idxCa] += -jCa*vCa*dx
  else:
    raise RuntimeError("Need to give a mode option") 

  
  
  # RyR release at R TT  
  #RHS+=  r1*q*ds(rMarker,domain=mesh)
  # RyR release at L TT  
  # replace w  j = (c/k)^n / (1 + (c/k)^n)
  #RHS+=  r3*q*ds(lMarker,domain=mesh) 
  # NCX extrusion at SL 
  # make function of u through Expression 
  #if dim==3:
  #  RHS+= -r2*u*q*ds(slMarker,domain=mesh) 
  
    
  ### LHS    
  dt = params.dt
  #LHS = (cCa_n-cCa_0)*vCa/dt * dx - RHSCa
  LHSis=[]
  L=0
  for i in range(nDOF):
    # TODO DBL check 
    #LHSi= (cns[i]-c0s[i])*vs[i]/dt * dx - RHSis[i]     
    LHSi= (u_n[i]-u_0[i])*vs[i]/dt * dx - RHSis[i]     
    LHSis.append(RHSi)
    L+= LHSi
  #i = 0 
  #LHSi0 = (cns[i]-c0s[i])*vs[i]/dt * dx - RHSis[i]
  #i = 1 
  #LHSi1 = (cns[i]-c0s[i])*vs[i]/dt * dx - RHSis[i]
  #i = 2 
  #LHSi2 = (cns[i]-c0s[i])*vs[i]/dt * dx - RHSis[i]
  #L = LHSi0 + LHSi1 + LHSi2
    
  # Compute directional derivative about u in the direction of du (Jacobian)
  # (for Newton iterations) 
  a = derivative(L, u_n, du)
  
  
  # Create nonlinear problem and Newton solver
  problem = MyEquation(a, L,bcs)
  solver = NewtonSolver()
  solver.parameters["linear_solver"] = "lu"
  solver.parameters["convergence_criterion"] = "incremental"
  solver.parameters["relative_tolerance"] = 1e-6
  
  # Output file
  if params.paraview:     
    file = File("output.pvd", "compressed")
    
  t = 0.
  concs=[]
  ts = []
  us = []
  if 0: 
      i=0
      plot(u_n.sub(i),rescale=True)
      interactive()	

  while (t < params.T):
      # advance 
      t0=t
      u_0.vector()[:] = u_n.vector()
      solver.solve(problem,u_n.vector())

      # store 
      file << (u_n,t) 

      
      # store
      #us.append(PrintSlice(mesh,up))

      # report on prev iter
      #uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      #area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      #conc = uds/area

      vol = assemble(Constant(1.)*dx(domain=mesh))
      assembles=np.zeros(nDOF)
      for i in range(nDOF):
        #assembles[i] = assemble(cns[i]*dx(domain=mesh))
        assembles[i] = assemble(u_n[i]*dx(domain=mesh))
      concis  = assembles/vol  
      if MPI.rank(mpi_comm_world())==0:
        print "Conc Ca", concis[idxCa]
 
      ts.append(t0)
      concs.append(concis)

      # update
      if mode=="sachse":
        iryr.t = t 
      if mode=="ode":
        # TODO run ODE model 
        # aggregate state variables
        # need to assemble concs over entire vol 
        #def UpdateBeardVals(nC,bM,mgatpi_uM,fadpi_uM ,fmgi_uM ):
        jCa.r  = 0. # TODO update with ode 
      t += dt


  #
  results = empty()
  results.mesh = mesh 
  results.concs = np.asarray(concs)
  #if MPI.rank(mpi_comm_world())==0:
  if 0: # Doesn't work inside ipython notebook 
    results.caSlice = PrintSlice(mesh,u_n[i])
  return results     


# In[19]:
def whosalive():
  print "I'm alive!"

def doit(debug=True,mode="",params=Params()):  
  tsolve(debug=debug,mode=mode,params=params)
#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################


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
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-debug"):
      #arg1=sys.argv[i+1] 
      doit(debug=True)     
    elif(arg=="-test0"):
      doit(debug=False)    
      quit()
    elif(arg=="-test1"):
      doit(debug=False,mode="ode")    
      quit()
    elif(arg=="-test2"):
      params = Params()
      params.T = 20
      doit(debug=False,params=params,mode="sachse") 
      quit()
  





  raise RuntimeError("Arguments not understood")




