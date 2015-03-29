#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
#import matplotlib.pylab as plt
class empty:pass 

idxCa = 0 # PDE 
idxBuff = 1
idxFluo = 2
nComp = 1 # compartments
nSpec = 3
nDOF=nComp*nSpec

## Units
nm_to_um = 1.e-3


class LeftTT(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmin[0]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class RightTT(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmax[0]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[2]- self.mmax[2]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

class Params():
  paraview = False 
  verbose=False 

  # time steps 
  steps = 500
  dt = 1.0   # [ms] 



  # diffusion params 
  #DAb   = Dbulk# [um^2/ms] Diff const within PDE (domain 1) 

  # init concs 
  cCainit =1.0  # [uM]
  cBuffinit =1.0  # [uM]
  cFluoinit =1.0  # [uM]


class InitialConditions(Expression):

  def eval(self, values, x):
    # edge  
    #if (x[0] < 0.5 and np.linalg.norm(x[1:2]) < 0.5):
    # corner 
    #oif (np.linalg.norm(x -np.zeros(3) ) < 0.5):
    if 1:
      values[0] = self.params.cCainit
      values[1] = self.params.cBuffinit
      values[2] = self.params.cFluoinit
  def value_shape(self):
    return (nDOF,)


def PrintSlice(mesh,u):  
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
    (gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                       mmin[0]:mmax[1]:(res*1j)]
    from scipy.interpolate import griddata
    img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy))
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

# TODO verify units 
def tsolve(fileName="sarcomere.xml",\
           outName="output.pvd",
	   debug=True):

  params = Params()
  if debug:
    params.T = 10

  # Create mesh 
  mesh = Mesh(fileName)   # not working with xml files I"ve been generatinf 
  mesh = UnitCubeMesh(16,16,16)  
  dim =mesh.ufl_cell().geometric_dimension()

  # function spaces
  V = FunctionSpace(mesh, "Lagrange", 1)
  ME = MixedFunctionSpace([V,V,V]) # Ca, Buffered Ca, Fluo Ca

  # Define trial and test functions
  du    = TrialFunction(ME)
  vCa,vBuff,vFluo  = TestFunctions(ME)

  # Define functions
  u_n   = Function(ME)  # current solution
  u_0  = Function(ME)  # solution from previous converged step

  # split mixed functions
  cCa_n,cBuff_n, cFluo_n = split(u_n)
  cCa_0,cBuff_0, cFluo_0 = split(u_0)

  ## mark boundaries 
  subdomains = MeshFunction("size_t",mesh,dim-1)
  boundary = LeftTT()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  boundary = RightTT()
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
  #bcs.append(bc)
   
  # define integrator 
  ds = Measure("ds")[subdomains]

  
  ## weak form 
  # params 
  D = Constant(Diff) 

  # three D diff const 
  #Dii  = Constant((problem.d_eff[0],problem.d_eff[1],problem.d_eff[2]))
  #Dij = diag(Dii)  # for now, but could be anisotropic
  #a =inner(Dij * grad(u), grad(v))*dx

  ## TODO add realistic rate constants, put into params()
  r0 = Constant(0.1)
  r1 = Constant(1.) 
  r2 = Constant(0.1)
  r3 = Constant(1.)   
  k3 = Constant(0.5)


  # diffusion 
  DCa = params.DCa
  RHS = -inner(DCa*grad(u), grad(q))*dx
  # reaction 
  #RHS+= -r0*u*q*dx
  # RyR release at R TT  
  #RHS+=  r1*q*ds(rMarker,domain=mesh)
  # RyR release at L TT  
  # replace w  j = (c/k)^n / (1 + (c/k)^n)
  #RHS+=  r3*q*ds(lMarker,domain=mesh) 
  # NCX extrusion at SL 
  # make function of u through Expression 
  #if dim==3:
  #  RHS+= -r2*u*q*ds(slMarker,domain=mesh) 
  
  dt = params.dt

  # might ultimately need to do operator splitting if we communicate with an ODE   model. 
  # however, perhaps it can all be communicated through an expression 
  ## TODO make into expression 
  #JRyr = delta0 * Iryr / (2*F*delta*V)
  ## TODO make into Form 
  #Jbuff = alphabuff * u * (Btot - BC) - betabuff*BC   # probably describe this as a  state
  #Jfluo = alphafluo * u * (Ftot - FC) - betafluo*FC   # probably describe this as a   state
  ## TODO Add NCX/SERCA 
  ## 1. through Shannon Model (just importing functions that I need) 
  ## 2. modified with my Markov model 
 

  #form += Dc*inner(grad(u),grad(v))*dx
  #form += JRyR*v*dx   # along boundary? (refer to enzyme kinetics ex)
  #form += -1*Jbuff*v*dx
  #form += -1*Jfluo*v*dx
  #form += DF*inner(grad(uF),grad(vF))*dx
  #form += Jfluo*vF*dx



  L = u_n*q*dx - u_0*q*dx - dt * RHS

  # Compute directional derivative about u in the direction of du (Jacobian)
  a = derivative(L, u_n, du)
  
  ## define solver,params 
  problem = MyEquation(a,L,bcs)
  solver = NewtonSolver()                
  solver.parameters["linear_solver"] = "gmres"
  solver.parameters["convergence_criterion"] = "incremental"
  solver.parameters["relative_tolerance"] = 1e-6
  file = File(outName, "compressed")
  
  ## Variables for storing answers (after projection) 
  # need to declare outside of loop 
  up = Function(V)
  u0p = Function(V)
  
  t = 0.
  concs=[]
  ts = []
  us = []
  while (t < params.T):
      # advance 
      t0=t
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem,u.vector())

      # remap to u from u' = e^+ * u (see above) 
      file << (u,t) 

      
      # store
      #us.append(PrintSlice(mesh,up))

      # report on prev iter
      #uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      #area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      #conc = uds/area
      #ts.append(t0)
      #concs.append(conc)
      plot(u,rescale=True)
      #interactive()	

  ts = np.asarray(ts)
  concs = np.asarray(concs)

  results = empty()
  results.us = us 
  results.ts = ts
  results.concs = concs 
  results.mesh = mesh 

  
  return (results)

  
  

def valid(fileName="sarcomere.xml"):
  (results) = tsolve(Diff=1.0,fileName=fileName)

import sys

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -valid" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):   
    if(arg=="-valid"):
      valid()
    if(arg=="-valid2D"):
      valid(fileName="sarcomere2D.xml")



