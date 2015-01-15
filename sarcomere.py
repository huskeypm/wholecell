#
# Time dependent solver, single species 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
class empty:pass 


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

class Sarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[2]- self.mmax[2]) < DOLFIN_EPS) 
    #print x[0], edge, on_boundary
    return on_boundary and edge

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


def tsolve(Diff=1.,fileName="sarcomere.xml",\
           outName="output.pvd",
           T=100):
  # Create mesh and define function space
  mesh = Mesh(fileName)     
  V = FunctionSpace(mesh, "Lagrange", 1)
  
  # Define trial and test functions
  du    = TrialFunction(V)
  q  = TestFunction(V)
  # Define functions
  u   = Function(V)  # current solution
  u0  = Function(V)  # solution from previous converged step



  ## mark boundaries 
  dim = 3 
  subdomains = MeshFunction("size_t",mesh,dim-1)
  boundary = LeftTT()
  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  boundary = RightTT()
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  rMarker = 3
  boundary.mark(subdomains,rMarker)
  boundary = Sarcolemma()
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  slMarker = 4
  boundary.mark(subdomains,slMarker)

  ## decide on source of probability density  
  Crest = 0.1 # Ca [uM] 
  ic = Expression("c",c=Crest)
  u.interpolate(ic)
  u0.interpolate(ic)

  bcs=[]
  #bcs.append(bc)
   
  # define integrator 
  ds = Measure("ds")[subdomains]

  
  ## weak form 
  # params 
  D = Constant(Diff) 
  r0 = Constant(0.1)
  r1 = Constant(1.) 
  r2 = Constant(0.1)
  r3 = Constant(1.)   
  k3 = Constant(0.5)
  t = 0.0
  dt=1.0   
  # diffusion 
  RHS = -inner(D*grad(u), grad(q))*dx
  # reaction 
  RHS+= -r0*u*q*dx
  # RyR release at R TT  
  RHS+=  r1*q*ds(rMarker,domain=mesh)
  # RyR release at L TT  
  # replace w  j = (c/k)^n / (1 + (c/k)^n)
  RHS+=  r3*q*ds(lMarker,domain=mesh) 
  # NCX extrusion at SL 
  # make function of u through Expression 
  RHS+= -r2*u*q*ds(slMarker,domain=mesh) 
  

  L = u*q*dx - u0*q*dx - dt * RHS

  # Compute directional derivative about u in the direction of du (Jacobian)
  a = derivative(L, u, du)
  
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
  
  concs=[]
  ts = []
  us = []
  while (t < T):
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


  ts = np.asarray(ts)
  concs = np.asarray(concs)

  results = empty()
  results.us = us 
  results.ts = ts
  results.concs = concs 
  results.mesh = mesh 

  
  return (results)

  
  

def valid():
  (results) = tsolve(Diff=1.0)

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



