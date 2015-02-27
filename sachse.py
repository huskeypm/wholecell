#
# Time dependent solver, single species 
# 
# TODO get concentrations of buffs, etc correct
#
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
class empty:pass 
import sys
sys.path.append("./siam/")

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

### For 4TT geometry
#class TopTT(SubDomain):
#  def inside(self,x,on_boundary):
#    # Define TT loc 
#    centroid1 = np.array([self.mmin[0],self.mmax[1]])
#    centroid2 = np.array([self.mmax[0],self.mmax[1]])
#
#    # check if point is nearby 
#    d1 = np.linalg.norm(centroid1-x[0:2])
#    d2 = np.linalg.norm(centroid2-x[0:2])
#    isTT1 = (d1 < (ttRad+eps))
#    isTT2 = (d2 < (ttRad+eps))
#    isTT = isTT1 or isTT2
#
#    #print x,centroid,d,isTT
#    #print x[0], edge, on_boundary
#    return on_boundary and isTT
#
#class BottomTT(SubDomain):
#  def inside(self,x,on_boundary):
#    # Define TT loc 
#    centroid1 = np.array([self.mmin[0],self.mmin[1]])
#    centroid2 = np.array([self.mmax[0],self.mmin[1]])
#
#    # check if point is nearby 
#    d1 = np.linalg.norm(centroid1-x[0:2])
#    d2 = np.linalg.norm(centroid2-x[0:2])
#    isTT1 = (d1 < (ttRad+eps))
#    isTT2 = (d2 < (ttRad+eps))
#    isTT = isTT1 or isTT2
#
#    #print x,centroid,d,isTT
#    #print x[0], edge, on_boundary
#    return on_boundary and isTT


## For 2TT geometry
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

# [CaB] = [Btot]/(KD/[Ca]+1)
def buffered(totB,freeCa,kFwd,kback):
    KD = kback/kFwd
    return totB/(KD/freeCa + 1) 

class Params():
  useMPI = True
  paraview = True  
  verbose=False 

  # time steps 
  T = 500
  T = 1  

  dt = 1.0   # [ms] 



  # diffusion params 
  #DAb   = Dbulk# [um^2/ms] Diff const within PDE (domain 1) 
  DCa = 0.39  # [um^2/ms] verified
  DBuff = Constant(0.) # TnC
  DFluo = 0.1 # [um^2/ms] verified 
  Ds = [DCa,DBuff,DFluo]


  # buffer (mostly TnC) 
  alphabuff=0.04 # kon verified [1/uMms] 
  Btot = 70. # [TnC] [uM] verified 
  betabuff = 0.04 # koff verfied [1/ms]
 
  # fluo
  alphafluo=0.23 # kon fluo  [1/uMms] verified 
  Ftot = 50. # [Fluo] [uM] verified  
  betafluo = 0.17 # koff fluo [1/ms] verified

  # RyR
  ryrAmp = 9.5 # [pA/pF] 
  ryrOffset = 5 # [ms]
  ryrTau = -50/np.log(1/2.) # half-max amp at 50 ms 
  ryrKm = 0.2 # uM (made this up)  

  # init concs 
  cCainit =0.1  # Verified [uM]
  cBuffinit = buffered(Btot,cCainit,alphabuff,betabuff)  # TnC unverified [uM]
  cFluoinit = buffered(Ftot,cCainit,alphafluo,betafluo)  # [uM]
  cInits = [cCainit,cBuffinit,cFluoinit]

class InitialConditions(Expression):

  def eval(self, values, x):
    # edge  
    #if (x[0] < 0.5 and np.linalg.norm(x[1:2]) < 0.5):
    # corner 
    #oif (np.linalg.norm(x -np.zeros(3) ) < 0.5):
    for i,j in enumerate(self.params.cInits):
            #print i 
            values[i] = j
#      values[0] = self.params.cCainit
#      values[1] = self.params.cBuffinit
#      values[2] = self.params.cFluoinit
  def value_shape(self):
    return (nDOF,)

 # TODO put in real info 


def InterpolateData(mesh,u,dims=3,mode="line",doplot=False):    
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
    #if dims==3:
    
    if mode=="slice":
      (gx,gy,gz) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                            mmin[1]:mmax[1]:(res*1j),  
                                  0:0:1j]
      img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))
      #print np.shape(img0)
      img0 = np.reshape(img0[:,:,0],[res,res])

      if doplot:
        plt.pcolormesh(img0)
        plt.colorbar()
        plt.gcf().savefig("test.png")
      return img0

    if mode=="line":
    
      yMid = 0.5*(mmax[1]+mmin[1])
      (gx,gy,gz) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                            yMid:yMid:1j,  
                            0:0:1j]
      line = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))
      #print np.shape(line)
      #img0 = np.reshape(line[:,0,0],[res])
      line = np.ndarray.flatten(line)
      #print np.shape(line)
      if doplot:
        plt.plot(line)
        plt.gcf().savefig("test1.png") 
      return line 




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
def tsolve(fileName="sarcomere2TT.xml",\
           outName="output.pvd",\
           mode="sachse2TT", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
	   debug=False): 
 
  if debug:
    params.T = 10

  ## Create mesh 
  if debug: 
    mesh = UnitCubeMesh(16,16,16)  
    mesh = UnitCubeMesh(5,5,5)

  if "2D" in mode:
    import sarcomere2DSSL
    if mode=="2D_SSL":
      1 # ME: R,R,Q  
      cm = sarcomere2DSSL.sarcomere2DSSL(mode="wSSL")
          
    if mode=="2D_noSSL":
      1 # ME: R,Q    
        # p/rxn for SSL 
      cm = sarcomere2DSSL.sarcomere2DSSL(mode="woSSL")

  elif "sachse" in mode: 
    if mode=="sachse4TT":
      import sarcomere4TT
      cm = sarcomere4TT.sarcomere4TT()
    else: 
      import sarcomere2TT
      cm = sarcomere2TT.sarcomere2TT()

  elif "satin" in mode:
    import sarcomereSatin
    cm = sarcomereSatin.sarcomereSatin()
    
    

  mesh = cm.GetMesh()
  quit()
  dim =mesh.ufl_cell().geometric_dimension()

  # function spaces
  V = FunctionSpace(mesh, "Lagrange", 1)
  ## PKH R = FunctionSpace(mesh,"R",0) # for each compartment 
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
 
  if mode!="sachse4TT":
    boundary = LeftTT()
  else:
    import sarcomere4TT as cm
    boundary = cm.BottomTT()

  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  lMarker = 2
  boundary.mark(subdomains,lMarker)
  
  if mode!="sachse4TT":
    boundary = RightTT()
  else: 
    boundary = cm.TopTT()

  boundary.mmin = np.min(mesh.coordinates(),axis=0)
  boundary.mmax = np.max(mesh.coordinates(),axis=0)
  rMarker = 3
  boundary.mark(subdomains,rMarker)

  if dim>2:
    boundary = cm.OuterSarcolemma()
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = 4
    boundary.mark(subdomains,slMarker)

  ## decide on source of probability density  
  init_cond = InitialConditions()
  init_cond.params = params
  u_n.interpolate(init_cond)
  u_0.interpolate(init_cond)

  bcs=[]
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
  # Test flux locations 
  # PKH verified 150216
  if mode=="bcs":
    ## Might ultimately put BC for mitochondria (top/bottom) here  

    #bc = DirichletBC(ME.sub(0),Constant(1.),subdomains,lMarker)
    #bcs.append(bc)
 
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
    1 

  if mode=="sachse4TT":
#    # TT radius 
#    xl= 0.25
#    sc = 0.01
#    #lhs = 1/(1+np.exp((xs-xl)/sc))    
#    xr = 2-0.25
#    #rhs = 1-1/(1+np.exp((xs-xr)/sc))    
#
#    # location of zline 
#    zLine = Expression("1/(1+exp((x[0]-xl)/sc)) + 1-1/(1+exp((x[0]-xr)/sc))", \
#      xl=xl,xr=xr,sc=sc)
#    # location of cytosol
#    cytosol = Expression("(1-1/(1+exp((x[0]-xl)/sc)))*(1/(1+exp((x[0]-xr)/sc)))", \
#      xl=xl,xr=xr,sc=sc)

    zLine = cm.zLine
    cytosol= cm.cytosol

    pBuff = zLine
    #plt.plot(xs,)
    #plt.plot(xs,lhs+rhs)
    #plt.plot(xs,1-(lhs+rhs))
  else:
    pBuff = Expression("1.")
   
  # works for either sachse2 or sachse4TT
  if "sachse" in mode:
    #- buffer 
    kp = params.alphabuff 
    btot = params.Btot 
    km = params.betabuff   # probably describe this as a state
    Jbuff = pBuff*(kp*cCa_n*(btot - cBuff_n) - km*cBuff_n)  # probably describe this as a state
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
    iryr = Expression("a*exp(-(t-to)/tau)",a=params.ryrAmp,\
                                    to=params.ryrOffset,\
                                    tau=params.ryrTau,\
                                    t=0)
    iryrDefect = Expression("a/(1+pow((Km/ca),n))",a=params.ryrAmp,
                                               Km=params.ryrKm,
                                               n=5,
                                               ca = 0.1)
    # Torres description 
    #delta0 = 1.  # TBD
    #delta = 1.
    #F = 1.
    #V = 1.
    #jRyR = delta0 * iryr / (2*F*delta*V)
   
    # converted into 3D using fluxes.pdf
    # TODO: warning: i think this needs to be updated according to the 
    # SA overwhich jRyr is applied. Check eqn 7 in fluxes.pdf
    jRyR = 0.1*iryr 
    jRyRDefect = 0.1*iryrDefect


    #jRyR = delta0 / (2*F*delta*V)
    #print "jRyR ", jRyR
    # Add left/right TTs
    ttConfig = "both"
    #RHSis[idxCa] += jRyR*vCa*ds(lMarker,domain=mesh)
    #RHSis[idxCa] += jRyR*vCa*ds(rMarker,domain=mesh)
    ttConfig = "normalL,defectR"
    RHSis[idxCa] += jRyR*vCa*ds(lMarker,domain=mesh)
    RHSis[idxCa] += jRyRDefect*vCa*ds(rMarker,domain=mesh)

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

  file << (u_n,t) 
  while (t < params.T):
      # advance 
      t0=t
      u_0.vector()[:] = u_n.vector()
      solver.solve(problem,u_n.vector())

      # store 
      file << (u_n,t) 

      
      # store

      # report on prev iter
      #uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      #area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      #conc = uds/area

      ## reporting 
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
      # doesn't work w mpi 
      #print "MPI mode", params.useMPI
      if params.useMPI==False:
        print "SDF"
        us.append(InterpolateData(mesh,u_n[idxCa]))

      # update
      if mode=="sachse":
        iryr.t = t 
        
        caTT= assemble(u_n[idxCa]*ds(rMarker,domain=mesh))
        vol = assemble(Constant(1.)*ds(rMarker,domain=mesh))
        caTT = caTT/vol 
        print "TT %f" % caTT
        iryrDefect.ca = caTT
      if mode=="ode":
        # TODO run ODE model 
        # aggregate state variables
        # need to assemble concs over entire vol 
        #def UpdateBeardVals(nC,bM,mgatpi_uM,fadpi_uM ,fmgi_uM ):
        jCa.r  = 0. # TODO update with ode 
      t += dt


  #
  results = empty()
  # not pickle-safe
  results.mesh = mesh 
  results.u_n = u_n  
  # pickle-safe
  results.concs = np.asarray(concs)
  if params.useMPI==False:
    results.us = us
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
    elif(arg=="-test2D"):
      doit(debug=False,mode="2D_SSL")    
      quit()
    elif(arg=="-testsatin"):
      doit(debug=False,mode="satin")    
      quit()
    elif(arg=="-test2"):
      params = Params()
      params.T = 20
      doit(debug=False,params=params,mode="sachse2TT")
      quit()
    elif(arg=="-test4"):
      params = Params()
      params.T = 20
      doit(debug=False,params=params,mode="sachse4TT")
      quit()
  





  raise RuntimeError("Arguments not understood")


