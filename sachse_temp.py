#
# Time dependent solver, single species 
# 
# TODO 
# - get correct SAs, volumes for all compartments
# - check that Ca is conserved between compartments 
# - get concentrations of buffs, etc correct
# - add diffusivity/buffer as D' = D B/(1+KD)
#
# Check README/Fields for SIAM models section for details on implementation 
#
# Validations
# - Ca flows between compartments, both ways 
#
# Standard units
# time [ms]
# volume fluxes [uM/ms]
# currents [A/F]?

from dolfin import *
import numpy as np
import matplotlib.pylab as plt
class empty:pass 
import sys
sys.path.append("./siam/")

eps = 0.05

## Units
nm_to_um = 1.e-3
ttRad = 0.25 # [um]


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
  DCa_SSL = 0.6 * DCa


  # buffer (mostly TnC) 
  alphabuff=0.04 # kon verified [1/uMms] 
  Btot = 70. # [TnC] [uM] verified 
  betabuff = 0.04 # koff verfied [1/ms]
 
  # buffer (imaginary for SSL) 
  alphabuff_SSL=0.04 # kon verified [1/uMms] 
  Btot_SSL = 70. # [TnC] [uM] verified 
  betabuff_SSL = 0.04 # koff verfied [1/ms]
    

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
  cCaSSLinit = 0.5
  cCaCleftinit = 1.0
  cInits = [cCainit,cBuffinit,cFluoinit,cCaCleftinit,cCaSSLinit]  # cCaCleft=cCainit, cCaSSL=cCainit at start 


  # compartments 
  Dcomp = 1.   # [um^2/ms] 
  dist = 1e-2 # [um]
  volume_SSL = 0.50 # [um^3]
  volume_Cleft = 0.70 # [um^3]

  # debug
  jConstant = 1 # [uM/ms]


 # TODO put in real info 


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
def tsolve(outName="output.pvd",\
           mode="sachse2TT", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
	   debug=False): 
 
  ## Create mesh 
  if debug: 
    params.T = 10
    mesh = UnitCubeMesh(5,5,5)

  if "2D" in mode:
    if mode=="2D_SSL":
      import sarcomere2DSSL
      cm = sarcomere2DSSL.Sarcomere2DSSL(params=params)
          
    if mode=="2D_noSSL":
        # p/rxn for SSL 
      import sarcomere2DwoSSL
      cm = sarcomere2DwoSSL.sarcomere2DwoSSL(params=params)

  elif "sachse" in mode: 
    if mode=="sachse4TT":
      import sarcomere4TT
      cm = sarcomere4TT.sarcomere4TT(params=params)
    else: 
      import sarcomere2TT
      cm = sarcomere2TT.sarcomere2TT(params=params)

  elif "satin" in mode:
    import sarcomereSatin
    cm = sarcomereSatin.sarcomereSatin(params=params)
    
  cm.params.nDOF = cm.nDOF # goes in master class 
  mesh = cm.GetMesh()
  dim =mesh.ufl_cell().geometric_dimension()
  cm.CalcGeomAttributes()
  print "MOVE ME" 
  params.volumeCytosol = cm.volume
  volume_fracCS = params.volumeCytosol/params.volume_SSL
  volume_fracCC = params.volumeCytosol/params.volume_Cleft   


  [idxCa,idxBuff,idxFluo,idxCaCleft,idxCaSSL]=cm.GetIndices()
  idxCa = 0 # PDE 
  idxBuff = 1
  idxFluo = 2
  idxCaCleft = 3
  idxCaSSL = 4  # not always used

  # function spaces
  V = FunctionSpace(mesh, "CG", 1)
  R = FunctionSpace(mesh,"R",0) # for each compartment 
#######
  ## Decide whether SSL is an explicit compartment or lumped into cyto 
  compartmentSSL = False
  if cm.nDOF == 4: 
    compartmentSSL = False
  elif cm.nDOF == 5: 
    compartmentSSL = True   
  else:
    raise RuntimeError("Not prepared to handle this") 




  Vs = []
  for i in range(cm.nDOF_Fields):
        #print i 
        Vs.append(V)
  for i in range(cm.nDOF_Scalars):
        #print i 
        Vs.append(R)
  ME = MixedFunctionSpace(Vs) # Ca, Buffered Ca, Fluo Ca

  # Define trial and test functions
  du    = TrialFunction(ME)
  if compartmentSSL: 
    vCa,vBuff,vFluo, vCaCleft,vCaSSL  = TestFunctions(ME)
    vs = [vCa,vBuff,vFluo,vCaCleft,vCaSSL]
  else:
    vCa,vBuff,vFluo, vCaCleft  = TestFunctions(ME)
    vs = [vCa,vBuff,vFluo,vCaCleft]

  # Define functions
  u_n   = Function(ME)  # current solution
  u_0  = Function(ME)  # solution from previous converged step

  # split mixed functions
  #print "WARNING: not sure how to generalize"       
  if compartmentSSL: 
    cCa_n,cBuff_n, cFluo_n,cCaCleft_n,cCaSSL_n = split(u_n)
    cns = [cCa_n,cBuff_n, cFluo_n,cCaCleft_n,cCaSSL_n]
    cCa_0,cBuff_0, cFluo_0, cCaCleft_0,cCaSSL_0 = split(u_0)
    c0s = [cCa_0,cBuff_0, cFluo_0, cCaCleft_0,cCaSSL_0]
  else: 
    cCa_n,cBuff_n, cFluo_n,cCaCleft_n = split(u_n)
    cns = [cCa_n,cBuff_n, cFluo_n,cCaCleft_n]
    cCa_0,cBuff_0, cFluo_0, cCaCleft_0 = split(u_0)
    c0s = [cCa_0,cBuff_0, cFluo_0, cCaCleft_0]

  ## mark boundaries 
  # problem specific          
  subdomains = MeshFunction("size_t",mesh,dim-1)
  subdomains.set_all(0)
  lMarker,rMarker,slMarker = cm.Boundaries(subdomains)

  ## decide on source of probability density  
  init_cond = cm.InitialConditions()
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
  for i in range(cm.nDOF_Fields):
    #RHSi= -inner(Ds[i]*grad(cns[i]), grad(vs[i]))*dx
    RHSi= -inner(Ds[i]*grad(u_n[i]), grad(vs[i]))*dx
    RHSis.append(RHSi)

  # no diffusion in scalar domain 
  for i in range(cm.nDOF_Scalars):
    j = i + cm.nDOF_Fields
    RHSi = Constant(0)*vs[j]*dx # for consistency
    RHSis.append(RHSi)


  ## buffreactions
    
  ### LHS    
  dt = params.dt
  ## Cytosol (PDE) domain 
  #LHS = (cCa_n-cCa_0)*vCa/dt * dx - RHSCa
  ##LHSis=[]
  L=0
  for i in range(cm.nDOF_Fields):
    # TODO DBL check 
    #LHSi= (cns[i]-c0s[i])*vs[i]/dt * dx - RHSis[i]     
    LHSi= (u_n[i]-u_0[i])*vs[i]/dt * dx - RHSis[i]     
    #LHSis.append(RHSi)
    L+= LHSi


  ## Cleft/SSL domains
  print "WARNING: need to use volfracCS for SSl, CC for Cleft"
  for i in range(cm.nDOF_Scalars):
    j = i + cm.nDOF_Fields
    LHSi = (u_n[j]-u_0[j])*vs[j]/(dt*volume_fracCS) * dx - RHSis[j]
    L+= LHSi

  ## Ca fluxes between compartments (no buffers)
  # verified flux is working 
  Dcomp = Constant(params.Dcomp*0.)
  dist  = params.dist 
  #if compartmentSSL: 
  #  # Flux from cleft to SSL - verified    
  #  L -= Dcomp*(cCaCleft_n - cCaSSL_n)*vCaCleft/dist*dx() #ds(markerCS) 
  #  L += Dcomp*(cCaCleft_n - cCaSSL_n)*vCaSSL/dist*dx() #ds(markerCS) 
  #  # Flux from SSL to cytosol 
  #  L -= Dcomp*(cCaSSL_n - cCa_n)*vCa/dist*ds(lMarker) 
  #  L += Dcomp*(cCaSSL_n - cCa_n)*vCaSSL/dist*ds(lMarker) 
  #  1
  #else:
  #  # Flux from Cleft to cytosol - verified (not calibrated)  
  #  L -= Dcomp*(cCaCleft_n - cCa_n)*vCa/dist*ds(lMarker) 
  #  L += Dcomp*(cCaCleft_n - cCa_n)*vCaCleft/dist*ds(lMarker) 
  #  1
#

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

  ## File IO
  ctr=0
  file << (u_n.sub(idxCa),t) 
  Hdf=HDF5File(mesh.mpi_comm(), hdfName, "w")
  Hdf.write(mesh, "mesh")
  # temp hack for storing volume 
  x = Function(V)
  x.vector()[0] = params.volume_SSL
  Hdf.write(x,"volume_SSL")
  x.vector()[0] = params.volume_Cleft
  Hdf.write(x,"volume_Cleft")
  params.volume_Cyto = assemble(Constant(1.)*dx(domain=mesh))
  x.vector()[0] = params.volume_Cyto
  Hdf.write(x,"volume_Cyto")
  print params.volume_SSL
  

  ## Loop 
  while (t < params.T):
      # advance 
      t0=t
      u_0.vector()[:] = u_n.vector()
      solver.solve(problem,u_n.vector())

      # store 
      file << (u_n.sub(idxCa),t) 

      
      # store

      # report on prev iter
      #uds = assemble(u0p*ds(rMarker,domain=mesh))#,mesh=mesh)
      #area = assemble(Constant(1.)*ds(rMarker,domain=mesh))#,mesh=mesh)
      #conc = uds/area

      ## reporting 
      vol = assemble(Constant(1.)*dx(domain=mesh))
      assembles=np.zeros(cm.nDOF)
      for i in range(cm.nDOF):
        #assembles[i] = assemble(cns[i]*dx(domain=mesh))
        assembles[i] = assemble(u_n[i]*dx(domain=mesh))

      concis = assembles
      for i in range(cm.nDOF_Fields):
        concis[i]  = assembles[i]/params.volume_Cyto
      #print "WARNING: must use correct volume" 
      #for i in range(cm.nDOF_Scalars):
      #  j = i = cm.nDOF_Fields
      if compartmentSSL: 
        concis[idxCaCleft]  = assembles[idxCaCleft]/params.volume_Cleft
        #print "ass ", assembles[idxCaCleft]
        concis[idxCaSSL]  = assembles[idxCaSSL]/params.volume_SSL
      else: 
        concis[idxCaCleft]  = assembles[idxCaCleft]/params.volume_Cleft

      if MPI.rank(mpi_comm_world())==0:
        if compartmentSSL: 
          print "Conc CaCleft", concis[idxCaCleft]
          print "Conc CaSSL", concis[idxCaSSL]
          print "Conc Ca", concis[idxCa]
        else:
          print "Conc CaCleft", concis[idxCaCleft]
          print "Conc Ca", concis[idxCa]
 
 
      ts.append(t0)
      concs.append(concis)
      # doesn't work w mpi 
      #print "MPI mode", params.useMPI
      #if params.useMPI==False:
      #  print "SDF"
      #  #us.append(InterpolateData(mesh,u_n[idxCa]))
 
      # store hdf 
      if compartmentSSL:
        uCa = project(u_n[idxCa],V)
        Hdf.write(uCa,"uCa",ctr)
        uCaSSL = project(u_n[idxCaSSL],R)
        Hdf.write(uCaSSL,"uCaSSL",ctr)
        uCaCleft = project(u_n[idxCaCleft],R)
        #print "ass2 ", assemble(uCaCleft*dx)
        Hdf.write(uCaCleft,"uCaCleft",ctr)
      else: 
        uCa = project(u_n[idxCa],FunctionSpace(mesh,"CG",1))
        Hdf.write(uCa,"uCa",ctr)
        uCaCleft = project(u_n[idxCa],FunctionSpace(mesh,"R",0)) 
        Hdf.write(uCaCleft,"uCaCleft",ctr)

      ctr+=1

      ## update
      if 0: # mode=="sachse":
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
  ##  End loop 


  #
  Hdf.close()
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

def doit(debug=True,mode="",params=Params(),hdfName="a.h5"):  
  tsolve(debug=debug,mode=mode,params=params,hdfName=hdfName)
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
  params = Params()
  params.T = 10 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-debug"):
      #arg1=sys.argv[i+1] 
      doit(debug=True)     
    elif(arg=="-test0"):
      doit(debug=False)    
      quit()
    elif(arg=="-test1"):
      tag = "ode"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2D"):
      tag = "2D_SSL"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2Dno"):
      tag = "2D_noSSL"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-testsatin"):
      tag = "satin"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2"):
      tag = "sachse2TT"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test4"):
      tag = "sachse4TT"
      doit(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()






  raise RuntimeError("Arguments not understood")


