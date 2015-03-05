"""
Diffusion between three compartments 

Verified for volumetric 'fluxes'
"""
from dolfin import *

class Boundary(SubDomain):
  def inside(self,x,on_boundary):
    left= x[0] <DOLFIN_EPS
    #print left 
    return left and on_boundary

class Params():
  useMPI = True
  paraview = True
  verbose=False

  Dsu = Constant(1.0) # [um^2/ms]
  Drs = Constant(1.0)
  dist = Constant(1.0) # distance between compartments [um] 


  # Init conditions
  cCaSSLInit = 0.1  # [uM] 
  cCaCleftInit = 1.0
  cCaInit = 10.

  volSSL = Constant(2.) # [um^3] *float(volCyto))
  volCleft = Constant(4.) # [um^3] *float(volSSL))


def tsolve(outName="output.pvd",\
           mode="sachse2TT", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
           debug=False):

  #mesh = UnitCubeMesh(8,8,8)
  mesh = UnitSquareMesh(5,5) # works for 2D too 
  dims = mesh.ufl_cell().geometric_dimension()



  subdomains = MeshFunction("size_t",mesh,dims-1)
  face_marker = 10
  boundary = Boundary()
  boundary.mark(subdomains,face_marker)
  ds = Measure("ds")[subdomains]
  dx = Measure("dx")
  ds = ds(domain=mesh,subdomain_data=subdomains)
  dx = dx(domain=mesh)

  # FiniteElements
  V = FunctionSpace(mesh, "CG", 1)
  R = FunctionSpace(mesh, "R", 0)
  ME = MixedFunctionSpace([V,R,R])
  
  vCa, vCaCleft, vCaSSL = TestFunctions(ME)
  
  U_n = Function(ME)
  U_0 = Function(ME)
  cCa_n, cCaCleft_n, cCaSSL_n = split(U_n)
  cCa_0, cCaCleft_0, cCaSSL_0 = split(U_0)
  
  # Params
  mode = "boundaryFlux" 
  Du = Constant(1.0)
  
  # intercompartment transport 
  Dsu = params.Dsu
  Drs = params.Drs
  dist = params.dist
  
  # iterator
  dt = Constant(1.0)  # [ms]  
  tstop = 20. 
  
  # volume/area info 
  volCyto = Constant(assemble(Constant(1.0)*dx()))
  #area_u = Constant(assemble(Constant(1.0)*ds(face_marker)))
  volSSL = params.volSSL
  volCleft= params.volCleft
  volFrac_CytoSSL= volCyto/volSSL
  volFrac_SSLCyto= 1/volFrac_CytoSSL
  volFrac_SSLClft= volSSL/volCleft
  volFrac_ClftSSL= 1/volFrac_SSLClft
  volFrac_ClftCyto= volCleft/volCyto
  volFrac_CytoClft= 1/volFrac_ClftCyto           
  
  
  # System
  
  ## PDE 
  # Time derivative and diffusion of field species
  F = ((cCa_n-cCa_0)*vCa/dt + Du*inner(grad(cCa_n), grad(vCa)))*dx()
  
  # Flux to mesh domain from scalar domain s 
  jFlux_SSLCyto = Dsu*(cCaCleft_n-cCa_n)/dist # note that this is multiplied by ds(face_marker) below
  F += -jFlux_SSLCyto*vCa*ds(face_marker)
  
  ## Scalar s
  # Time derivative of scalar species 
  F += (cCaCleft_n-cCaCleft_0)/dt*vCaCleft*dx() 
  
  # flux from scalar domain s to mesh 
  F += jFlux_SSLCyto*volFrac_CytoSSL*vCaCleft*ds(face_marker)
  
  # flux to scalar domain r to s
  jFlux_ClftSSL =Drs*(cCaSSL_n-cCaCleft_n)/dist 
  F += -jFlux_ClftSSL*volFrac_CytoSSL*vCaCleft*dx()
  
  
  ## Scalar r
  F += (cCaSSL_n-cCaSSL_0)/dt*vCaSSL*dx()
  F +=  jFlux_ClftSSL*volFrac_CytoClft*vCaSSL*dx()
  
  # Init conditions
  U_0.interpolate(Constant((params.cCaInit, params.cCaCleftInit, params.cCaSSLInit)))
  
  def conservation():
      ## PKH test
      acCaSSL_n = assemble(cCaSSL_n*volFrac_ClftCyto*dx())
      ccCaSSL_n = assemble(cCaSSL_n*volFrac_ClftCyto/volCleft*dx())
      acCaCleft_n = assemble(cCaCleft_n*volFrac_SSLCyto*dx())
      ccCaCleft_n = assemble(cCaCleft_n*volFrac_SSLCyto/volSSL*dx())
      acCa_n = assemble(cCa_n*dx())
      ccCa_n = assemble(cCa_n/volCyto*dx())
      print "cCa_n " , acCa_n, "conc ", ccCa_n
      print "cCaCleft_n " , acCaCleft_n, "conc ", ccCaCleft_n
      print "cCaSSL_n " , acCaSSL_n, "conc ", ccCaSSL_n
      conserved = acCaCleft_n+acCaSSL_n+acCa_n
      print "CONSERVATION:", conserved
      return conserved

  
  
      
  U_n.vector()[:] = U_0.vector()[:]
  conserved_t0 = conservation()
  
  
  
  t=0
  U_plot = Function(V)
  #plot(U_plot, interactive=True, range_min=min_value, range_max=max_value)
  while t < tstop+DOLFIN_EPS:
      solve(F==0, U_n)
      U_plot.assign(U_n.sub(0, deepcopy=True))
      #plot(U_plot, interactive=False, range_min=min_value, range_max=max_value)
      #print "CONSERVATION:", assemble(cCa_n*dx())+assemble(cCaCleft_n/volFrac_CytoSSL*dx())
      #print "CONSERVATION:", assemble(cCa_n*dx())+assemble(cCaCleft_n*volFrac_SSLCyto*dx())
      print "############### ", t
      conserved_ti = conservation()
  
      U_0.assign(U_n)
      t += float(dt)
  
  
  message = "%f | %f " % (conserved_ti , conserved_t0)
  assert( abs(conserved_ti - conserved_t0) < 1e-3 ), message
  #assert( abs(acCaCleft_n+acCaSSL_n+acCa_n - 60000) < 1e-5)
  print "PASS"
  
  #plot(U_plot, interactive=True, range_min=min_value, range_max=max_value)

tsolve()
