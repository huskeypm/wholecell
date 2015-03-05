"""
"""
import sys
sys.path.append("./siam/")
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

  D_SSLCyto = Constant(1.0) # [um^2/ms]
  D_ClftSSL = Constant(1.0)
  dist = Constant(1.0) # distance between compartments [um] 


  # Init conditions
  cCaCleftInit = 0.1  # [uM] 
  cCaSSLInit = 1.0
  cCaInit = 10.

  volSSL = 2. #Constant(2.) # [um^3] *float(volCyto))
  volCleft = 4. #  Constant(4.) # [um^3] *float(volSSL))


def tsolve(outName="output.pvd",\
           mode="2D_SSL", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
           doAssert = True,  
           debug=False):

  #mesh = UnitCubeMesh(8,8,8)
  if "2D" in mode:
    if mode=="2D_SSL":
      import sarcomere2DSSL
      cm = sarcomere2DSSL.Sarcomere2DSSL(params=params)

    if mode=="2D_noSSL":
      raise RuntimeError("not tet implemented") 
        # p/rxn for SSL 
      import sarcomere2DwoSSL
      cm = sarcomere2DwoSSL.sarcomere2DwoSSL(params=params)


  cm.params.nDOF = cm.nDOF # goes in master class 
  mesh = cm.GetMesh()
  dims =mesh.ufl_cell().geometric_dimension()
  cm.CalcGeomAttributes()
  params.volCyto= cm.volume

  ## 
  print 'WARNING: will need to updated these'
  idxCa = 0 # PDE 
  #idxBuff = 1
  #idxFluo = 2
  idxCaSSL = 1
  idxCaCleft = 2  # not always used


  # Boundaries and measures 
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
  
  vCa, vCaSSL, vCaCleft = TestFunctions(ME)
  
  U_n = Function(ME)
  U_0 = Function(ME)
  cCa_n, cCaSSL_n, cCaCleft_n = split(U_n)
  cCa_0, cCaSSL_0, cCaCleft_0 = split(U_0)
  
  # Params
  mode = "boundaryFlux" 
  Du = Constant(1.0)
  
  # intercompartment transport 
  D_SSLCyto = params.D_SSLCyto
  D_ClftSSL = params.D_ClftSSL
  dist = params.dist
  
  # iterator
  dt = Constant(1.0)  # [ms]  
  tstop = 20. 
  
  # volume/area info 
  volCyto = params.volCyto
  #volCyto = Constant(assemble(Constant(1.0)*dx()))
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
  jFlux_SSLCyto = D_SSLCyto*(cCaSSL_n-cCa_n)/dist # note that this is multiplied by ds(face_marker) below
  F += -jFlux_SSLCyto*vCa*ds(face_marker)
  
  ## Scalar s
  # Time derivative of scalar species 
  F += (cCaSSL_n-cCaSSL_0)/dt*vCaSSL*dx() 
  
  # flux from scalar domain s to mesh 
  F += jFlux_SSLCyto*volFrac_CytoSSL*vCaSSL*ds(face_marker)
  
  # flux to scalar domain r to s
  jFlux_ClftSSL =D_ClftSSL*(cCaCleft_n-cCaSSL_n)/dist 
  F += -jFlux_ClftSSL*volFrac_CytoSSL*vCaSSL*dx()
  
  
  ## Scalar r
  F += (cCaCleft_n-cCaCleft_0)/dt*vCaCleft*dx()
  F +=  jFlux_ClftSSL*volFrac_CytoClft*vCaCleft*dx()
  
  # Init conditions
  U_0.interpolate(Constant((params.cCaInit, params.cCaSSLInit, params.cCaCleftInit)))
  
  def conservation():
      ## PKH test
      acCaCleft_n = assemble(cCaCleft_n*volFrac_ClftCyto*dx())
      ccCaCleft_n = assemble(cCaCleft_n*volFrac_ClftCyto/volCleft*dx())
      acCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto*dx())
      ccCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto/volSSL*dx())
      acCa_n = assemble(cCa_n*dx())
      ccCa_n = assemble(cCa_n/volCyto*dx())
      print "cCa_n " , acCa_n, "conc ", ccCa_n
      print "cCaSSL_n " , acCaSSL_n, "conc ", ccCaSSL_n
      print "cCaCleft_n " , acCaCleft_n, "conc ", ccCaCleft_n
      conserved = acCaSSL_n+acCaCleft_n+acCa_n
      print "CONSERVATION:", conserved
      return conserved

  
  
      
  U_n.vector()[:] = U_0.vector()[:]
  conserved_t0 = conservation()

  ## File IO
  t=0.
  ctr=0
  file = File("output.pvd", "compressed")
  file << (U_n.sub(idxCa),t)     
  hdf=HDF5File(mesh.mpi_comm(), hdfName, "w")
  hdf.write(mesh, "mesh")
  # temp hack for storing volume 
  x = Function(V)
  x.vector()[0] = params.volSSL
  hdf.write(x,"volSSL")
  x.vector()[0] = params.volCleft
  hdf.write(x,"volCleft")
  params.volCyto = assemble(Constant(1.)*dx(domain=mesh))
  x.vector()[0] = params.volCyto
  hdf.write(x,"volCyto")
  
  
  
  while t < tstop+DOLFIN_EPS:
      ## solve 
      solve(F==0, U_n)

      ## store 
      file << (U_n.sub(idxCa),t)
      # store hdf 
      #if compartmentSSL:
      if 1: 
        uCa = project(U_n[idxCa],V)
        hdf.write(uCa,"uCa",ctr)
        uCaSSL = project(U_n[idxCaSSL],R)
        hdf.write(uCaSSL,"uCaSSL",ctr)
        uCaCleft = project(U_n[idxCaCleft],R)
        #print "ass2 ", assemble(uCaCleft*dx)
        hdf.write(uCaCleft,"uCaCleft",ctr)

      # report 
      print "############### ", t
      conserved_ti = conservation()

  
      ## update 
      ctr+=1
      U_0.assign(U_n)
      t += float(dt)

  # close 
  hdf.close()
  
  
  if doAssert:
    message = "%f | %f " % (conserved_ti , conserved_t0)
    assert( abs(conserved_ti - conserved_t0) < 1e-3 ), message
    #assert( abs(acCaSSL_n+acCaCleft_n+acCa_n - 60000) < 1e-5)
    print "PASS"
  
  #plot(U_plot, interactive=True, range_min=min_value, range_max=max_value)

tsolve()
