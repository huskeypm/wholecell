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
  D_CleftSSL = Constant(1.0)
  D_CleftCyto= Constant(1.0)
  dist = Constant(1.0) # distance between compartments [um] 


  # Init conditions
  cCaCleftInit = 0.1  # [uM] 
  cCaSSLInit = 1.0
  cCaInit = 10.

  volSSL = 100. #Constant(2.) # [um^3] *float(volCyto))
  volCleft = 1. #  Constant(4.) # [um^3] *float(volSSL))


def tsolve(outName="output.pvd",\
           mode="2D_SSL", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
           doAssert = True,  
           debug=False):

  #mesh = UnitCubeMesh(8,8,8)
  ssl=True
  if "2D" in mode:
    import sarcomere2DSSL
    cm = sarcomere2DSSL.Sarcomere2DSSL(params=params)

  if "noSSL" in mode:
    ssl=False

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
  D_CleftSSL = params.D_CleftSSL
  D_CleftCyto = params.D_CleftCyto
  dist = params.dist
  
  # iterator
  dt = Constant(1.0)  # [ms]  
  tstop = 1.  
  
  # volume/area info 
  volCyto = params.volCyto
  #volCyto = Constant(assemble(Constant(1.0)*dx()))
  #area_u = Constant(assemble(Constant(1.0)*ds(face_marker)))
  volSSL = params.volSSL
  volCleft= params.volCleft
  volFrac_CytoSSL= volCyto/volSSL
  volFrac_SSLCyto= 1/volFrac_CytoSSL
  volFrac_SSLCleft= volSSL/volCleft
  volFrac_CleftSSL= 1/volFrac_SSLCleft
  volFrac_CleftCyto= volCleft/volCyto
  volFrac_CytoCleft= 1/volFrac_CleftCyto           

# System
  ssl = False
  #ssl = True
  
  # System
  
  ## PDE 
  # Time derivative and diffusion of field species
  F = ((cCa_n-cCa_0)*vCa/dt + Du*inner(grad(cCa_n), grad(vCa)))*dx()
  
  if ssl:
    # Flux to mesh domain from scalar domain s 
    jFlux_SSLCyto = D_SSLCyto*(cCaSSL_n-cCa_n)/dist # note that this is multiplied by ds(face_marker) below
    F += -jFlux_SSLCyto*vCa*ds(face_marker)
  else: 
    # Flux to mesh domain from scalar domain s 
    jFlux = D_CleftCyto*(cCaCleft_n-cCa_n)/dist

    F += -jFlux*vCa*ds(face_marker)


  
  ## Scalar SSL (even if ssl=False, we still define a PDE just ot keep the code more simple) 
  F += (cCaSSL_n-cCaSSL_0)/dt*vCaSSL*dx() 
  
  if ssl:
    # flux from scalar domain s to mesh 
    F += jFlux_SSLCyto*volFrac_CytoSSL*vCaSSL*ds(face_marker)
  
    # flux to scalar domain r to s
    jFlux_CleftSSL =D_CleftSSL*(cCaCleft_n-cCaSSL_n)/dist 
    F += -jFlux_CleftSSL*volFrac_CytoSSL*vCaSSL*dx()
  else:
    1 # do nothing 
  
  
  ## Scalar cleft
  F += (cCaCleft_n-cCaCleft_0)/dt*vCaCleft*dx()
  if ssl:
    # flux between SSL/cleft
    F +=  jFlux_CleftSSL*volFrac_CytoCleft*vCaCleft*dx()
  else: 
    # flux between SSL/cleft
    F += volFrac_CytoCleft*jFlux*vCaCleft*ds(face_marker)


  
  # Init conditions
  U_0.interpolate(Constant((params.cCaInit, params.cCaSSLInit, params.cCaCleftInit)))
  
  def conservation():
      ## PKH test
      acCaCleft_n = assemble(cCaCleft_n*volFrac_CleftCyto*dx())
      #ccCaCleft_n = assemble(cCaCleft_n*volFrac_CleftCyto/volCleft*dx())
      # below is shorthand, but i left comment above so its clear where it originated 
      ccCaCleft_n = assemble(cCaCleft_n*volCyto*dx())
      acCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto*dx())
      #ccCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto/volSSL*dx())
      ccCaSSL_n = assemble(cCaSSL_n*volCyto*dx())
      acCa_n = assemble(cCa_n*dx())
      ccCa_n = assemble(cCa_n/volCyto*dx())
      conserved = acCaSSL_n+acCaCleft_n+acCa_n
    
      if MPI.rank(mpi_comm_world())==0:
        print "cCa_n " , acCa_n, "conc ", ccCa_n
        print "cCaSSL_n " , acCaSSL_n, "conc ", ccCaSSL_n
        print "cCaCleft_n " , acCaCleft_n, "conc ", ccCaCleft_n
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

      # report 
      print "############### ", t
      conserved_ti = conservation()

      ## store 
      file << (U_n.sub(idxCa),t)
      # store hdf 
      #if compartmentSSL:
      if 1: 
        uCa = project(cCa_n,V)
        hdf.write(uCa,"uCa",ctr)

        # see conservation routine above for why we multiply by volCyto 
        uCaSSL = project(cCaSSL_n,R)      
        uCaSSL.vector()[:]*=volCyto
        #print assemble(uCaSSL*dx())
        hdf.write(uCaSSL,"uCaSSL",ctr)

        uCaCleft = project(cCaCleft_n,R)         
        uCaCleft.vector()[:]*=volCyto
        #print assemble(uCaCleft*dx())
        hdf.write(uCaCleft,"uCaCleft",ctr)


  
      ## update 
      ctr+=1
      U_0.assign(U_n)
      t += float(dt)

  # close 
  hdf.close()
  
  
  # Check that the gods of transport are happy 
  if doAssert:
    message = "%f | %f " % (conserved_ti , conserved_t0)
    assert( abs(conserved_ti - conserved_t0) < 1e-3 ), message
    #assert( abs(acCaSSL_n+acCaCleft_n+acCa_n - 60000) < 1e-5)
    print "PASS"
  
  #plot(U_plot, interactive=True, range_min=min_value, range_max=max_value)

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



