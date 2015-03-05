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

if 1: 
  #mesh = UnitCubeMesh(8,8,8)
  mesh = UnitSquareMesh(5,5) # works for 2D too 
  dims = mesh.ufl_cell().geometric_dimension()
  subdomains = MeshFunction("size_t",mesh,dims-1)
  face_marker = 10
  boundary = Boundary()
  boundary.mark(subdomains,face_marker)
  ds = ds(domain=mesh,subdomain_data=subdomains)
  dx = dx(domain=mesh)

# FiniteElements
V = FunctionSpace(mesh, "CG", 1)
R = FunctionSpace(mesh, "R", 0)
ME = MixedFunctionSpace([V,R,R])

v, p, n = TestFunctions(ME)

U_n = Function(ME)
U_0 = Function(ME)
u_n, s_n, r_n = split(U_n)
u_0, s_0, r_0 = split(U_0)

# Params
mode = "volumeSource"
mode = "boundaryFlux" 
Du = Constant(1.0)

# intercompartment transport 
if mode == "boundaryFlux":
  Dsu = Constant(1.0)
  Drs = Constant(1.0)
  dist = Constant(1.0) # distance between compartments 
elif mode == "volumeSource":
  jv_us = 0.0 # volume 'flux' rate into u from s [uM/ms] (e.g. applied over dx) 
  jv_sr = 0.0 # flux rate into s from r [uM/ms] 

# iterator
dt = Constant(1.0)   
tstop = 20. 
# Init conditions
scalar_w = 0.1
scalar_0 = 1.0
field_0 = 10.

# volume/area info 
volume_u = Constant(assemble(Constant(1.0)*dx()))
area_u = Constant(assemble(Constant(1.0)*ds(face_marker)))
volume_s = Constant(2.*float(volume_u))
volume_r = Constant(4.*float(volume_s))
volFrac_us= volume_u/volume_s
volFrac_su= 1/volFrac_us
volFrac_sr= volume_s/volume_r
volFrac_rs= 1/volFrac_sr
volFrac_ru= volume_r/volume_u
volFrac_ur= 1/volFrac_ru           


# System

## PDE 
# Time derivative and diffusion of field species
F = ((u_n-u_0)*v/dt + Du*inner(grad(u_n), grad(v)))*dx()
if mode=="volumeSource":
  F+= -jv_us*v*dx()

# Flux to mesh domain from scalar domain s 
jFlux_su = Dsu*(s_n-u_n)/dist # note that this is multiplied by ds(face_marker) below
if mode=="boundaryFlux":
  F += -jFlux_su*v*ds(face_marker)

## Scalar s
# Time derivative of scalar species 
F += (s_n-s_0)/dt*p*dx() 

# flux from scalar domain s to mesh 
if mode=="volumeSource":
  F += jv_us*volFrac_us* p*dx()
  F += -jv_sr*volFrac_us* p*dx()

if mode=="boundaryFlux":
  F += jFlux_su*volFrac_us*p*ds(face_marker)

# flux to scalar domain r to s
jFlux_rs =Drs*(r_n-s_n)/dist 
if mode=="boundaryFlux":
  F += -jFlux_rs*volFrac_us*p*dx()


## Scalar r
F += (r_n-r_0)/dt*n*dx()
if mode=="volumeSource":
  F +=  jv_sr*volFrac_ur* n*dx()
if mode=="boundaryFlux":
  F +=  jFlux_rs*volFrac_ur*n*dx()

# Init conditions
U_0.interpolate(Constant((field_0, scalar_0, scalar_w)))

def conservation():
    ## PKH test
    ar_n = assemble(r_n*volFrac_ru*dx())
    cr_n = assemble(r_n*volFrac_ru/volume_r*dx())
    as_n = assemble(s_n*volFrac_su*dx())
    cs_n = assemble(s_n*volFrac_su/volume_s*dx())
    au_n = assemble(u_n*dx())
    cu_n = assemble(u_n/volume_u*dx())
    print "u_n " , au_n, "conc ", cu_n
    print "s_n " , as_n, "conc ", cs_n
    print "r_n " , ar_n, "conc ", cr_n
    conserved = as_n+ar_n+au_n
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
    #print "CONSERVATION:", assemble(u_n*dx())+assemble(s_n/volFrac_us*dx())
    #print "CONSERVATION:", assemble(u_n*dx())+assemble(s_n*volFrac_su*dx())
    print "############### ", t
    conserved_ti = conservation()

    U_0.assign(U_n)
    t += float(dt)


message = "%f | %f " % (conserved_ti , conserved_t0)
assert( abs(conserved_ti - conserved_t0) < 1e-3 ), message
#assert( abs(as_n+ar_n+au_n - 60000) < 1e-5)
print "PASS"

#plot(U_plot, interactive=True, range_min=min_value, range_max=max_value)
