#
# Time dependent solver, single species 
# 
# TODO 
# - get concentrations of buffs, etc correct
# - add diffusivity/buffer as D' = D B/(1+KD)
# Merge in sachse_temp, sachse_REPLACEMENT?
#- meaningful iRyr
#- defining SSL region within cyto #
# Check README/Fields for SIAM models section for details on implementation 
#
# Validations
# - see validations function 
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

class Boundary(SubDomain):
  def inside(self,x,on_boundary):
    left= x[0] <DOLFIN_EPS
    #print left 
    return left and on_boundary

class Params():
  #useMPI = True
  #paraview = True
  verbose=False

  T = 3.  # Total simulation time [ms]
  dt = 1. # time-step size [ms] 

  D_SSLCyto = Constant(1.0) # Diffusion rate between SSL/Cyto compartments (if SSL exists) [um^2/ms]
  D_CleftSSL = Constant(1.0 ) #  Diffusion rate between Cleft/SSL compartments (if SSL exists)
  D_CleftCyto= Constant(1.0) #  Diffusion rate between Cleft/Cyto compartments (if SSL does not exist)
  dist = Constant(1.0) # distance between compartments [um] 


#  # Init conditions
#  cCaCleftInit = 0.1  # Initial Ca concentration in Cleft compartment [uM] 
#  cCaSSLInit = 1.0
#  cCaInit = 10.

  # diffusion params 
  #DAb   = Dbulk# [um^2/ms] Diff const within PDE (domain 1) 
  DCa = 0.39  # [um^2/ms] verified
  DBuff = Constant(0.) # TnC
  DFluo = 0.1 # [um^2/ms] verified 
  Ds = [DCa,DBuff,DFluo]
  DCa_SSL = 0.6 * DCa

  # geom 
  volSSL = 1. # Volume of SSL domain [um^3]
  volCleft = 0.1 # [um^3] 

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
  cCaInit =0.1  # Initial Ca concentration in Cleft compartment [uM] 
  cBuffInit = buffered(Btot,cCaInit,alphabuff,betabuff)  # TnC unverified [uM]
  cFluoInit = buffered(Ftot,cCaInit,alphafluo,betafluo)  # [uM]
  cCaSSLInit = cCaInit
  cCaCleftInit = cCaInit
  cInits = [cCaInit,cBuffInit,cFluoInit,cCaCleftInit,cCaSSLInit]  # cCaCleft=cCainit, cCaSSL=cCainit at start 


  jTest = 0.1 # Generic flux applied to cytosolic domain [uM/ms]


# J [uM/ms] --> j[uM*um/ms]
# volDomain - domain in which flux is applied
# areaDomainFlux - area over which flux is applied 
def wholeCellJ_to_jSurf(wholeCellJ,volDomain,areaDomainFlux):
  return wholeCellJ*volDomain/areaDomainFlux

# Converts whole cell current I [A/F] into surface flux 
# For rabbit (C/V=4.5), 1 [A/F]  --> ~0.02 [uM/ms]
def wholeCellI_to_j(wholeCellI):
  # from flux notes on sharelatex
  nIons = 2 # number of ions per atom (2 for Ca, +/- 1 for NaCa2+ exchange)
  Fc = 9.6e-2 # Faraday [umol/C]
  Cm_o_Vm = 4.5 # Capacitance over cell volume [Farad/Liter]
  s_to_ms = 1e3
  
  jWholeCell = 1/(nIons*Fc) * Cm_o_Vm * wholeCellI * (1/s_to_ms)
  return jWholeCell
  

def wholeCellJ_to_jCleft(wholeCellJ, params):
  cleftRescale = (params.volCleft + params.volSSL + params.volCyto)
  cleftRescale/=params.volCleft
  return wholeCellJ*cleftRescale # [uM/ms]

def validation_Conversions():
  def dotest(cCaFinal,cCaInit):
    refConcChange = params.T*params.jTest # [uM]
    concChange = cCaFinal - cCaInit
    msg = "%f != %f " % ( concChange,refConcChange)
    assert( abs(concChange - refConcChange) < 1e-5 ), msg
    print "Passed flux test"

  idxCa = 0
  idxCaCleft = 4 # NEED TO NOT HARDCODE

  ## Validate A/F conversion 
  # apply flux to cytosol domain 
  print "WARNING: need to enable" 
  iFlux = 1 # [A/F]
  params = Params()
  params.jTest = wholeCellI_to_j(iFlux) 
  params.jTest*=4  # to make about 0.1 uM/ms
  params.T = 10
  params.D_SSLCyto = 0
  params.D_CleftSSL = 0
  params.D_CleftCyto = 0
  concsFinal = tsolve(doAssert="fluxTest_jSurfCyto", params=params)    
  cCaFinal = concsFinal[idxCa] # for cCa
  dotest(cCaFinal,params.cCaInit) 

  volCyto = params.volCyto # store here so dont have to load model

  ## apply to cleft
  # goal here is to find a 10 ms flux that raises the -entire-
  # cell concentration by 1.0 uM
  # 1) We first pump up the cleft for 10 ms, 
  # 2) We use the elevated CaCleft as an init conc, then show the total conc is correct in all compartments 
  ## 1) 
  params = Params()
  params.T = 10
  params.dt = 1 # ms 
  params.D_SSLCyto = 0. 
  params.D_CleftSSL = 0. 
  params.D_CleftCyto = 0. 
  chgConc = 0.5 # [uM] 
  params.jTest = chgConc/(params.dt*params.T) # [uM/ms]
  concsFinal= tsolve(doAssert="fluxTest_jCleft", params=params)    
  
  caCleftInit = concsFinal[idxCaCleft]
 
  ## 2) 
  params.T = 100
  params.dt = 10
  params.cInits[idxCaCleft] = caCleftInit
  #params.jTest = 0.
  params.D_SSLCyto = 1e2
  params.D_CleftSSL = 1e2
  params.D_CleftCyto = 1e2
  concsFinal= tsolve(doAssert="conservation", params=params)    
 
  # check that the change in conc. we anticipated is reflected
  # in cytosol
  assert(abs((params.cCaInit + chgConc) - concsFinal[idxCa]) < 1e-4)
  print "Pass!"
  quit()
  #params.jTest = wholeCellI_to_j(iFlux) 
  
  concsFinal= tsolve(doAssert="fluxTest_jSurfAllCompartments", params=params)    
  cCaCleftFinal = concsFinal[4] # for cCaCleft
  dotest(cCaCleftFinal) 
  assert(1==0), "Not done" 


  ##Validate Phi conversion
  params = Params()
  params.T = 10
  params.D_SSLCyto = 0
  params.D_CleftSSL = 0 
  params.D_CleftCyto = 0
  params.cCaInit = 11.

  cCaFinal = tsolve(doAssert="fluxTest_jSurfCyto", params=params)    
  dotest(cCaFinal) 
  
def validation():
  validation_Conversions()

  ## Reactions
  params = Params()
  params.D_SSLCyto = 0
  params.D_CleftSSL = 0 
  params.D_CleftCyto = 0
  tsolve(doAssert="conservation",params=params, buffers = True)
  
  assert(1==0), "WARNING: NOT FULLY VALIDATED"

  quit()
  
  ## conservation 
  tsolve(doAssert="conservation")
  print "Passed conservation test"

def tsolve(pvdName="output.pvd",\
           mode="2D_SSL", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
           doAssert = None,  
           reactions = None,
           buffers = None,
           debug=False):

  #mesh = UnitCubeMesh(8,8,8)
  ssl=True
  if "2D" in mode:
    import sarcomere2DSSL
    cm = sarcomere2DSSL.Sarcomere2DSSL(params=params)
    if "noSSL" in mode:
      cm.ssl = True
      ssl=False

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

  if ssl!=True:
    params.volSSL = 0.   # needed for validaton stuff
    
  cm.Init()
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
  # problem specific          
  subdomains = MeshFunction("size_t",mesh,dims-1)
  subdomains.set_all(0)
  lMarker,rMarker,slMarker = cm.Boundaries(subdomains)

  ds = Measure("ds")[subdomains]
  dx = Measure("dx")
  ds = ds(domain=mesh,subdomain_data=subdomains)
  dx = dx(domain=mesh)

  # FiniteElements
  V = FunctionSpace(mesh, "CG", 1)
  R = FunctionSpace(mesh, "R", 0)
  ME = MixedFunctionSpace([V,V,V,R,R])
  
  vCa, vCaBuff, vCaFluo, vCaSSL, vCaCleft = TestFunctions(ME)
  
  U_n = Function(ME)
  U_0 = Function(ME)
  cCa_n, cCaBuff_n, cCaFluo_n, cCaSSL_n, cCaCleft_n = split(U_n)
  cCa_0, cCaBuff_0, cCaFluo_0, cCaSSL_0, cCaCleft_0 = split(U_0)
  
  # Params
  mode = "boundaryFlux" 
  Du = Constant(1.0)
  
  # intercompartment transport 
  D_SSLCyto = params.D_SSLCyto
  D_CleftSSL = params.D_CleftSSL
  D_CleftCyto = params.D_CleftCyto
  dist = params.dist
  
  # iterator
  dt = params.dt 
  tstop = params.T
  
  # volume/area info 
  volCyto = params.volCyto
  #volCyto = Constant(assemble(Constant(1.0)*dx()))
  areaCyto = Constant(assemble(Constant(1.0)*ds(lMarker)))
  volSSL = params.volSSL
  volCleft= params.volCleft
  volFrac_CytoSSL= volCyto/volSSL
  volFrac_SSLCyto= 1/volFrac_CytoSSL
  volFrac_SSLCleft= volSSL/volCleft
  volFrac_CleftSSL= 1/volFrac_SSLCleft
  volFrac_CleftCyto= volCleft/volCyto
  volFrac_CytoCleft= 1/volFrac_CleftCyto           

  # System
  
  ## PDE 
  # Time derivative and diffusion of field species
  F = ((cCa_n-cCa_0)*vCa/dt + params.DCa*inner(grad(cCa_n), grad(vCa)))*dx()
  #Ds  = params.Ds
  #RHSis=[]
  #for i in range(cm.nDOF_Fields):
  # #RHSi= -inner(Ds[i]*grad(cns[i]), grad(vs[i]))*dx
  #  RHSi= -inner(Ds[i]*grad(u_n[i]), grad(vs[i]))*dx
  #  RHSis.append(RHSi)
  F+= ((cCaBuff_n-cCaBuff_0)*vCaBuff/dt + params.DBuff*inner(grad(cCaBuff_n), grad(vCaBuff)))*dx()
  F+= ((cCaFluo_n-cCaFluo_0)*vCaFluo/dt + params.DFluo*inner(grad(cCaFluo_n), grad(vCaFluo)))*dx()


  #add in SSL region 
  #what is minimal example Ican get in place for caitlin 
  
  if ssl:
    # Flux to mesh domain from scalar domain s 
    jFlux_SSLCyto = D_SSLCyto*(cCaSSL_n-cCa_n)/dist # note that this is multiplied by ds(lMarker) below
    F += -jFlux_SSLCyto*vCa*ds(lMarker)
  else: 
    # Flux to mesh domain from scalar domain s 
    jFlux = D_CleftCyto*(cCaCleft_n-cCa_n)/dist

    F += -jFlux*vCa*ds(lMarker)


  
  ## Scalar SSL (even if ssl=False, we still define a PDE just ot keep the code more simple) 
  F += (cCaSSL_n-cCaSSL_0)/dt*vCaSSL*dx() 
  
  if ssl:
    # flux from scalar domain s to mesh 
    F += jFlux_SSLCyto*volFrac_CytoSSL*vCaSSL*ds(lMarker)
  
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
    F += volFrac_CytoCleft*jFlux*vCaCleft*ds(lMarker)


  ### Fluxes
  RHS = 0
  if doAssert=="fluxTest_jSurfCyto":
    jSurf = wholeCellJ_to_jSurf(params.jTest,volCyto,areaCyto)
    print assemble(jSurf*ds(lMarker))
    RHS = jSurf*vCa*ds(lMarker)
  if doAssert=="fluxTest_jCleft":
  #  totVol = volCyto + volCleft + volSSL # total of open comparmtnets
  #  areaCleft= 1.  # use 1 for scalar compartment ?
  #  #jSurf = wholeCellJ_to_jSurf(params.jTest,volCleft,areaCleft)
  #  jSurf = params.jTest # since applied to a sclaar? 
  #  #jSurf = Expression("conditional( And( lt(u, 0.5), gt(u, -0.5) ), 1.0, 0.0 )",u=1)
    jCleft = wholeCellJ_to_jCleft(params.jTest,params)
    RHS = jCleft*vCaCleft*dx()
  F += -RHS 


  # adding in fluxes from Torres paper  
  RHS = 0
  if reactions !=None:
    if reactions =="torres":
      import torres
      reactions = torres.Torres()
      reactions.Init(params)
      print "Need to renormalize to agree w whole cell"
      print "Put into module?"
      #RHS = reactions.iryr*(volCleft/areaCleft)*vCaCleft*dx()
      rescaleFactor = 1.
      RHS = reactions.iryr*(rescaleFactor)*vCaCleft*dx()

    if reactions == "simple":
      import simple
      reactions = simple.Simple()
      reactions.Init(params)
      rescaleFactor = 1.
      RHS = reactions.iryr*(rescaleFactor)*vCaCleft*dx()
      rescaleFactor = 1.
      RHS+= reactions.jSERCA*(rescaleFactor)*vCa*dx()
       
    # apply
    F += - RHS

  ## Reactions
  if buffers !=None:
    RHS = 0 
    RCaBuff = -params.alphabuff*(params.Btot - cCaBuff_n)*cCa_n + params.betabuff*cCaBuff_n 
    RHS+=  RCaBuff*vCa*dx()
    RHS+= -RCaBuff*vCaBuff*dx()

    RCaFluo = -params.alphafluo*(params.Ftot - cCaFluo_n)*cCa_n + params.betafluo*cCaFluo_n 
    RHS+=  RCaFluo*vCa*dx()
    RHS+= -RCaFluo*vCaFluo*dx()
    
    F += -RHS
    

    


  
  # Init conditions
  #U_0.interpolate(Constant((params.cCaInit,params.cCaSSLInit, params.cCaCleftInit)))
  print params.cInits
  print "WARNING: concentrations aren't correct for scalar domainsw"
  U_0.interpolate(Constant(params.cInits))
  
  def conservation(t=0):
      ## scalars     
      acCaCleft_n = assemble(cCaCleft_n*volFrac_CleftCyto*dx())
      #ccCaCleft_n = assemble(cCaCleft_n*volFrac_CleftCyto/volCleft*dx())
      # below is shorthand, but i left comment above so its clear where it originated 
      ccCaCleft_n = assemble(cCaCleft_n/volCyto*dx())

      acCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto*dx())
      #ccCaSSL_n = assemble(cCaSSL_n*volFrac_SSLCyto/volSSL*dx())
      ccCaSSL_n = assemble(cCaSSL_n/volCyto*dx())

      ## cyto
      acCa_n = assemble(cCa_n*dx())
      ccCa_n = assemble(cCa_n/volCyto*dx())

      acCaBuff_n = assemble(cCaBuff_n*dx())
      ccCaBuff_n = assemble(cCaBuff_n/volCyto*dx())

      acCaFluo_n = assemble(cCaFluo_n*dx())
      ccCaFluo_n = assemble(cCaFluo_n/volCyto*dx())


      caconserved= acCaSSL_n+acCaCleft_n+acCa_n
      conserved=caconserved
      conserved+= acCaBuff_n+acCaFluo_n 
    
      if MPI.rank(mpi_comm_world())==0:
        print "############### ", t
        print "cCa_n " , acCa_n, "conc ", ccCa_n
        print "cCaSSL_n " , acCaSSL_n, "conc ", ccCaSSL_n
        print "cCaCleft_n " , acCaCleft_n, "conc ", ccCaCleft_n
        print "caconserved " , caconserved
        print "cCaBuff_n " , acCaBuff_n, "conc ", ccCaBuff_n
        print "cCaFluo_n " , acCaFluo_n, "conc ", ccCaFluo_n
        print "CONSERVATION:", conserved

      return conserved

  
  
      
  U_n.vector()[:] = U_0.vector()[:]
  conserved_t0 = conservation()

  ## File IO
  t=0.
  ctr=0
  pvdFile = File(pvdName,"compressed")
  pvdFile << (U_n.sub(idxCa),t)     
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
      conserved_ti = conservation(t=t)

      ## store 
      pvdFile << (U_n.sub(idxCa),t)
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
      ## update
      if reactions!=None:
        reactions.Update(t) 
      if 0: # mode=="sachse":
        1
        #iryr.t = t 
        
        #caTT= assemble(u_n[idxCa]*ds(rMarker,domain=mesh))
        #vol = assemble(Constant(1.)*ds(rMarker,domain=mesh))
        #caTT = caTT/vol 
        #print "TT %f" % caTT
        #iryrDefect.ca = caTT

      if mode=="ode":
        # TODO run ODE model 
        # aggregate state variables
        # need to assemble concs over entire vol 
        #def UpdateBeardVals(nC,bM,mgatpi_uM,fadpi_uM ,fmgi_uM ):
        jCa.r  = 0. # TODO update with ode 

  # close 
  hdf.close()
  
  
  # Check that the gods of transport are happy 
  if doAssert == "conservation":
    message = "%f | %f " % (conserved_ti , conserved_t0)
    assert( abs(conserved_ti - conserved_t0) < 1e-3 ), message
    #assert( abs(acCaSSL_n+acCaCleft_n+acCa_n - 60000) < 1e-5)
    print "PASS"

  # returning for validaton
  concsFinal = np.zeros(cm.nDOF) 
  for i in range(cm.nDOF):
    concsFinal[i] = assemble(U_n[i]/volCyto*dx())
  print concsFinal 
  return concsFinal 


# This case examines the effect of delaying diffusion through the SSL domain 
def mytest():
  params.T = 1000
  params.dt = 5 


  tag = "2D_SSL_torres"
  tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag,reactions="torres")

  tag = "2D_SSL_torres_slowSSL"
  params.D_SSLCyto = Constant(0.3) # Diffusion rate between SSL/Cyto compartments (if SSL exists) [um^2/ms]
  params.D_CleftSSL = Constant(0.3 ) #  Diffusion rate between Cleft/SSL compartments (if SSL exists)
  params.D_CleftCyto= Constant(0.3)
  tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag,reactions="torres")

def mytest2():
  params.T = 1000
  params.dt = 5 


  tag = "2D_noSSL_simple" 
  tsolve(debug=False,params=params,pvdName = "noSSL.pvd", hdfName=tag+".h5",\
    mode=tag,reactions="simple",buffers=True)
  tag = "2D_SSL_simple" 
  tsolve(debug=False,params=params,pvdName = "SSL.pvd", hdfName=tag+".h5",\
    mode=tag,reactions="simple",buffers=True)


  

# In[19]:
def whosalive():
  print "I'm alive!"

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
    # calls 'tsolve' with the next argument following the argument '-validation'
    if(arg=="-debug"):
      #arg1=sys.argv[i+1] 
      tsolve(debug=True)     
    elif(arg=="-test0"):
      tsolve(debug=False)    
      quit()
    elif(arg=="-test1"):
      tag = "ode"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2D_torres"):
      tag = "2D_SSL_torres"
      params.T = 1000
      params.dt = 5 
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag,reactions="torres",buffers=True)
      quit()
    elif(arg=="-test2D"):
      tag = "2D_SSL"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2Dno"):
      tag = "2D_noSSL"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-testsatin"):
      tag = "satin"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test2"):
      tag = "sachse2TT"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-test4"):
      tag = "sachse4TT"
      tsolve(debug=False,params=params,hdfName=tag+".h5",mode=tag)
      quit()
    elif(arg=="-mytest"):
      mytest()
      quit() 
    elif(arg=="-mytest2"):
      mytest2()
      quit() 
    elif(arg=="-validation"):
      validation()
      quit()






  raise RuntimeError("Arguments not understood")



