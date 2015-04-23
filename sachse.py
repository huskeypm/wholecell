#
# Time dependent solver, multiple species 
# 
# TODO 
# rescale simple.iRyR correctly 
# add real SERCA model to simple  
#- defining SSL region within cyto #
# - add diffusivity/buffer as D' = D B/(1+KD)
# DONE - get concentrations of buffs, etc correct
# DONE Merge in sachse_temp, sachse_REPLACEMENT?
# DONE  meaningful iRyr
#
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



# [CaB] = [Btot]/(KD/[Ca]+1)
def buffered(totB,freeCa,kFwd,kback):
    KD = kback/kFwd
    return totB/(KD/freeCa + 1) 

class Boundary(SubDomain):
  def inside(self,x,on_boundary):
    left= x[0] <DOLFIN_EPS
    #print left 
    return left and on_boundary

class Params(object):
  def __init__(self): 
    #useMPI = True
    #paraview = True
    self.verbose=False
  
    self.T = 200.  # Total simulation time [ms]
    self.dt = 1. # time-step size [ms] 
  
    self.DCompart = 1 # [um/ms]
    self.CaitlinsNewVariabvle = 1 # [INITS!!]
    self.D_SSLCyto = self.DCompart # Diffusion rate between SSL/Cyto compartments (if SSL exists) [um/ms]
    self.D_CleftSSL = self.DCompart #  Diffusion rate between Cleft/SSL compartments (if SSL exists)
    self.D_CleftCyto= self.DCompart #  Diffusion rate between Cleft/Cyto compartments (if SSL does not exist)
    self.dist = 0.05 # distance between compartments [um] 
  
    # diffusion params 
    self.DCa = 0.39  # Diff const. within cytosol [um^2/ms] verified
    self.DBuff = Constant(0.) # TnC
    self.DFluo = 0.1 # [um^2/ms] verified 
    #Ds = [DCa,DBuff,DFluo]
    self.DCa_SSL = 0.6 * self.DCa # Diff const within SSL region (of cytosol) 
  
    # geom 
    self.volSSL = 1. # Volume of SSL domain [um^3] (overridden by sarcomere files)
    self.volCleft = 0.1 # [um^3] 
  
    # buffer (mostly TnC) 
    self.alphabuff=0.04 # kon verified [1/uMms] 
    self.Btot = 70. # [TnC] [uM] verified 
    self.betabuff = 0.04 # koff verfied [1/ms]
   
    # buffer (imaginary for SSL) 
    self.alphabuff_SSL=0.04 # kon verified [1/uMms] 
    self.Btot_SSL = 70. # [TnC] [uM] verified 
    self.betabuff_SSL = 0.04 # koff verfied [1/ms]
      
  
    # fluo
    self.alphafluo=0.23 # kon fluo  [1/uMms] verified 
    self.Ftot = 50. # [Fluo] [uM] verified  
    self.betafluo = 0.17 # koff fluo [1/ms] verified
  
    # RyR
    # from Torres
    # self.ryrAmp = 9.5 # [pA/pF] #CES
    self.ryrAmp = 4.5 # [pA/pF] #CES
    self.ryrOffset = 5 # [ms]
##    self.ryrTau = -50/np.log(1/2.) # half-max amp at 50 ms  #CES
    self.ryrTau = 274.243089 #CES
    self.ryrKm = 0.2 # uM (made this up)  

    # from Soeller
    self.sercaVmax = 0.2 # uM/ms
    self.sercaKmf = 0.184 # uM
    self.sercaH = 4.  # Hill coeff 
  
    # init concs 
    self.cCaInit =0.1  # Initial Ca concentration in Cleft compartment [uM] 
    self.cBuffInit = buffered(self.Btot,self.cCaInit,self.alphabuff,self.betabuff)  # TnC unverified [uM]
    self.cFluoInit = buffered(self.Ftot,self.cCaInit,self.alphafluo,self.betafluo)  # [uM]
    self.cCaSSLInit = self.cCaInit
    self.cCaCleftInit = self.cCaInit
    self.cInits = np.array([self.cCaInit,self.cBuffInit,self.cFluoInit,self.cCaCleftInit,self.cCaSSLInit])
  
  
    self.jTest = 0.1 # Generic flux applied to cytosolic domain [uM/ms]


# J [uM/ms] --> j[uM*um/ms]
# volDomain - domain in which flux is applied
# areaDomainFlux - area over which flux is applied 
def wholeCellJ_to_jSurf(wholeCellJ,params,volDomain,areaDomainFlux):
  jCyto = wholeCellJ_to_jCyto(wholeCellJ,params)
  return jCyto*volDomain/areaDomainFlux


# Converts whole cell current I [A/F] into volume flux 
# For rabbit (C/V=4.5), 1 [A/F]  --> ~0.02 [uM/ms]
def wholeCellI_to_wholeCellJ(wholeCellI):
  # from flux notes on sharelatex
  nIons = 2 # number of ions per atom (2 for Ca, +/- 1 for NaCa2+ exchange)
  Fc = 9.6e-2 # Faraday [umol/C]
  Cm_o_Vm = 4.5 # Capacitance over cell volume [Farad/Liter]
  s_to_ms = 1e3
  
  jWholeCell = 1/(nIons*Fc) * Cm_o_Vm * wholeCellI * (1/s_to_ms)
  return jWholeCell
  
# assumes all compartments are open
def wholeCellJ_to_jCyto(wholeCellJ, params):
  cytoRescale = (params.volCleft + params.volSSL + params.volCyto)
  cytoRescale/=params.volCyto 
  return wholeCellJ*cytoRescale # [uM/ms]

def wholeCellJ_to_jCleft(wholeCellJ, params):
  cleftRescale = (params.volCleft + params.volSSL + params.volCyto)
  cleftRescale/=params.volCleft
  return wholeCellJ*cleftRescale # [uM/ms]


def tsolve(pvdName=None,\
           mode="2D_SSL", # sachse2TT, sachse4TT, ode, bcs \ 
           params = Params(),\
           hdfName = "out.h5",
           doAssert = None,  
           reactions = "simple", #None,
           existsCleft = True,
           existsSSL = True,
           buffers = True, # buffers = False,
           debug=False):

  #mesh = UnitCubeMesh(8,8,8)
  ssl=True
  if "2D" in mode:
    import sarcomere2DSSL
    cm = sarcomere2DSSL.Sarcomere2DSSL(params=params)
    if "noSSL" in mode:
      ssl=False

  elif "sachse" in mode: 
    if "sachse4TT" in mode:
      import sarcomere4TT
      cm = sarcomere4TT.sarcomere4TT(params=params,geom="2D")

    else: 
      import sarcomere2TT
      cm = sarcomere2TT.sarcomere2TT(params=params)

  elif "satin" in mode:
    import sarcomereSatin
    cm = sarcomereSatin.sarcomereSatin(params=params)



  cm.ssl = ssl
  if ssl==False:
    existsSSL=False

  cm.Init()
  cm.params.nDOF = cm.nDOF # goes in master class 
  mesh = cm.GetMesh()
  dims =mesh.ufl_cell().geometric_dimension()
  #print cm.volume
  cm.CalcGeomAttributes()
  params.volCyto= cm.volume

  ## 
  print 'WARNING: will need to updated these'
  idxCa = 0 # PDE 
  idxBuff = 1
  idxFluo = 2
  idxCaSSL = 1
  idxCaCleft = 2  # not always used


  # rescale the IMMOBILE buffer by volumeFrac cyto (when SSL/cyto are compiled) 
  if buffers==False:
    params.Btot = 0.
    params.cInits[idxBuff] = 0.
    params.Ftot = 0.
    params.cInits[idxFluo] = 0.

  if ssl==False:
    params.Btot = params.Btot*cm.volFracCyto
    if "2D" not in mode:
      raise RuntimeError("Have not written support for this yet") 



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
  
  
  # iterator
  dt = params.dt 
  tstop = params.T
  
  # volume/area info 
  volCyto = params.volCyto
  #volCyto = Constant(assemble(Constant(1.0)*dx()))
  areaCyto = Constant(assemble(Constant(1.0)*ds(lMarker)))
  volSSL = params.volSSL
  volCleft= params.volCleft
  volFrac_CytoSSL= volCyto/np.max([1e-9,volSSL])
  volFrac_SSLCyto= 1/volFrac_CytoSSL
  volFrac_SSLCleft= volSSL/np.max([1e-9,volCleft])
  volFrac_CleftSSL= 1/np.max([1e-9,volFrac_SSLCleft])
  volFrac_CleftCyto= volCleft/volCyto
  volFrac_CytoCleft= 1/np.max([1e-9,volFrac_CleftCyto])
  if existsCleft==False:
    params.D_CleftCyto=0.
    volFrac_CleftCyto= 0.
    volFrac_CytoCleft= 0.
    params.volCleft = 1e-9
  if existsSSL==False:
    volFrac_CytoSSL= 0.
    volFrac_SSLCyto= 0.
    volFrac_SSLCleft= 0.
    volFrac_CleftSSL= 0.
    params.volSSL = 1e-9
    params.D_SSLCyto = 0.
    params.D_CleftSSL = 0.

  if ssl!=True:
    params.volSSL = 1e-9   # needed for validaton stuff


  # intercompartment transport 
  D_SSLCyto = params.D_SSLCyto
  D_CleftSSL = params.D_CleftSSL
  D_CleftCyto = params.D_CleftCyto
  dist = params.dist

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
    jSurf = wholeCellJ_to_jSurf(params.jTest,params,volCyto,areaCyto)
    RHS = jSurf*vCa*ds(lMarker)
  elif doAssert=="fluxTest_jCleft":
    jCleft = wholeCellJ_to_jCleft(params.jTest,params)
    RHS = jCleft*vCaCleft*dx()
  elif doAssert=="fluxTest_jVolCyto":
    jVol = wholeCellJ_to_jCyto(params.jTest,params)
    RHS = jVol*vCa*dx()

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

    elif reactions == "simple":
      import simple
      reactions = simple.Simple()
      #reactions = simple.Simple(noSERCA=True)
      reactions.Init(params)
      itoJ = wholeCellI_to_wholeCellJ(1.) # rescale 1 [A/F] --> [uM/ms]
      itoj = wholeCellJ_to_jCleft(itoJ,params) 
      RHS += reactions.iryr*itoj*vCaCleft*dx()
      # SERCA
      Jtoj = wholeCellJ_to_jCyto(1.0,params) # 1 [uM/ms] --> [uM um /ms]
      RHS+= reactions.jSERCA*(Jtoj)*vCa*dx()
      #RHS+= -0.1*(Jtoj)*vCa*dx()


    elif "ryrOnly" in reactions:
      import simple
      if reactions=="ryrOnlySwitch":
        reactions = simple.Simple(ryrOnlySwitch=True)
      else:
        reactions = simple.Simple()
      reactions.Init(params)
      itoJ = wholeCellI_to_wholeCellJ(1.) # rescale 1 [A/F] --> [uM/ms]
      itoj = wholeCellJ_to_jCleft(itoJ,params) 
      RHS += reactions.iryr*itoj*vCaCleft*dx()

    elif reactions=="caitlinSERCA": #CES
      import simple 
      reactions = simple.Simple(caitlinSERCA=True) #CES
      reactions.Init(params)
      # SERCA
      Jtoj = wholeCellJ_to_jCyto(1.0,params) # 1 [uM/ms] --> [uM um /ms]
      RHS+= reactions.jSERCA*(Jtoj)*vCa*dx()

      
    # apply
    F += - RHS


  ## Reactions
  if buffers !=False: 
    RHS = 0 
    RCaBuff = -params.alphabuff*(params.Btot - cCaBuff_n)*cCa_n + params.betabuff*cCaBuff_n 
    RHS+=  RCaBuff*vCa*dx()
    RHS+= -RCaBuff*vCaBuff*dx()

    RCaFluo = -params.alphafluo*(params.Ftot - cCaFluo_n)*cCa_n + params.betafluo*cCaFluo_n 
    RHS+=  RCaFluo*vCa*dx()
    RHS+= -RCaFluo*vCaFluo*dx()
    
    F += -RHS
  else:
    params.Btot = 0.
    params.Ftot = 0.

    


  
  ## Apply Init conditions
  #U_0.interpolate(Constant((params.cCaInit,params.cCaSSLInit, params.cCaCleftInit)))
  U_0.interpolate(Constant(params.cInits))


  # here we revise the init conds
  if ssl==False:
    f = Function(V)
    #f.interpolate(Constant(params.cInits[idxBuff]))
    pFunc = cm.pCyto
    idxBuff = cm.idxBuff  
    idxs = [cm.idxBuff,cm.idxFluo]
    for i,idxi in enumerate(idxs):
      pFunc.m = params.cInits[ idxi    ] 
      f.interpolate(pFunc)
      #f.interpolate(Expression("x[0]"))
      assign(U_0.sub(idxi)   ,f)
      pvdFile = File("init.pvd","compressed")
      pvdFile << U_0.sub(idxi)        
 
  # left hand side of box has Ca2+, Ca2+/Buff present 
  #cInits = params.cInits
  #class InitialConditions(Expression):
  #  def eval(self, values, x):
  #    values[0] = cInits[0]
  #    values[1] = cInits[1]
  #    values[2] = cInits[2]
  #    values[3] = Expression("1") 
  #    values[4] = cInits[4]
  #  def value_shape(self):
  #    return (5,)
  #ic = InitialConditions()
  #U_0.interpolate(ic)
  
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
  if pvdName!=None:
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
      if pvdName != None:
        pvdFile << (U_n.sub(idxCa),t)
      # store hdf 

      writeHdf=True
      if writeHdf:
        uCa = project(cCa_n,V)
        uCa.vector()[:]/=volCyto
        #print assemble(uCa*dx())
        hdf.write(uCa,"uCa",ctr)
      
        # see conservation routine above for why we multiply by volCyto 
        uCaSSL = project(cCaSSL_n,R)      
        uCaSSL.vector()[:]/=volCyto
        #print assemble(uCaSSL*dx())
        hdf.write(uCaSSL,"uCaSSL",ctr)

        uCaCleft = project(cCaCleft_n,R)         
        uCaCleft.vector()[:]/=volCyto
        #print assemble(uCaCleft*dx())
        hdf.write(uCaCleft,"uCaCleft",ctr)


  
      ## update 
      ctr+=1
      U_0.assign(U_n)
      t += float(dt)
      ## update
      if reactions!=None:
        caiAvg = assemble(cCa_n*dx())/volCyto
        reactions.Update(t,caiAvg)  
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

def test4TT():
  params = Params()
  params.T = 500
  params.dt = 1.0
  params.T = 150  
  params.dt = 5.0  
  params.ryrOffset = 0. # maybe why nonzero fluxes all over seen?

  params.ryrAmp = 30 
  params.sercaVmax = 0.5

  reactions = "simple"

  ## REAPID 
  suffix = "_prod"
  h5Tag = "_rapid"+suffix
  rapid = 1e2
  params.D_SSLCyto = rapid # Can't go much faster than this   
  params.D_CleftSSL= rapid
  params.D_CleftCyto = rapid

  tag = "sachse4TT_SSL_%s" % reactions
  tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
    hdfName=tag+h5Tag+".h5",\
    mode=tag,reactions=reactions,buffers=True)

  #tag = "2D_noSSL_%s" % reactions
  #tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
  #  hdfName=tag+h5Tag+".h5",\
  #  mode=tag,reactions=reactions,buffers=True)

  ## REDUCED 
  #h5Tag = "_red"+suffix

def simpleCompare():
  params = Params()
  params.T = 1000
  params.dt = 5 


  params.ryrAmp = 30 
  params.sercaVmax = 0.5

  reactions = "ryrOnly" 
  reactions = "simple"

  ## REAPID 
  suffix = "_prod"
  h5Tag = "_rapid"+suffix
  rapid = 1e2
  params.D_SSLCyto = rapid # Can't go much faster than this   
  params.D_CleftSSL= rapid
  params.D_CleftCyto = rapid

  tag = "2D_SSL_%s" % reactions
  tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
    hdfName=tag+h5Tag+".h5",\
    mode=tag,reactions=reactions,buffers=True)

  tag = "2D_noSSL_%s" % reactions
  tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
    hdfName=tag+h5Tag+".h5",\
    mode=tag,reactions=reactions,buffers=True)

  ## REDUCED 
  h5Tag = "_red"+suffix
  red   = 1e-1
  params.D_SSLCyto = red   # Can't go much faster than this   
  params.D_CleftSSL= red  
  params.D_CleftCyto = red  

  tag = "2D_SSL_%s" % reactions
  tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
    hdfName=tag+h5Tag+".h5",\
    mode=tag,reactions=reactions,buffers=True)

  tag = "2D_noSSL_%s" % reactions
  tsolve(debug=False,params=params,pvdName=tag+h5Tag+".pvd",
    hdfName=tag+h5Tag+".h5",\
    mode=tag,reactions=reactions,buffers=True)



  

# Verifying that we get the correct free Calcium if we merege 
# the SSL/Cyto into a single domain 
# TODO understand why there is still some numerical error. Did this at 5 am 
# Total conservation is not exact; might have to do with interpoaltion of the diract function
# used for defining the SSL region 
def validationMergingSSLCyto():
  params = Params()
  params.T = 1000   
  params.dt = 500

  idxCa = 0
  idxFluo= 2 # s.b. generalized 
  idxCaCleft = 4 # s.b. generalized 

  reactions = "ryrOnlySwitch"
  reactions = "ryrOnly"
  reactions = None     

  params.D_SSLCyto = 1e-1 # Can't go much faster than this   
  params.D_CleftSSL= 1e-1
  params.D_CleftCyto = 0.#1e-1
  params.dist = 0.01

  # kill Fluo
  params.Ftot = 0.
  params.cInits[idxFluo] = 0. 

  buffers = True  
  if buffers==False:
    params.Btot = 0.
    params.Ftot = 0.
 
  # separate SSL/Cyto compartments 
  case = "new"
  if 1: 
    mode = "2D_SSL"
    threeComps = tsolve(mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
      reactions = reactions,buffers=buffers) 
 
  # merged compartments 
  print "###\n###\n####\n"
  mode = "2D_noSSL"
  twoComps = tsolve(pvdName="test.pvd",mode=mode,params=params,hdfName="%s_%s.h5"%(mode,case),
    reactions = reactions,buffers=buffers) 
  msg = "%f != %f " %( threeComps[idxCa] ,twoComps[idxCa])
  # NOTE: more generous with this error 
  assert(abs(threeComps[idxCa] - twoComps[idxCa]) < 1e-2), msg
  print "PASSED (but worth a couple check)" + msg 

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
      raise RuntimeError("Not suported") 
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
    elif(arg=="-testSimple"): 
      params.T = 1000
      params.dt = 5 
      tag = "2D_SSL_simple" 
      tsolve(debug=False,params=params,pvdName = "SSL.pvd", hdfName=tag+".h5",\
        mode=tag,reactions="simple",buffers=True)
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
    elif(arg=="-test4TT"):
      test4TT()
      quit()
    elif(arg=="-mytest"):
      mytest()
      quit() 
    elif(arg=="-simpleCompare"):
      simpleCompare()
      quit() 
    elif(arg=="-validation"):
      raise RuntimeError("RUN python validations.py INSTEAD\nCODE NOT VALIDATED!!!") 
      quit()






  raise RuntimeError("Arguments not understood")



