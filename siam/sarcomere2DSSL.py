from dolfin import * 
import numpy as np

from sarcomereBase import SarcomereBase



## var
eps = 0.05



## Boundary defs
class TT(SubDomain):
  def inside(self,x,on_boundary):
    isTT = np.abs(self.mmin[0]-x[0]) < eps

    #print x[0], edge, on_boundary
    return on_boundary and isTT

class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    isSL = np.abs(self.mmin[1]-x[1]) < eps

    #print x[0], edge, on_boundary
    return on_boundary and isSL



class Sarcomere2DSSL(SarcomereBase):             
  def __init__(self,params=""):
    SarcomereBase.__init__(self)
    self.params = params 
    self.TT = TT
    self.OuterSarcolemma = OuterSarcolemma
    self.ssl = False

    self.TTHeight = 6. # TT height [um]
    self.sarcWidth = 1. # sarcomere witdth [um]
    self.SSLWidth= 0.25 # define ssl as x=0..xSSL [nm]
    self.params.volSSL = self.TTHeight*self.sarcWidth*self.SSLWidth
    

  def GetMesh(self):
    #V = FunctionSpace(self.mesh,"CG",1)
    #self.area = assemble(Constant(1.)*ds(domain=self.mesh))
    #self.volume  = assemble(Constant(1.)*dx(domain=self.mesh))
    return self.mesh 


  def Init(self): 
    ## Get mesh 
    mesh = UnitSquareMesh(16,16)
    c = mesh.coordinates()[:]

    #c[:,1]*= 6. # TT height [um]

    # if in no ssl mode, rescale cytosol to subsume SSL compartment 
    if self.ssl:
      c[:,0]*= self.TTHeight #  [um]
      c[:,1]*= self.sarcWidth #  [um]
    else:
      c[:,0]*= self.TTHeight #  [um]
      c[:,1]*= self.SSLWidth + self.sarcWidth #  [um]

    mesh.coordinates()[:] = c 


    self.mesh = mesh 
    self.mmin = np.min(mesh.coordinates(),axis=0)
    self.mmax = np.max(mesh.coordinates(),axis=0)


    ## Define SSL region 
    if self.ssl==False:
      return 1 

    # below-listing things are broken for now 
    # rigorously unit-test before checking in 
    return 1

    # params 
    xSSL= 0.25 # define ssl as x=0..xSSL [nm]


    ## resize cyto domain to include SSL 
    print "NEEDS to be made kosher with SSL compart volume"
    xRange = self.mmax[0] - self.mmin[0]
    xNew = xSSL+xRange
    self.mesh.coordinates()[:]*=xNew/xRange

    # replace Constant diffusion constant with expression 
    sc = 0.01

    # using a 'rescaled' dirac function to interpolate between 
    # lower diffusion region and higher diffusion region 
    # Deff = dirac*(DCa-DSSL) + DSSL
    params = self.params
    self.DCaBase = params.DCa
    params.DCa = Expression("1/(1+exp((x[0]-xSSL)/sc))*(DCa-DCa_SSL) + DCa_SSL ",\
          xSSL=xSSL,DCa = self.DCaBase, DCa_SSL=params.DCa_SSL,sc=sc)


    # location of SSL region 
    self.pSSL = Expression("1/(1+exp((x[0]-xSSL)/sc)) ",\
          xSSL=xSSL,sc=sc)

  def Boundaries(self,subdomains):
    mesh = self.mesh

    boundary = self.TT()
    boundary.mmin = self.mmin
    boundary.mmax = self.mmax
    lMarker = 2
    boundary.mark(subdomains,lMarker)
  
    rMarker = -1
  
    boundary = self.OuterSarcolemma()
    boundary.mmin = self.mmin
    boundary.mmax = self.mmax
    slMarker = 4
    boundary.mark(subdomains,slMarker)

    return lMarker,rMarker,slMarker

    
  
  

