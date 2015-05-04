from dolfin import * 
import numpy as np

from sarcomereBase import SarcomereBase

## var
ttRad = 0.25 # [um]
eps = 0.05



## Boundary defs
class TopTT(SubDomain):
  def inside(self,x,on_boundary):
    # Define TT loc 
    centroid1 = np.array([self.mmin[0],self.mmax[1]])
    centroid2 = np.array([self.mmax[0],self.mmax[1]])

    # check if point is nearby 
    d1 = np.linalg.norm(centroid1-x[0:2])
    d2 = np.linalg.norm(centroid2-x[0:2])
    isTT1 = (d1 < (ttRad+eps))
    isTT2 = (d2 < (ttRad+eps))
    isTT = isTT1 or isTT2

    #print x,centroid,d,isTT
    #print x[0], edge, on_boundary
    return on_boundary and isTT

class BottomTT(SubDomain):
  def inside(self,x,on_boundary):
    # Define TT loc 
    centroid1 = np.array([self.mmin[0],self.mmin[1]])
    centroid2 = np.array([self.mmax[0],self.mmin[1]])

    # check if point is nearby 
    d1 = np.linalg.norm(centroid1-x[0:2])
    d2 = np.linalg.norm(centroid2-x[0:2])
    isTT1 = (d1 < (ttRad+eps))
    isTT2 = (d2 < (ttRad+eps))
    isTT = isTT1 or isTT2

    #print x,centroid,d,isTT
    #print x[0], edge, on_boundary
    return on_boundary and isTT

# This is where SSL is 
class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[2]- self.mmax[2]) < DOLFIN_EPS)
    #print x[0], edge, on_boundary
    return on_boundary and edge


class sarcomere4TT(SarcomereBase):
  def __init__(self,params="",mode="",geom="2D"):
    SarcomereBase.__init__(self)
    self.mode = mode
    if geom=="2D":
      self.fileName = "siam/sarcomere4TT_2D.xml"
      self.dim = 2
    else:
      self.fileName = "siam/sarcomere4TT.xml"
      self.dim = 3
      raise RuntimeError("Not supported") 
    self.params = params 
    self.distributions()

    # Geometric considerations 
    self.TTHeight = 6. # TT height [um]
    self.sarcWidth = 1. # sarcomere witdth [um]
    #print "WARNING: Need to change"


  def GetMesh(self):
    self.mesh = Mesh(self.fileName)
    return self.mesh 

  def CalcSSLVol(self):
    self.params.volSSL = self.TTHeight*self.sarcWidth*self.params.SSLWidth



  def distributions(self):
    # TT radius 
    xl= 0.25
    sc = 0.01
    #lhs = 1/(1+np.exp((xs-xl)/sc))    
    xr = 2-0.25
    #rhs = 1-1/(1+np.exp((xs-xr)/sc))    
    # location of zline 
    self.zLine = Expression("1/(1+exp((x[0]-xl)/sc)) + 1-1/(1+exp((x[0]-xr)/sc))", \
          xl=xl,xr=xr,sc=sc)
    # location of cytosol
    self.cytosol = Expression("(1-1/(1+exp((x[0]-xl)/sc)))*(1/(1+exp((x[0]-xr)/sc)))", \
          xl=xl,xr=xr,sc=sc)

  def Boundaries(self,subdomains,ttMode="allfunctional"):
    mesh = self.mesh

    boundary = BottomTT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    lMarker = self.params.lMarker
    boundary.mark(subdomains,lMarker)
    #aboundary = DirichletBC(FunctionSpace(mesh,"CG",1),Constant(1.),boundary)
  
    boundary = TopTT()
    # validation aboundary = DirichletBC(FunctionSpace(mesh,"CG",1),Constant(1.),boundary)
  
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    rMarker = self.params.rMarker

    if ttMode=="allfunctional":
      boundary.mark(subdomains,lMarker)
    else:
      boundary.mark(subdomains,rMarker)
      raise RuntimeError("Not supported yet") 

  
  
    boundary = OuterSarcolemma()
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = self.params.slMarker
    if self.dim==3:
      boundary.mark(subdomains,slMarker)
      raise RuntimeError("Not supported yet") 


    # for debugging 
    #f = Function(FunctionSpace(mesh,"CG",1))
    #aboundary.apply(f.vector())
    #File("boundary.pvd")<<f
    #quit()
  
    return lMarker,rMarker,slMarker

  
    
    
  
  
