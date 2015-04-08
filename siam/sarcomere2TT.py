from dolfin import * 
import numpy as np

from sarcomereBase import SarcomereBase

## var
ttRad = 0.25 # [um]
eps = 0.05


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


# This is where SSL is 
class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[2]- self.mmax[2]) < DOLFIN_EPS)
    #print x[0], edge, on_boundary
    return on_boundary and edge

##  ADDED BY CES
#class sarcomere2TT(SarcomereBase):
#  def __init__(self,params="",mode="",geom="2D"):
#    SarcomereBase.__init__(self)
#    self.mode = mode
#    if geom=="2D":
#      self.fileName = "siam/sarcomere2TT_2D.xml"
#      self.dim = 2
#    else:
#      self.fileName = "./siam/sarcomere2TT.xml"
#      self.dim = 3
#    self.params = params
#    self.distributions() 

class sarcomere2TT(SarcomereBase):
  def __init__(self,params="",mode=""):
    SarcomereBase.__init__(self)
    self.mode = mode
    self.fileName = "./siam/sarcomere2TT.xml"
    self.params = params

  def GetMesh(self):
    self.mesh = Mesh(self.fileName)
    return self.mesh 

  def Boundaries(self,subdomains):
    mesh = self.mesh

    boundary = LeftTT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    lMarker = 2
    boundary.mark(subdomains,lMarker)
  
    boundary = RightTT()
  
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    rMarker = 3
    boundary.mark(subdomains,rMarker)
  
    boundary = OuterSarcolemma()
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = 4
    boundary.mark(subdomains,slMarker)


    return lMarker,rMarker,slMarker

    
  
  

  
  
  
