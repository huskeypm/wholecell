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



class sarcomere2DwoSSL(SarcomereBase):             
  def __init__(self,params=""):
    SarcomereBase.__init__(self)

    self.nDOF_Fields= 3
    self.nDOF_Scalars= 1  # cleft
    self.nDOF = self.nDOF_Fields + self.nDOF_Scalars
    self.params = params 

  def GetMesh(self):
    self.mesh = UnitSquareMesh(16,16)
    return self.mesh 

  def Boundaries(self,subdomains):
    mesh = self.mesh

    boundary = TT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    lMarker = 2
    boundary.mark(subdomains,lMarker)
  
    rMarker = -1
  
    boundary = OuterSarcolemma()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = 4
    boundary.mark(subdomains,slMarker)

    return lMarker,rMarker,slMarker

  # Need to manually put in DOF for now  
  class InitialConditions(Expression):
    def eval(self, values, x):
      for i in range(self.params.nDOF):
              #print i 
              values[i] = self.params.cInits[i]
    def value_shape(self):
      return (4,)             
    
    
  
  

