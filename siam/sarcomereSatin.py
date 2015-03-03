from dolfin import * 
import numpy as np

from sarcomereBase import *


## var
eps = 0.05


class TopTT(SubDomain):
  def inside(self,x,on_boundary):
    # Define TT loc 
    if (
      (x[0]>(eps+self.mmin[0]) and x[0]<(-eps+self.mmax[0]) and \
       x[1]>(eps+self.mmin[1]) and x[1]<(-eps+self.mmax[1]))):
      isIn = True
    else:
      isIn = False

    isAbove = (x[1]>(-eps+self.mmax[1]/2.))

    #print x[0], edge, on_boundary
    result = on_boundary and isIn and isAbove
    #print result 
    return result 

class BottomTT(SubDomain):
  def inside(self,x,on_boundary):
    # Define TT loc 
    if (
      (x[0]>(eps+self.mmin[0]) and x[0]<(-eps+self.mmax[0]) and \
       x[1]>(eps+self.mmin[1]) and x[1]<(-eps+self.mmax[1]))):
      isIn = True
    else:
      isIn = False

    isBelow = (x[1]<(-eps+self.mmax[1]/2.))

    #print x[0], edge, on_boundary
    result = on_boundary and isIn and isBelow
    #print result 
    return result 

# This is where SSL is 
class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmax[0]) < DOLFIN_EPS)
    #print x[0], edge, on_boundary
    return on_boundary and edge

class sarcomereSatin(sarcomereBase):
  def __init__(self,params="",mode=""):
    self.mode = mode
    self.fileName = "./siam/sarco.xml" 
    self.nDOF_Fields= 3
    self.nDOF_Scalars= 1 # cleft 
    self.nDOF = self.nDOF_Fields + self.nDOF_Scalars
    self.params = params 

  def GetMesh(self):
    self.mesh = Mesh(self.fileName)
    return self.mesh 

  def Boundaries(self,subdomains):
    mesh = self.mesh

    boundary = BottomTT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    lMarker = 2
    boundary.mark(subdomains,lMarker)

    boundary = TopTT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    rMarker = 3
    boundary.mark(subdomains,rMarker)

    boundary = OuterSarcolemma()
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
    
    
  
  


