from dolfin import * 
import numpy as np


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

def GetMesh():
  fileName = "sarcomere4TT.xml"
  mesh = Mesh(fileName)
  return mesh 


class sarcomere4TT():
  def __init__(self,mode=""):
    self.mode = mode
    self.fileName = "sarcomere2TT.xml"
    self.distributions()

  def GetMesh(self):
    self.mesh = Mesh(self.fileName)
    return self.mesh 

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


