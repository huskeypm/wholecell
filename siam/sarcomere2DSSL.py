from dolfin import * 
import numpy as np


## var
eps = 0.05



## Boundary defs
class Outersarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    isTT = np.abs(self.mmax[0]-x[0]) < eps

    #print x[0], edge, on_boundary
    return on_boundary and isTT



class sarcomere2DSSL():
  def __init__(self,mode="wSSL"):
    self.mode = mode

  def GetMesh(self):
    self.mesh = UnitSquareMesh(16,16)
    return self.mesh 

