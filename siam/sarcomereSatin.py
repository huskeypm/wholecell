from dolfin import * 
import numpy as np


## var
eps = 0.05



class sarcomereSatin():
  def __init__(self,mode=""):
    self.mode = mode
    self.fileName = "./siam/sarco.xml" 

  def GetMesh(self):
    self.mesh = Mesh(self.fileName)
    return self.mesh 

