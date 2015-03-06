from dolfin import *

class SarcomereBase(object):
  def __init__(self):
    self.fileName = "NONE"
    self.idxCa = 0 # PDE 
    self.idxBuff = 1
    self.idxFluo = 2
    self.idxCaCleft = 3
    self.idxCaSSL = 4  # not always used

  def GetMesh(self):
    return 1 
 
  def GetReactions(self):
    return 1  

  def CalcGeomAttributes(self):
    V = FunctionSpace(self.mesh,"CG",1)
    self.area = assemble(Constant(1.)*ds(domain=self.mesh))
    self.volume  = assemble(Constant(1.)*dx(domain=self.mesh))

  def GetIndices(self):
    return [self.idxCa,self.idxBuff,self.idxFluo,\
            self.idxCaCleft,self.idxCaSSL]

