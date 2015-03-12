from dolfin import *

class SarcomereBase(object):
  def __init__(self):
    self.fileName = "NONE"
    self.idxCa = 0 # PDE 
    self.idxBuff = 1
    self.idxFluo = 2
    self.idxCaCleft = 3
    self.idxCaSSL = 4  # not always used
    self.pSSL= Constant(0.) # defines where SSL region is within cytosol (0 for none, P(x<Val)=1 otherwise)  
    self.pCyto= Constant(1.) # defines where SSL region is within cytosol (0 for none, P(x<Val)=1 otherwise)  
    self.nDOF_Fields= 3
    self.nDOF_Scalars= 2  # cleft, SSL 
    self.nDOF = self.nDOF_Fields + self.nDOF_Scalars


  def Init(self):
    return 1 

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

