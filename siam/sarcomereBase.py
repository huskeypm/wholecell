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

    # Follow sharelatex notes (Geometric considerations) for this 
    # Should really go in the sachse parameters section 
    self.nSarcomeres = 1. 

  def Init(self):
    return 1 

  def GetMesh(self):
    return 1 
 
  def GetReactions(self):
    return 1  

  #set by individual sarcomere objecst 
  def CalcSSLVol(self):
    raise RuntimeError("Need to overload this") 
    self.params.volSSL = False
    

  def CalcGeomAttributes(self):
    #V = FunctionSpace(self.mesh,"CG",1)
    # pass in measures, markers 
    ds = self.ds
    dx = self.dx
    lMarker = self.params.lMarker
    rMarker = self.params.rMarker
    slMarker = self.params.slMarker
    
    #self.area = assemble(Constant(1.)*ds(domain=self.mesh))
    self.area = assemble(
      Constant(1.)*ds(lMarker)+\
#      Constant(0.)*ds(rMarker)+\
      Constant(1.)*ds(slMarker)\
    )
    #self.volume  = assemble(Constant(1.)*dx(domain=self.mesh))
    self.volume  = assemble(Constant(1.)*dx)

    print "(Active) Area ", self.area, " Vol " , self.volume

    # compute total volume of CRUs
    self.nCRUs = self.params.nCRUs_per_sarco * self.nSarcomeres
    self.params.volCleft = self.nCRUs * self.params.volPerCRU 

    self.CalcSSLVol()

  def GetIndices(self):
    return [self.idxCa,self.idxBuff,self.idxFluo,\
            self.idxCaCleft,self.idxCaSSL]

