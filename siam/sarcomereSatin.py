from dolfin import * 
import numpy as np

from sarcomereBase import SarcomereBase


## var
eps = 0.05


#class TopTT(SubDomain):
class TT(SubDomain):
  def inside(self,x,on_boundary):
    # Define TT loc 
    if (
      (x[0]>(eps+self.mmin[0]) and x[0]<(-eps+self.mmax[0]) and \
       x[1]>(eps+self.mmin[1]) and x[1]<(-eps+self.mmax[1]))):
      isIn = True
    else:
      isIn = False

    #isAbove = (x[1]>(-eps+self.mmax[1]/2.))
    isAbove = True # want all to be true 

    #print x[0], edge, on_boundary
    result = on_boundary and isIn and isAbove
    #print result 
    return result 

#class BottomTT(SubDomain):
#  def inside(self,x,on_boundary):
#    # Define TT loc 
#    if (
#      (x[0]>(eps+self.mmin[0]) and x[0]<(-eps+self.mmax[0]) and \
#       x[1]>(eps+self.mmin[1]) and x[1]<(-eps+self.mmax[1]))):
#      isIn = True
#    else:
#      isIn = False
#
#    isBelow = (x[1]<(-eps+self.mmax[1]/2.))
#
#    #print x[0], edge, on_boundary
#    result = on_boundary and isIn and isBelow
#    #print result 
#    return result 

# This is where SSL is 
class OuterSarcolemma(SubDomain):
  def inside(self,x,on_boundary):
    edge = (np.abs(x[0]- self.mmax[0]) < DOLFIN_EPS)
    #print x[0], edge, on_boundary
    return on_boundary and edge

class sarcomereSatin(SarcomereBase):
  def __init__(self,params="",mode=""):
    SarcomereBase.__init__(self)
    self.mode = mode
    self.fileName = "./siam/sarco.xml" 
    self.fileName = "/home/AD/pmke226/labscripts/mesh/test.xml"
    self.params = params 
   
    # geometric considerations 
    # WARNING: need to make sure these preserve the aspect ratio 
    # of the original image 
    self.nLongitudinalSarc = 2.25   # sarcomeres in longitudinal direction (by eyeball)  
    self.nAxialSarc = 3   # sarcomeres in axial direction (by eyeball)  
    self.lenSarc = 2.0 # sarcomere leng [um]
    self.widSarc = 1.0659 # sarcomere width [um] (to fit aspect ratio)

  def GetMesh(self):
    return self.mesh 



  def RescaleMesh(self):
    # computee dim 
    self.width = self.widSarc*self.nAxialSarc
    self.length= self.lenSarc*self.nLongitudinalSarc

    # normalize 
    c = self.mesh.coordinates()
    #print np.shape(c)
    delta = np.max(c,axis=0) - np.min(c,axis=0)
    print "Aspect ratio: %f" % (delta[0]/delta[1])
    #print np.max(c,axis=0)
    c = c-np.min(c,axis=0)
    c = c/np.max(c,axis=0)
    #print np.max(c,axis=0)

    # rescale 
    c = c*np.array([self.width,self.length]) 
    delta = np.max(c,axis=0) - np.min(c,axis=0)
    print "NewAspect ratio: %f" % (delta[0]/delta[1])
    #print np.max(c,axis=0)
    
    self.mesh.coordinates()[:] = c 
    quit()

  def Init(self):
    self.mesh = Mesh(self.fileName)
    self.RescaleMesh()      

## NOTE: for now we'll just label everything inside the geometry as 
## TTubule 
  def Boundaries(self,subdomains):
    mesh = self.mesh

    boundary = TT()
    boundary.mmin = np.min(mesh.coordinates(),axis=0)
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    lMarker = 2
    boundary.mark(subdomains,lMarker)

    rMarker = False 

    boundary = OuterSarcolemma()
    boundary.mmax = np.max(mesh.coordinates(),axis=0)
    slMarker = 4
    boundary.mark(subdomains,slMarker)

    # To double check BCs
    #export PYTHONPATH=/home/AD/pmke226/labscripts/dolfin/:$PYTHONPATH
    #import view
    #view.PrintSubdomains(mesh,subdomains)
    #quit()

    return lMarker,rMarker,slMarker
 
    
  
  


