
from dolfin import *


class Torres(): 
  def __init__(self):
    1
    
  def Init(self,params):
    self.params = params
    self.iryr = Expression("a*exp(-(t-to)/tau)",a=params.ryrAmp,\
                                    to=params.ryrOffset,\
                                    tau=params.ryrTau,\
                                    t=0)

    print "Borrow NCX, SERCA, etc from Subcell" 
    self.NCX = Expression("0.")
    self.SERCA= Expression("0.")

  def Update(self,t):
    self.iryr.t = t 

    


