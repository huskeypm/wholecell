
from dolfin import *


class Torres(): 
  def __init__(self):
    1
    
  def Init(self,params):
    self.params = params
    # RyR [A/F]???
    # term 'before' a is the dirac function, which turns on the ryr term
    # continuously after time to. The term 'after' a is the ryr term 
    self.iryr = Expression("1/(1+exp(-(t-to)/v))*a*exp(-(t-to)/tau)",
                             v = 0.05,\
                             a=params.ryrAmp,\
                             to=params.ryrOffset,\
                             tau=params.ryrTau,\
                             t=0)
    #self.iryr = Expression("t*0", t=0)

    # Action potential [mV]
    self.V = Expression("0") 

    print "Borrow NCX, SERCA, etc from Subcell" 
    # NCX [A/F]???
    self.NCX = Expression("0.")

    # SERCA [uM/ms]
    self.SERCA= Expression("0.")

  def Update(self,t):
    self.iryr.t = t 

    # Need hooks for V, NCX, SERCA 

    


