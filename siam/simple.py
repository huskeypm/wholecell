
from dolfin import *
from reactionsBase import ReactionsBase


class Simple(ReactionsBase): 
  def __init__(self):
    ReactionsBase.__init__(self)
    1
    
  def Init(self,params):
    self.params = params
    # RyR [A/F]???
    # term 'before' a is the dirac function, which turns on the ryr term
    # continuously after time to. The term 'after' a is the ryr term 
    dirac = "1/(1+exp(-(t-to)/v))*"
    self.iryr = Expression(dirac+"a*exp(-(t-to)/tau)",
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
    self.iNCX = Expression("0.")

    # SERCA [uM/ms]
    self.jSERCA= Expression(dirac+"-r",t=0,r=0.005,to=params.ryrOffset,v=0.05)

  def Update(self,t):
    self.iryr.t = t 

    # Need hooks for V, NCX, SERCA 
    self.jSERCA.t = t 

    


