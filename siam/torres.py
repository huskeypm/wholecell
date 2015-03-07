
from dolfin import *
from reactionsBase import ReactionsBase

class Torres(ReactionsBase): 
  def __init__(self):
    ReactionsBase.__init__(self)
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
  def Update(self,t):
    self.iryr.t = t 

