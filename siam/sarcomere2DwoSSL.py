from dolfin import * 
import numpy as np

from sarcomere2DSSL import Sarcomere2DSSL 


## var
eps = 0.05




class Sarcomere2DwoSSL(Sarcomere2DSSL):             
  def __init__(self,params=""):
    Sarcomere2DSSL.__init__(self)
    self.params = params 

