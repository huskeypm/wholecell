
from dolfin import *
from reactionsBase import ReactionsBase


# TODO implmenet me 

# from SIAM notebook
if 0: 
    print "WARNING: need to rewrite as A/F"
    self.ryrAmp = 0.82 # [uM/ms] 
    self.ryrOffset = 5 # [ms]
    self.ryrTau = 50 # [ms]  



# uM
def SERCAExpression(cai=0.1, casr=500):
  Vmax=286. # uM/s
  Kmf = 0.246 # uM
  Kmr = 1.7e3 # mu (note: Shannon uses 1.7 mM, which I believe is for total Ca. However, free SR Ca seems to track exactly with total Ca)
  H = 1.787
  #Q10 = 2.6

  # PKH params 
  #Vmax=320. # uM/s (by setting Vmax to zero and rescaling so that maxcai/mincasr gives J=205 uM/s)
  Vmax = 3000 
  Kmf = 0.100 # uM # if set to the normal value, the pump goes in reverse
  Kmr = 0.55e3 # uM (this is my adjusted value based on free Ca)
  H = 1.0
  Q10=1 # PKH

  # cai = 0.45 #uM
  # casr = 0.55 #mM
  caiK = cai/Kmf
  caiK = caiK**H
  casrK = casr/Kmr
  casrK = casrK**H
  J = (caiK - casrK) / (1 + caiK + casrK)
  J = Q10 * Vmax * J

  # reverse so flowing INTO SR
  J = -1 * J
  #str = "-1 * Q10*Vmax*((cai/Kmf)**H - (casr/Kmr)**H)/(1+(cai/Kmf)**H - (casr/Kmr)**H)"
  str = "-1 * Q10*Vmax*(pow((cai*s/Kmf),H) - casrK)/(1+pow((cai*s/Kmf),H) + casrK)"
  str = "-1 * is_to_ims * Vmax*(cai/Kmf - casrK)/(1+cai/Kmf + casrK)"
  #str = "-1 * is_to_ims * Vmax*(x[0]/Kmf - casrK)/(1+x[0]/Kmf + casrK)"
  #str = "-1 * Q10*Vmax*(pow((cai/Kmf),H))"
  is_to_ims= 1e-3 # 1/s --> 1/ms 
  J = Expression(str,Q10=Q10,Vmax=Vmax,casrK=casrK,Kmf =Kmf,H=H,\
    cai=cai,is_to_ims = is_to_ims)
  return J



class Simple(ReactionsBase): 
  def __init__(self,noSERCA=False):
    ReactionsBase.__init__(self)
    self.noSERCA=noSERCA
    
  def Init(self,params):
    self.params = params
    # RyR [A/F]???
    # term 'before' a is the dirac function, which turns on the ryr term
    # continuously after time to. The term 'after' a is the ryr term 
    dirac = "1/(1+exp(-(t-to)/v))*"
    print params.ryrAmp
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
    self.jSERCA= SERCAExpression()
    
    if self.noSERCA:
      self.jSERCA = Expression("t*0",t=0)

  def Update(self,t,cai):
    self.iryr.t = t 

    # Need hooks for V, NCX, SERCA 
    self.jSERCA.t = t 
    self.jSERCA.cai=cai

    


