
from dolfin import *
from reactionsBase import ReactionsBase


# TODO implmenet me 

# based on Soeller 
# [1]	C. Soeller, I. D. Jayasinghe, P. Li, A. V. Holden, and M. B. Cannell, 2009
def SERCAExpression(cai=0.1,   
  Vmax = 200,# uM/Hs
  Kmf = 0.184,# uM
  H = 4.  # Hill coeff 
  ):
  mystr = "-Vmax * pow(cai,H) / (Kmf + pow(cai,H))"
  #mystr = "-0.1+0*cai"
  expr = Expression(mystr,
    Vmax=Vmax,
    Kmf=Kmf,
    cai=cai,
    H = H)
  return expr

# uM/mso
# cai [uM]
def SERCAShannonExpression(cai=0.1, casr=500):
  Vmax=286. # uM/s
  Kmf = 0.246 # uM
  Kmr = 1.7e3 # mu (note: Shannon uses 1.7 mM, which I believe is for total Ca. However, free SR Ca seems to track exactly with total Ca)
  H = 1.787
  #Q10 = 2.6

  # PKH params 
  #Vmax=320. # uM/s (by setting Vmax to zero and rescaling so that maxcai/mincasr gives J=205 uM/s)
  Vmax = 300
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


  #J = Expression("-0.00+t*0", t=0)
  return J



class Simple(ReactionsBase): 
  def __init__(self,noSERCA=False,ryrOnlySwitch=False,caitlinSERCA=False):
    ReactionsBase.__init__(self)
    self.noSERCA=noSERCA
    self.ryrOnlySwitch=ryrOnlySwitch
    self.caitlinSERCA=caitlinSERCA #CES

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

    # shuts off after time interval (see Update) 
    if self.ryrOnlySwitch:
      self.iryr = Expression("a+t*0",a=params.ryrAmp,t=0)
  #  elif self.caitlinSerca:
  #    self.iryr = Expression("0+t*0")
  #    params.flux = 0.2 # [uM/ms]
  #    self.jSERCA = Expression("a+t*0",a=params.flux,t=0)

    # Action potential [mV]
    self.V = Expression("0") 

    print "Borrow NCX, SERCA, etc from Subcell" 
    # NCX [A/F]???
    self.iNCX = Expression("0.")

    # SERCA [uM/ms]
    if self.caitlinSERCA:
      self.iryr = Expression("0+t*0",t=0)
      self.jSERCA = Expression("a+t*0",a=params.sercaVmax,t=0)
    else:
      self.jSERCA= SERCAExpression(Vmax=params.sercaVmax,\
                                   H=params.sercaH,Kmf=params.sercaKmf)
    
    if self.noSERCA:
      self.jSERCA = Expression("t*0",t=0)

  def Update(self,t,cai):
    self.iryr.t = t 

    # apply  10 ms pulse 
    if self.ryrOnlySwitch:
      if t>10: 
        self.iryr.a= 0.

    # Need hooks for V, NCX, SERCA 
    #self.jSERCA.t = t 
    self.jSERCA.cai=cai

    


