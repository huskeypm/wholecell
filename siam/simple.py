
from dolfin import *
from reactionsBase import ReactionsBase


# TODO implmenet me 

# based on Soeller but revised to match with shannon (see fitting.ipynb) 
# [1]	C. Soeller, I. D. Jayasinghe, P. Li, A. V. Holden, and M. B. Cannell, 2009
mM_to_uM = 1e3
def SERCAExpression(cai=0.1,   
  # CellML values
  Vmax= 5.311e-3*mM_to_uM,# uM/ms
  Kmf = 0.246,# uM
  #Kmr = 1.7e3 # mu (note: Shannon uses 1.7 mM, which I believe is for total Ca. However, free SR Ca seems to track exactly with total Ca)
  H = 1.787,
  casrTerm = 0.112265722393 #  math.pow(Ca_SR/Kmr, H_Jpump)  --> 0.112265722393 
  ): 
  #Q10 = 2.6
    
  # Grabbing cellML stuff   
  #T = 310
  #Ca_SR = casr  
  #Cai = cai   
  #H_Jpump = H  
  Q_SRCaP = 1 # math.pow(Q10, -31 + T/10.)
  #casr=500  

 # print casrTerm
  #j_pump_SR = Vmax*(-casrTerm +math.pow(cai/Kmf, H))*Q_SRCaP/(1 + casrTerm + math.pow(cai/Kmf, H))
  #j_pump_SR = Vmax*(-casrTerm +math.pow(cai/Kmf, H))/(1 + casrTerm + math.pow(cai/Kmf, H))

  #mystr = "-Vmax * pow(cai,H) / (pow(Kmf,H) + pow(cai,H))"
  j_pump_SR = "-Vmax*(-casrTerm +pow(cai/Kmf, H))/(1 + casrTerm + pow(cai/Kmf, H))"
  expr = Expression(j_pump_SR,
    Vmax=Vmax,
    casrTerm=casrTerm,
    Kmf=Kmf,
    cai=cai,
    H = H)
  return expr

# uM/mso
# cai [uM]
# Removed, since simplified according to strategy in fitting.ipynb
def SERCAShannonExpression(cai=0.1, casr=500):
  raise RuntimeError("Antiquated") 



class Simple(ReactionsBase): 
  def __init__(self,noSERCA=False,ryrOnlySwitch=False,caitlinSERCA=False):
    ReactionsBase.__init__(self)
    self.noSERCA=noSERCA
    self.ryrOnlySwitch=ryrOnlySwitch
    self.ryrTerminate = 10
    self.caitlinSERCA=caitlinSERCA #CES

  def Init(self,params):
    self.params = params
    # RyR [A/F]???
    # term 'before' a is the dirac function, which turns on the ryr term
    # continuously after time to. The term 'after' a is the ryr term 
    def RyRExpression(v = 0.05,a=9.5,to=1,tau=50,t=0):
      dirac = "1/(1+exp(-(t-to)/v))*"
      myiryr = Expression(dirac+"a*exp(-(t-to)/tau)",
                             v = v,a=a, to=to,tau=tau,t=t)
      return myiryr

    self.iryr = RyRExpression(a=params.ryrAmp,\
                             to=params.ryrOffset,\
                             tau=params.ryrTau,\
                             t=0)

    # shuts off after time interval (see Update) 
    if self.ryrOnlySwitch:
      self.iryr = Expression("a+t*0",a=params.ryrAmp,t=0)
      self.ryrTerminate = params.ryrTerminate
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
                                   H=params.sercaH,Kmf=params.sercaKmf,casrTerm=params.sercaCaSRTerm)
    
    if self.noSERCA:
      self.jSERCA = Expression("t*0",t=0)

  def Update(self,t,cai):
    self.iryr.t = t 

    # apply  10 ms pulse 
    if self.ryrOnlySwitch:
      if t>self.ryrTerminate: 
        self.iryr.a= 0.

    # Need hooks for V, NCX, SERCA 
    #self.jSERCA.t = t 
    self.jSERCA.cai=cai

    


