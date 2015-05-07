import sys
import matplotlib.pylab import plt
sys.path.append("../")
sys.path.append("../siam/")
from sachse import *
def satinGeom(tag="satin"):
  params = Params()
  idxCaCleft = 4
  idxCaSSL   = 3
  idxCa = 0
  idxCaBuff = 1
  idxCaFluo = 2

  params.ryrOffset = 10.

  withBuffers = True
  params.T = 200
  #params.T = 50
  params.dt = 0.25 # ms

  params.D_SSLCyto = 5e1
  params.D_CleftSSL = 1e-5
  params.ryrAmp*=7.5
  reactions = "simple"

  pvdName = tag+".pvd"
  hdfName = tag+".h5"   
  concsFinal= tsolve(mode=tag,pvdName=pvdName,hdfName=hdfName,   
                       params=params,reactions=reactions,buffers=withBuffers)
  return hdfName 

def plotit(hdfName):
  import analyzePDE
  class empty:pass

  job1 = empty() # create container
  job1.hdfName = hdfName
  job1.ts,job1.concsCaCleft,job1.concsCaSSL,job1.concsCa = analyzePDE.ReadHdf(\
    job1.hdfName,ssl=True)#,verbose=True)

  plt.subplot(1,2,1)
  plt.plot(job1.ts,job1.concsCaCleft,'r-',label="CaCleft")
  plt.plot(job1.ts,job1.concsCaSSL,'b-',label="CaSSL")
  plt.plot(job1.ts,job1.concsCa,'g-',label="CaCytosol") 
  plt.xlabel("t [ms]")
  plt.legend(loc=0)
  plt.ylabel("[Ca] [uM]")
  #plt.xlim([0,100])
  
  plt.subplot(1,2,2)
  #plt.plot(job1.ts,job1.concsCaCleft,'r-',label="CaCleft")
  plt.plot(job1.ts,job1.concsCaSSL,'b-',label="CaSSL")
  plt.plot(job1.ts,job1.concsCa,'g-',label="Ca") 
  plt.xlabel("t [ms]")
  plt.tight_layout()
  
  plt.gcf().savefig("satin.png",dpi=300)


  



def doit():
  hdfName = "satin.h5"
  #hdfName = satinGeom(tag="satin")
  plotit(hdfName)

doit()
