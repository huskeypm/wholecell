import matplotlib.pylab as plt
import analyze
import matplotlib
from matplotlib import pyplot
class empty:pass

job1 = empty() # create container
#job1.hdfName = "2D_SSL.h5"
#job1.hdfName = "2D_SSL_torres.h5"
job1.hdfName = "2D_SSL_simple_red_prod.h5"
#job1.hdfName = "out.h5"
job1.ts,job1.concsCaCleft,job1.concsCaSSL,job1.concsCa = analyze.ReadHdf(\
  job1.hdfName,ssl=True)


job2 = empty() # create container
#job1.hdfName = "2D_SSL.h5"
#job2.hdfName = "2D_SSL_torres_slowSSL.h5"
#job2.hdfName = "2D_noSSL_torres.h5"
job2.hdfName = "2D_noSSL_simple_red_prod.h5"
job2.ts,job2.concsCaCleft,job2.concsCaSSL,job2.concsCa = analyze.ReadHdf(\
  job2.hdfName,ssl=True)

plt.figure(num=None, figsize=(4, 3), dpi=600, facecolor='w', edgecolor='k')
plt.plot(job1.ts,job1.concsCaCleft,'r-',label="CaCleft")
plt.plot(job1.ts,job1.concsCaSSL,'b-',label="CaSSL")
print job1.concsCa[-1]
plt.plot(job1.ts,job1.concsCa,'g-',label="Ca") 

plt.ylabel("[Ca] [$\mu $m]")
plt.xlabel("time [ms]")
font = {
        'weight' : 'bold',
        'size'   : 12}
pyplot.locator_params(nbins=4)
matplotlib.rc('font', **font)
plt.legend()
plt.tight_layout()
plt.gcf().savefig("transients2D.png",dpi=600)

