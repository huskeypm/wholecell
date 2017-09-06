
import numpy as np
import matplotlib.pylab as plt 


import simpleRyR as model

idxR=model.state_indices("R")
idxO=model.state_indices("O")
idxI=model.state_indices("I")
idxjSR = 0 # model.monitor_indices("j_rel_SR")




s=model.init_state_values()
p=model.init_parameter_values()
t=0; dt=300; dtn=dt/1000.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
from scipy.integrate import odeint


ca_SR = 0.5542
ca_jct = 1.7e-4 # [mM]

ca_SR = 1e-9
ca_jct = 5e-2 # [mM]

import gotranUtil as gu

def doit(ca_jct=1.7e-4,# [mM]
         ca_SR = 0.5542 #[mM]
         ):
    s=model.init_state_values()
    p=model.init_parameter_values()
    p[model.parameter_indices("Ca_SR")] = ca_SR
    p[model.parameter_indices("Ca_jct")] = ca_jct
    mM_to_uM = 1e3
    #print "%f [uM]"%(mM_to_uM*p[model.parameter_indices("Ca_SR")])
    #print "%f [uM]"%(mM_to_uM*p[model.parameter_indices("Ca_jct")])
    states = odeint(model.rhs,s,tsteps,(p,))

    monitors = ["j_rel_SR"]
    sp,monitored = gu.RunnerMonitor(model,tsteps,s,p,monitors)
    return states,p,monitored
    
    
mM_to_uM=1e3
#jcts = np.linspace(1.7e-4,1.e-1,8)
jcts = np.linspace(1.7e-4,0.6e-1,8)
maxes = np.zeros(np.shape(jcts))
p0avgs= np.zeros(np.shape(jcts))
#jcts = np.concatenate([jcts,[1.0]])
for i,jcti in enumerate(jcts):
    #print jcti
    si,pi,mi = doit(ca_jct=jcti,# [mM]
                ca_SR = 0.5542 #[mM]
                )
    #plt.title("jsr")
    #plt.plot(tsteps[1:],mi[:,0],label="jct=%5.1f"%(mM_to_uM*jcti))
    
    plt.title("pO")
    plt.plot(tsteps,si[:,idxO],label="jct=%5.1f [uM]"%(mM_to_uM*jcti))    
    plt.xlim([0,50])
    
    # get p0 max 
    maxes[i] = np.max(si[:,idxO])

    # get integrated values
    idxtMax = np.argwhere(tsteps < 10 )[-1]
    print idxtMax,tsteps[idxtMax]    
    p0avg = np.mean( si[0:idxtMax,idxO])
    p0avgs[i] = p0avg


plt.legend()
plt.ylim([-0.005,0.03])
plt.gcf().savefig("p0.png",dpi=300)


plt.figure()
plt.plot(jcts,maxes)
plt.title("Cleft Ca vs p0 (max)") 
plt.gcf().savefig("p0max.png",dpi=300)

plt.figure()
plt.plot(jcts,p0avgs)
plt.title("Cleft Ca vs p0 (max)") 
plt.gcf().savefig("p0avg.png",dpi=300)

# should use pandas for this, but oh well
data = np.zeros([np.shape(jcts)[0],3])
data[:,0] = jcts
data[:,1] = maxes
data[:,2] = p0avgs
np.savetxt("simpleRyR.txt",data) 



#plt.figure()
#plt.title("ECC Gain")
#plt.plot(jcts,maxes)
#plt.gcf().savefig("/net/share/shared/papers/150817_satin/ecc.png",dpi=300)


