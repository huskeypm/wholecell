import numpy as np 
# In[39]:

paramList=["h","j","m","Xr","Xs","X_tos","Y_tos","R_tos","X_tof","Y_tof","d","f",
           "fCaB_SL","fCaB_jct1","R","O","I",
           "Ca_TroponinC","Ca_TroponinC_Ca_Mg",
           "Mg_TroponinC_Ca_Mg","Ca_Calmodulin","Ca_Myosin","Mg_Myosin","Ca_SRB",
           "Na_jct1_buf","Na_SL_buf","Na_jct1","Na_SL","Nai","Ca_Calsequestrin",
           "Ca_SLB_SL","Ca_SLB_jct1","Ca_SLHigh_SL","Ca_SLHigh_jct1","Ca_SR","Ca_jct1",
           "Ca_SL","Cai","V"]
paramListXX=["h","j","m","Xr","Xs","X_tos","Y_tos","R_tos","X_tof","Y_tof","d","f",
           "fCaB_SL","fCaB_jct1","R1","O1","I1",
           "Ca_TroponinC","Ca_TroponinC_Ca_Mg",
           "Mg_TroponinC_Ca_Mg","Ca_Calmodulin","Ca_Myosin","Mg_Myosin","Ca_SRB",
           "Na_jct1_buf","Na_SL_buf","Na_jct1","Na_SL","Nai","Ca_Calsequestrin",
           "Ca_SLB_SL","Ca_SLB_jct1","Ca_SLHigh_SL","Ca_SLHigh_jct1","Ca_SR","Ca_jct1",
           "Ca_SL","Cai","V"]


# In[43]:

import gotranUtil as gu
t=0; dt=400; dtn=dt/2000.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
#t=0; dt=4; dtn=dt/2000.; tsteps = np.linspace(t, t+dt, (dt)/dtn+1)
from scipy.integrate import odeint

# In[49]:

import shannon_2004_hack as model
stim_period = 400
stim_period_pIdx = model.parameter_indices("stim_period")
#states = odeint(model.rhs,s,tsteps,(p,),mxstep=1000,hmax=.03,rtol=1e-12, atol=1e-12)



## Testing 1st model 
s=model.init_state_values()
p=model.init_parameter_values()
s1 = model.rhs(s,0,p)
s1p=[]
for i,statei in enumerate(paramList):
  s1p.append( s1[model.state_indices(statei)] )

s1p = np.asarray(s1p)
#print s1p


## Testing spllit model (jct1) 
import shannon_splitcleft as msplit
s=msplit.init_state_values()
p=msplit.init_parameter_values()
stim_period_pIdx = msplit.parameter_indices("stim_period")
p[stim_period_pIdx]=stim_period
#states = odeint(msplit.rhs,s,tsteps,(p,),mxstep=1,hmax=.03,rtol=1e-12, atol=1e-12)

s2 = msplit.rhs(s,0,p)
print s2
s2p=[]
for i,statei in enumerate(paramListXX):
  s2p.append( s2[msplit.state_indices(statei)] )

s2p = np.asarray(s2p)
#print s2p

#print "___"
#print (s2p-s1p)/s1p
#bad = np.abs(s2p-s1p)/s1p>1e-6
#z = np.arange(np.shape(paramListXX)[0])[bad]
#print paramListXX[z[0]]
#print msplit.state_indices(paramListXX[z[0]])
#print paramListXX[z]

## Switch params between jct1/jct2
changedList = [
"Fx_Na_jct1",
"Fx_Na_jct2",
"Fx_NaBk_jct1",
"Fx_NaBk_jct2",
"Fx_NaK_jct1",
"Fx_NaK_jct2",
"Fx_Ks_jct1",
"Fx_Ks_jct2",
"Fx_Cl_jct1",
"Fx_Cl_jct2",
"Fx_ICaL_jct1",
"Fx_ICaL_jct2",
"Fx_NCX_jct1",
"Fx_NCX_jct2",
"Fx_SLCaP_jct1",
"Fx_SLCaP_jct2",
"Fx_CaBk_jct1",
"Fx_CaBk_jct2",
"ks1",
"ks2",
"KSRleak1",
"KSRleak2",
"Bmax_SLB_jct1",
"Bmax_SLB_jct2",
"Bmax_SLHigh_jct1",
"Bmax_SLHigh_jct2",
"Bmax_jct1",
"Bmax_jct2"]

v = np.zeros(np.shape(changedList)[0])
for i,statei in enumerate(changedList):
  v[i]=p[msplit.parameter_indices(statei)] 
  #print p[msplit.parameter_indices(statei)] 

switched = v*0
switched[1::2]=v[::2]
switched[::2] =v[1::2]
#print v
#print switched


import shannon_splitcleft as msplit
s=msplit.init_state_values()
p=msplit.init_parameter_values()
for i,statei in enumerate(changedList):
  p[msplit.parameter_indices(statei)] = switched[i]
  #print p[msplit.parameter_indices(statei)] 

s3 = msplit.rhs(s,0,p)
print s3






print "___"
print (s2-s3)
bad = np.abs(s2-s3)/s3 >1e-6

z = np.arange(np.shape(bad)[0])[bad]
print z
print ((s2-s3)/s3)[31]
#print paramListXX[z[0]]
#print msplit.state_indices(paramListXX[z[0]])
#print paramListXX[z]




