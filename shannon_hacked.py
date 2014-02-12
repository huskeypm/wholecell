# Gotran generated code for: shannon_2004

## All hacked in entries are marked with 'HACK'
## HACK
from __future__ import division
import sys
sys.path.append("/home/huskeypm/bin/Computational_Tools/serca/") 

# from Fig 10? of Cantilina paper from which Trieber model is derived     
class paramsCantilina:
    # Cantilina
    k1= 100000000.0 # 1/sec (Af)
    kn1= 400.  # 1/sec
    k2= 30.    # 1/sec (Bf) 
    kn2= 40.   # 1/sec
    k3= 400000000.0 # 1/sec (Cf) 
    kn3= 16.   # 1/sec

    K1 = kn1/k1
    K2 = kn2/k2
    K3 = kn3/k3        
    
class paramsSS:            
    K1 = paramsCantilina.kn1/paramsCantilina.k1
    K2 = paramsCantilina.kn2/paramsCantilina.k2
    K3 = paramsCantilina.kn3/paramsCantilina.k3

# mM_to_M
import sercaInesi as sI
vfunc = sI.InesiSS() # expects [M] input 
mM_to_M = 1e-3
def srMarkov(Cai_mM):
  cai_M = Cai_mM*mM_to_M
  cai_o_kmf = vfunc(cai_M,paramsCantilina)
  # not needed vals *= vmaxSI
  return cai_o_kmf 





def init_values(**values):
    """
    Init values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # h=0.9867005, j=0.991562, m=0.001405627, Xr=0.008641386, Xs=0.005412034,
    # X_tos=0.004051574, Y_tos=0.9945511, R_tos=0.9946,
    # X_tof=0.004051574, Y_tof=0.9945511, d=7.175662e-06, f=1.000681,
    # fCaB_SL=0.01452605, fCaB_jct=0.02421991, I=1.024274e-07,
    # O=8.156628e-07, R=0.8884332, Ca_Calmodulin=0.0002911916,
    # Ca_Myosin=0.001298754, Ca_SRB=0.002143165, Ca_TroponinC=0.008773191,
    # Ca_TroponinC_Ca_Mg=0.1078283, Mg_Myosin=0.1381982,
    # Mg_TroponinC_Ca_Mg=0.01524002, Na_SL=8.80733, Na_SL_buf=0.7720854,
    # Na_jct=8.80329, Na_jct_buf=3.539892, Nai=8.80853,
    # Ca_Calsequestrin=1.242988, Ca_SL=0.0001031812, Ca_SLB_SL=0.1110363,
    # Ca_SLB_jct=0.009566355, Ca_SLHigh_SL=0.07297378,
    # Ca_SLHigh_jct=0.007347888, Ca_SR=0.5545201, Ca_jct=0.0001737475,
    # Cai=8.597401e-05, V=-85.56885
    init_values = np.array([0.9867005, 0.991562, 0.001405627, 0.008641386,\
        0.005412034, 0.004051574, 0.9945511, 0.9946, 0.004051574, 0.9945511,\
        7.175662e-06, 1.000681, 0.01452605, 0.02421991, 1.024274e-07,\
        8.156628e-07, 0.8884332, 0.0002911916, 0.001298754, 0.002143165,\
        0.008773191, 0.1078283, 0.1381982, 0.01524002, 8.80733, 0.7720854,\
        8.80329, 3.539892, 8.80853, 1.242988, 0.0001031812, 0.1110363,\
        0.009566355, 0.07297378, 0.007347888, 0.5545201, 0.0001737475,\
        8.597401e-05, -85.56885], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict(h=(0, Range()), j=(1, Range()), m=(2, Range()), Xr=(3,\
        Range()), Xs=(4, Range()), X_tos=(5, Range()), Y_tos=(6, Range()),\
        R_tos=(7, Range()), X_tof=(8, Range()), Y_tof=(9, Range()), d=(10,\
        Range()), f=(11, Range()), fCaB_SL=(12, Range()), fCaB_jct=(13,\
        Range()), I=(14, Range()), O=(15, Range()), R=(16, Range()),\
        Ca_Calmodulin=(17, Range()), Ca_Myosin=(18, Range()), Ca_SRB=(19,\
        Range()), Ca_TroponinC=(20, Range()), Ca_TroponinC_Ca_Mg=(21,\
        Range()), Mg_Myosin=(22, Range()), Mg_TroponinC_Ca_Mg=(23, Range()),\
        Na_SL=(24, Range()), Na_SL_buf=(25, Range()), Na_jct=(26, Range()),\
        Na_jct_buf=(27, Range()), Nai=(28, Range()), Ca_Calsequestrin=(29,\
        Range()), Ca_SL=(30, Range()), Ca_SLB_SL=(31, Range()),\
        Ca_SLB_jct=(32, Range()), Ca_SLHigh_SL=(33, Range()),\
        Ca_SLHigh_jct=(34, Range()), Ca_SR=(35, Range()), Ca_jct=(36,\
        Range()), Cai=(37, Range()), V=(38, Range()))

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{{0}} is not a state.".format(state_name))
        ind, range = state_ind[state_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(state_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def default_parameters(**values):
    """
    Parameter values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Param values
    # Cao=1.8, Cli=15, Clo=150, Cm=1.381e-10, F=96485, Ki=135, Ko=5.4, Mgi=1,
    # Nao=140, Rgas=8314.3, T=310, cell_length=100, cell_radius=10.25,
    # Fx_Na_SL=0.89, Fx_Na_jct=0.11, G_INa=16, Fx_NaBk_SL=0.89,
    # Fx_NaBk_jct=0.11, G_NaBk=0.000297, Fx_NaK_SL=0.89, Fx_NaK_jct=0.11,
    # H_NaK=4, I_NaK_max=1.90719, Km_Ko=1.5, Km_Nai=11, Fx_Ks_SL=0.89,
    # Fx_Ks_jct=0.11, pKNa=0.01833, g_Kp=0.001, G_tos=0.06, G_tof=0.02,
    # Fx_Cl_SL=0.89, Fx_Cl_jct=0.11, G_Cl=0.109625, Kd_ClCa=0.1,
    # G_ClBk=0.009, Fx_NCX_SL=0.89, Fx_NCX_jct=0.11, HNa=3,
    # K_mCai=0.00359, K_mCao=1.3, K_mNai=12.29, K_mNao=87.5,
    # Kd_act=0.000256, Q10_NCX=1.57, V_max=9, eta=0.35, ksat=0.27,
    # Fx_SLCaP_SL=0.89, Fx_SLCaP_jct=0.11, H=1.6, Km=0.0005,
    # Q10_SLCaP=2.35, V_maxAF=0.0673, Fx_CaBk_SL=0.89, Fx_CaBk_jct=0.11,
    # G_CaBk=0.0002513, EC50_SR=0.45, HSR=2.5, Max_SR=15, Min_SR=1,
    # kiCa=0.5, kim=0.005, koCa=10, kom=0.06, ks=25,
    # KSRleak=5.348e-06, H_Jpump=1.787, Kmf=0.000246, Kmr=1.7,
    # Q10_SRCaP=2.6, V_max_Jpump=0.0053114, Bmax_Calmodulin=0.024,
    # Bmax_Myosin_Ca=0.14, Bmax_Myosin_Mg=0.14, Bmax_SRB=0.0171,
    # Bmax_TroponinC=0.07, Bmax_TroponinC_Ca_Mg_Ca=0.14,
    # Bmax_TroponinC_Ca_Mg_Mg=0.14, koff_Calmodulin=0.238,
    # koff_Myosin_Ca=0.00046, koff_Myosin_Mg=5.7e-05, koff_SRB=0.06,
    # koff_TroponinC=0.0196, koff_TroponinC_Ca_Mg_Ca=3.2e-05,
    # koff_TroponinC_Ca_Mg_Mg=0.00333, kon_Calmodulin=34,
    # kon_Myosin_Ca=13.8, kon_Myosin_Mg=0.0157, kon_SRB=100,
    # kon_TroponinC=32.7, kon_TroponinC_Ca_Mg_Ca=2.37,
    # kon_TroponinC_Ca_Mg_Mg=0.003, Fx_ICaL_SL=0.1, Fx_ICaL_jct=0.9,
    # PCa=0.00054, PK=2.7e-07, PNa=1.5e-08, Q10_CaL=1.8,
    # gamma_Cai=0.341, gamma_Cao=0.341, gamma_Ki=0.75, gamma_Ko=0.75,
    # gamma_Nai=0.75, gamma_Nao=0.75, Bmax_SL=1.65, Bmax_jct=7.561,
    # koff=0.001, kon=0.0001, Bmax_Calsequestrin=0.14,
    # Bmax_SLB_SL=0.0374, Bmax_SLB_jct=0.0046, Bmax_SLHigh_SL=0.0134,
    # Bmax_SLHigh_jct=0.00165, koff_Calsequestrin=65, koff_SLB=1.3,
    # koff_SLHigh=0.03, kon_Calsequestrin=100, kon_SL=100,
    # stim_amplitude=9.5, stim_duration=5, stim_period=1000,
    # stim_start=100
    param_values = np.array([1.8, 15, 150, 1.381e-10, 96485, 135, 5.4, 1,\
        140, 8314.3, 310, 100, 10.25, 0.89, 0.11, 16, 0.89, 0.11, 0.000297,\
        0.89, 0.11, 4, 1.90719, 1.5, 11, 0.89, 0.11, 0.01833, 0.001, 0.06,\
        0.02, 0.89, 0.11, 0.109625, 0.1, 0.009, 0.89, 0.11, 3, 0.00359, 1.3,\
        12.29, 87.5, 0.000256, 1.57, 9, 0.35, 0.27, 0.89, 0.11, 1.6, 0.0005,\
        2.35, 0.0673, 0.89, 0.11, 0.0002513, 0.45, 2.5, 15, 1, 0.5, 0.005,\
        10, 0.06, 25, 5.348e-06, 1.787, 0.000246, 1.7, 2.6, 0.0053114, 0.024,\
        0.14, 0.14, 0.0171, 0.07, 0.14, 0.14, 0.238, 0.00046, 5.7e-05, 0.06,\
        0.0196, 3.2e-05, 0.00333, 34, 13.8, 0.0157, 100, 32.7, 2.37, 0.003,\
        0.1, 0.9, 0.00054, 2.7e-07, 1.5e-08, 1.8, 0.341, 0.341, 0.75, 0.75,\
        0.75, 0.75, 1.65, 7.561, 0.001, 0.0001, 0.14, 0.0374, 0.0046, 0.0134,\
        0.00165, 65, 1.3, 0.03, 100, 100, 9.5, 5, 1000, 100],\
        dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict(Cao=(0, Range()), Cli=(1, Range()), Clo=(2, Range()),\
        Cm=(3, Range()), F=(4, Range()), Ki=(5, Range()), Ko=(6, Range()),\
        Mgi=(7, Range()), Nao=(8, Range()), Rgas=(9, Range()), T=(10,\
        Range()), cell_length=(11, Range()), cell_radius=(12, Range()),\
        Fx_Na_SL=(13, Range()), Fx_Na_jct=(14, Range()), G_INa=(15, Range()),\
        Fx_NaBk_SL=(16, Range()), Fx_NaBk_jct=(17, Range()), G_NaBk=(18,\
        Range()), Fx_NaK_SL=(19, Range()), Fx_NaK_jct=(20, Range()),\
        H_NaK=(21, Range()), I_NaK_max=(22, Range()), Km_Ko=(23, Range()),\
        Km_Nai=(24, Range()), Fx_Ks_SL=(25, Range()), Fx_Ks_jct=(26,\
        Range()), pKNa=(27, Range()), g_Kp=(28, Range()), G_tos=(29,\
        Range()), G_tof=(30, Range()), Fx_Cl_SL=(31, Range()), Fx_Cl_jct=(32,\
        Range()), G_Cl=(33, Range()), Kd_ClCa=(34, Range()), G_ClBk=(35,\
        Range()), Fx_NCX_SL=(36, Range()), Fx_NCX_jct=(37, Range()), HNa=(38,\
        Range()), K_mCai=(39, Range()), K_mCao=(40, Range()), K_mNai=(41,\
        Range()), K_mNao=(42, Range()), Kd_act=(43, Range()), Q10_NCX=(44,\
        Range()), V_max=(45, Range()), eta=(46, Range()), ksat=(47, Range()),\
        Fx_SLCaP_SL=(48, Range()), Fx_SLCaP_jct=(49, Range()), H=(50,\
        Range()), Km=(51, Range()), Q10_SLCaP=(52, Range()), V_maxAF=(53,\
        Range()), Fx_CaBk_SL=(54, Range()), Fx_CaBk_jct=(55, Range()),\
        G_CaBk=(56, Range()), EC50_SR=(57, Range()), HSR=(58, Range()),\
        Max_SR=(59, Range()), Min_SR=(60, Range()), kiCa=(61, Range()),\
        kim=(62, Range()), koCa=(63, Range()), kom=(64, Range()), ks=(65,\
        Range()), KSRleak=(66, Range()), H_Jpump=(67, Range()), Kmf=(68,\
        Range()), Kmr=(69, Range()), Q10_SRCaP=(70, Range()),\
        V_max_Jpump=(71, Range()), Bmax_Calmodulin=(72, Range()),\
        Bmax_Myosin_Ca=(73, Range()), Bmax_Myosin_Mg=(74, Range()),\
        Bmax_SRB=(75, Range()), Bmax_TroponinC=(76, Range()),\
        Bmax_TroponinC_Ca_Mg_Ca=(77, Range()), Bmax_TroponinC_Ca_Mg_Mg=(78,\
        Range()), koff_Calmodulin=(79, Range()), koff_Myosin_Ca=(80,\
        Range()), koff_Myosin_Mg=(81, Range()), koff_SRB=(82, Range()),\
        koff_TroponinC=(83, Range()), koff_TroponinC_Ca_Mg_Ca=(84, Range()),\
        koff_TroponinC_Ca_Mg_Mg=(85, Range()), kon_Calmodulin=(86, Range()),\
        kon_Myosin_Ca=(87, Range()), kon_Myosin_Mg=(88, Range()),\
        kon_SRB=(89, Range()), kon_TroponinC=(90, Range()),\
        kon_TroponinC_Ca_Mg_Ca=(91, Range()), kon_TroponinC_Ca_Mg_Mg=(92,\
        Range()), Fx_ICaL_SL=(93, Range()), Fx_ICaL_jct=(94, Range()),\
        PCa=(95, Range()), PK=(96, Range()), PNa=(97, Range()), Q10_CaL=(98,\
        Range()), gamma_Cai=(99, Range()), gamma_Cao=(100, Range()),\
        gamma_Ki=(101, Range()), gamma_Ko=(102, Range()), gamma_Nai=(103,\
        Range()), gamma_Nao=(104, Range()), Bmax_SL=(105, Range()),\
        Bmax_jct=(106, Range()), koff=(107, Range()), kon=(108, Range()),\
        Bmax_Calsequestrin=(109, Range()), Bmax_SLB_SL=(110, Range()),\
        Bmax_SLB_jct=(111, Range()), Bmax_SLHigh_SL=(112, Range()),\
        Bmax_SLHigh_jct=(113, Range()), koff_Calsequestrin=(114, Range()),\
        koff_SLB=(115, Range()), koff_SLHigh=(116, Range()),\
        kon_Calsequestrin=(117, Range()), kon_SL=(118, Range()),\
        stim_amplitude=(119, Range()), stim_duration=(120, Range()),\
        stim_period=(121, Range()), stim_start=(122, Range()))

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{{0}} is not a param".format(param_name))
        ind, range = param_ind[param_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(param_name,\
                range.format_not_in(value)))

        # Assign value
        param_values[ind] = value

    return param_values

def rhs(states, time, parameters, dy=None):
    """
    Compute right hand side
    """
    # Imports
    import numpy as np
    import math
    from math import pow, sqrt, log

    # Assign states
    assert(len(states) == 39)
    h, j, m, Xr, Xs, X_tos, Y_tos, R_tos, X_tof, Y_tof, d, f, fCaB_SL,\
        fCaB_jct, I, O, R, Ca_Calmodulin, Ca_Myosin, Ca_SRB, Ca_TroponinC,\
        Ca_TroponinC_Ca_Mg, Mg_Myosin, Mg_TroponinC_Ca_Mg, Na_SL, Na_SL_buf,\
        Na_jct, Na_jct_buf, Nai, Ca_Calsequestrin, Ca_SL, Ca_SLB_SL,\
        Ca_SLB_jct, Ca_SLHigh_SL, Ca_SLHigh_jct, Ca_SR, Ca_jct, Cai, V =\
        states

    # Assign parameters
    assert(len(parameters) == 123)
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length,\
        cell_radius, Fx_Na_SL, Fx_Na_jct, G_INa, Fx_NaBk_SL, Fx_NaBk_jct,\
        G_NaBk, Fx_NaK_SL, Fx_NaK_jct, H_NaK, I_NaK_max, Km_Ko, Km_Nai,\
        Fx_Ks_SL, Fx_Ks_jct, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct,\
        G_Cl, Kd_ClCa, G_ClBk, Fx_NCX_SL, Fx_NCX_jct, HNa, K_mCai, K_mCao,\
        K_mNai, K_mNao, Kd_act, Q10_NCX, V_max, eta, ksat, Fx_SLCaP_SL,\
        Fx_SLCaP_jct, H, Km, Q10_SLCaP, V_maxAF, Fx_CaBk_SL, Fx_CaBk_jct,\
        G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa, kom, ks,\
        KSRleak, H_Jpump, Kmf, Kmr, Q10_SRCaP, V_max_Jpump, Bmax_Calmodulin,\
        Bmax_Myosin_Ca, Bmax_Myosin_Mg, Bmax_SRB, Bmax_TroponinC,\
        Bmax_TroponinC_Ca_Mg_Ca, Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin,\
        koff_Myosin_Ca, koff_Myosin_Mg, koff_SRB, koff_TroponinC,\
        koff_TroponinC_Ca_Mg_Ca, koff_TroponinC_Ca_Mg_Mg, kon_Calmodulin,\
        kon_Myosin_Ca, kon_Myosin_Mg, kon_SRB, kon_TroponinC,\
        kon_TroponinC_Ca_Mg_Ca, kon_TroponinC_Ca_Mg_Mg, Fx_ICaL_SL,\
        Fx_ICaL_jct, PCa, PK, PNa, Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki,\
        gamma_Ko, gamma_Nai, gamma_Nao, Bmax_SL, Bmax_jct, koff, kon,\
        Bmax_Calsequestrin, Bmax_SLB_SL, Bmax_SLB_jct, Bmax_SLHigh_SL,\
        Bmax_SLHigh_jct, koff_Calsequestrin, koff_SLB, koff_SLHigh,\
        kon_Calsequestrin, kon_SL, stim_amplitude, stim_duration,\
        stim_period, stim_start = parameters

    # Model parameters
    Vol_Cell = 3.141592654e-15*(cell_radius*cell_radius)*cell_length
    Vol_SR = 0.035*Vol_Cell
    Vol_SL = 0.02*Vol_Cell
    Vol_jct = 0.000539*Vol_Cell
    Vol_myo = 0.65*Vol_Cell

    # Reversal potentials
    E_Na_jct = Rgas*T*math.log(Nao/Na_jct)/F
    E_Na_SL = Rgas*T*math.log(Nao/Na_SL)/F
    E_Ca_jct = Rgas*T*math.log(Cao/Ca_jct)/(2.0*F)
    E_Ca_SL = Rgas*T*math.log(Cao/Ca_SL)/(2.0*F)
    E_K = Rgas*T*math.log(Ko/Ki)/F
    E_Cl = Rgas*T*math.log(Cli/Clo)/F

    # Ina
    openProb = (m*m*m)*h*j
    i_Na_jct = (V - E_Na_jct)*Fx_Na_jct*G_INa*openProb
    i_Na_SL = (-E_Na_SL + V)*Fx_Na_SL*G_INa*openProb
    i_Na = i_Na_SL + i_Na_jct

    # Ina h gate
    alpha_h = (1.04951082542696e-6*math.exp(-0.147058823529412*V) if V <\
        -40.0 else 0.0)
    beta_h = (3.56*math.exp(0.079*V) + 310000.0*math.exp(0.35*V) if V < -40.0 else\
        1.0/(0.13 + 0.0497581410839387*math.exp(-0.0900900900900901*V)))

    # Ina j gate
    alpha_j = ((37.78 + V)*(-3.474e-5*math.exp(-0.04391*V) -\
        127140.0*math.exp(0.2444*V))/(1.0 +\
        50262745825.954*math.exp(0.311*V)) if V < -40.0 else 0.0)
    beta_j = (0.1212*math.exp(-0.01052*V)/(1.0 +\
        0.00396086833990426*math.exp(-0.1378*V)) if V < -40.0 else\
        0.3*math.exp(-2.535e-7*V)/(1.0 +\
        0.0407622039783662*math.exp(-0.1*V)))

    # Ina m gate
    alpha_m = (15.0816 + 0.32*V)/(1.0 - 0.00897780373069724*math.exp(-0.1*V))
    beta_m = 0.08*math.exp(-V/11.0)

    # Inab
    i_Nab_jct = (V - E_Na_jct)*Fx_NaBk_jct*G_NaBk
    i_Nab_SL = (-E_Na_SL + V)*Fx_NaBk_SL*G_NaBk
    i_Nab = i_Nab_jct + i_Nab_SL

    # Inak
    sigma = -1/7 + math.exp(0.0148588410104012*Nao)/7.0
    f_NaK = 1.0/(1.0 + 0.0365*sigma*math.exp(-F*V/(Rgas*T)) +\
        0.1245*math.exp(-0.1*F*V/(Rgas*T)))
    i_NaK_jct = Fx_NaK_jct*I_NaK_max*Ko*f_NaK/((1.0 + math.pow(Km_Nai/Na_jct,\
        H_NaK))*(Km_Ko + Ko))
    i_NaK_SL = Fx_NaK_SL*I_NaK_max*Ko*f_NaK/((1.0 + math.pow(Km_Nai/Na_SL,\
        H_NaK))*(Km_Ko + Ko))
    i_NaK = i_NaK_jct + i_NaK_SL

    # Ikr xr gate
    Xr_infinity = 1.0/(1.0 +\
        0.00127263380133981*math.exp(-0.133333333333333*V))
    tau_Xr = 1.0/((0.0061 + 0.00061*V)/(-1.0 +\
        4.26311451516882*math.exp(0.145*V)) + (0.00966 + 0.00138*V)/(1.0 -\
        0.422739131745963*math.exp(-0.123*V)))

    # Ikr rr gate
    Rr = 1.0/(1.0 + 4.36323731688614*math.exp(0.0446428571428571*V))

    # Iks
    pCa_jct = 3.0 - math.log(Ca_jct)
    pCa_SL = 3.0 - math.log(Ca_SL)
    G_Ks_jct = 0.00399 + 0.0133/(1.0 +\
        6.14421235332821e-6*math.exp(1.66666666666667*pCa_jct))
    G_Ks_SL = 0.00399 + 0.0133/(1.0 +\
        6.14421235332821e-6*math.exp(1.66666666666667*pCa_SL))
    E_Ks = Rgas*T*math.log((Ko + Nao*pKNa)/(Nai*pKNa + Ki))/F
    i_Ks_jct = (Xs*Xs)*(V - E_Ks)*Fx_Ks_jct*G_Ks_jct
    i_Ks_SL = (Xs*Xs)*(V - E_Ks)*Fx_Ks_SL*G_Ks_SL
    i_Ks = i_Ks_jct + i_Ks_SL

    # Iks xs gate
    Xs_infinity = 1.0/(1.0 + 1.09397774310453*math.exp(-0.0598802395209581*V))
    tau_Xs = 1.0/((0.002157 + 7.19e-5*V)/(1.0 -\
        0.0117959385197516*math.exp(-0.148*V)) + (0.00393 + 0.000131*V)/(-1.0 +\
        7.85381970442166*math.exp(0.0687*V)))

    # Ikp
    i_Kp = (-E_K + V)*g_Kp/(1.0 +\
        1786.47556537862*math.exp(-0.167224080267559*V))

    # Itos
    i_tos = (-E_K + V)*(0.5*R_tos + Y_tos)*G_tos*X_tos

    # Itos x gate
    X_tos_infinity = 1.0/(1.0 + math.exp(-1/5 - V/15.0))
    tau_X_tos = 0.5 + 9.0/(1.0 + math.exp(1/5 + V/15.0))

    # Itos y gate
    Y_tos_infinity = 1.0/(1.0 + 28.5027336437673*math.exp(V/10.0))
    tau_Y_tos = 30.0 + 3000.0/(1.0 + math.exp(6.0 + V/10.0))

    # Itos r gate
    R_tos_infinity = 1.0/(1.0 + 28.5027336437673*math.exp(V/10.0))
    tau_R_tos = 220.0 + 2800.0/(1.0 + math.exp(6.0 + V/10.0))

    # Itof
    i_tof = (-E_K + V)*G_tof*X_tof*Y_tof

    # Itof x gate
    X_tof_infinity = 1.0/(1.0 + math.exp(-1/5 - V/15.0))
    tau_X_tof = 1.5 + 3.5*math.exp(-(V*V)/900.0)

    # Itof y gate
    Y_tof_infinity = 1.0/(1.0 + 28.5027336437673*math.exp(V/10.0))
    tau_Y_tof = 20.0 + 20.0/(1.0 + 28.5027336437673*math.exp(V/10.0))

    # Ik1 k1 gate
    alpha_K1 = 1.02/(1.0 + 7.35454251046446e-7*math.exp(-0.2385*E_K +\
        0.2385*V))
    beta_K1 = (0.762624006506308*math.exp(0.08032*V - 0.08032*E_K) +\
        1.15340563518656e-16*math.exp(-0.06175*E_K + 0.06175*V))/(1.0 +\
        0.0867722941576933*math.exp(0.5143*E_K - 0.5143*V))
    K1_infinity = alpha_K1/(beta_K1 + alpha_K1)

    # Icl ca
    i_Cl_Ca = (-E_Cl + V)*(Fx_Cl_SL/(1.0 + Kd_ClCa/Ca_SL) + Fx_Cl_jct/(1.0 +\
        Kd_ClCa/Ca_jct))*G_Cl

    # Iclb
    i_Clb = (-E_Cl + V)*G_ClBk

    # Ical d gate
    d_infinity = 1.0/(1.0 + 0.0892185174092601*math.exp(-V/6.0))
    tau_d = (1.0 - 0.0892185174092601*math.exp(-V/6.0))*d_infinity/(0.5075 +\
        0.035*V)

    # Ical f gate
    f_infinity = 0.6/(1.0 + math.exp(5/2 - V/20.0)) + 1.0/(1.0 +\
        16964.6812589854*math.exp(0.277777777777778*V))
    tau_f = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 + 0.0337*V)*(0.48865 +\
        0.0337*V))))

    # Ical fca gate
    fCa_SL = 1.0 - fCaB_SL
    fCa_jct = 1.0 - fCaB_jct

    # Inaca
    temp_jct = (-math.pow(Nao, HNa)*Ca_jct*math.exp((-1.0 +\
        eta)*F*V/(Rgas*T)) + math.pow(Na_jct,\
        HNa)*Cao*math.exp(F*V*eta/(Rgas*T)))/(1.0 + ksat*math.exp((-1.0 +\
        eta)*F*V/(Rgas*T)))
    temp_SL = (math.pow(Na_SL, HNa)*Cao*math.exp(F*V*eta/(Rgas*T)) -\
        math.pow(Nao, HNa)*Ca_SL*math.exp((-1.0 + eta)*F*V/(Rgas*T)))/(1.0 +\
        ksat*math.exp((-1.0 + eta)*F*V/(Rgas*T)))
    Q_NCX = math.pow(Q10_NCX, -31.0 + T/10.0)
    Ka_SL = 1.0/(1.0 + (Kd_act*Kd_act*Kd_act)/(Ca_SL*Ca_SL*Ca_SL))
    Ka_jct = 1.0/(1.0 + (Kd_act*Kd_act*Kd_act)/(Ca_jct*Ca_jct*Ca_jct))
    i_NaCa_jct = Fx_NCX_jct*Ka_jct*Q_NCX*V_max*temp_jct/(math.pow(Na_jct,\
        HNa)*K_mCao + math.pow(Nao, HNa)*(1.0 + math.pow(Na_jct/K_mNai,\
        HNa))*K_mCai + math.pow(Na_jct, HNa)*Cao + math.pow(K_mNao, HNa)*(1.0 +\
        Ca_jct/K_mCai)*Ca_jct + math.pow(Nao, HNa)*Ca_jct)
    i_NaCa_SL = Fx_NCX_SL*Ka_SL*Q_NCX*V_max*temp_SL/(math.pow(Na_SL,\
        HNa)*K_mCao + math.pow(Nao, HNa)*(1.0 + math.pow(Na_SL/K_mNai,\
        HNa))*K_mCai + math.pow(Na_SL, HNa)*Cao + math.pow(K_mNao, HNa)*(1.0 +\
        Ca_SL/K_mCai)*Ca_SL + math.pow(Nao, HNa)*Ca_SL)
    i_NaCa = i_NaCa_jct + i_NaCa_SL

    # Icap
    Q_SLCaP = math.pow(Q10_SLCaP, -31.0 + T/10.0)
    i_Cap_jct = Fx_SLCaP_jct*Q_SLCaP*V_maxAF/(1.0 + math.pow(Km/Ca_jct, H))
    i_Cap_SL = Fx_SLCaP_SL*Q_SLCaP*V_maxAF/(1.0 + math.pow(Km/Ca_SL, H))
    i_Cap = i_Cap_SL + i_Cap_jct

    # Icab
    i_Cab_jct = (-E_Ca_jct + V)*Fx_CaBk_jct*G_CaBk
    i_Cab_SL = (-E_Ca_SL + V)*Fx_CaBk_SL*G_CaBk
    i_Cab = i_Cab_jct + i_Cab_SL

    # Jrel sr
    kCaSR = Max_SR - (-Min_SR + Max_SR)/(1.0 + math.pow(EC50_SR/Ca_SR, HSR))
    koSRCa = koCa/kCaSR
    kiSRCa = kCaSR*kiCa
    RI = 1.0 - I - R - O
    j_rel_SR = (Ca_SR - Ca_jct)*O*ks

    # Jleak sr
    j_leak_SR = (Ca_SR - Ca_jct)*KSRleak

    # Jpump sr
    Q_SRCaP = math.pow(Q10_SRCaP, -31.0 + T/10.0)
    j_pump_SR = (-math.pow(Ca_SR/Kmr, H_Jpump) + math.pow(Cai/Kmf,\
        H_Jpump))*Q_SRCaP*V_max_Jpump/(1.0 + math.pow(Ca_SR/Kmr, H_Jpump) +\
        math.pow(Cai/Kmf, H_Jpump))

    ## HACK 
    j_pump_SR = (-math.pow(Ca_SR/Kmr, H_Jpump) + 0*math.pow(Cai/Kmf,\
        H_Jpump))*Q_SRCaP*V_max_Jpump/(1.0 + math.pow(Ca_SR/Kmr, H_Jpump) +\
        math.pow(Cai/Kmf, H_Jpump))
    cm7n = srMarkov(Cai)
    j_pump_SR += Q_SRCaP*V_max_Jpump*cm7n 
    

    # Ion diffusion
    J_Na_jct_SL = -1.8313e-14*Na_SL + 1.8313e-14*Na_jct
    J_Na_SL_myo = -1.6386e-12*Nai + 1.6386e-12*Na_SL
    J_Ca_jct_SL = -8.2413e-13*Ca_SL + 8.2413e-13*Ca_jct
    J_Ca_SL_myo = -3.7243e-12*Cai + 3.7243e-12*Ca_SL

    # Cytosolic ca buffer
    dCa_TroponinC = (Bmax_TroponinC - Ca_TroponinC)*Cai*kon_TroponinC -\
        Ca_TroponinC*koff_TroponinC
    dCa_TroponinC_Ca_Mg = -Ca_TroponinC_Ca_Mg*koff_TroponinC_Ca_Mg_Ca +\
        (Bmax_TroponinC_Ca_Mg_Ca - Ca_TroponinC_Ca_Mg -\
        Mg_TroponinC_Ca_Mg)*Cai*kon_TroponinC_Ca_Mg_Ca
    dMg_TroponinC_Ca_Mg = -Mg_TroponinC_Ca_Mg*koff_TroponinC_Ca_Mg_Mg +\
        (Bmax_TroponinC_Ca_Mg_Mg - Ca_TroponinC_Ca_Mg -\
        Mg_TroponinC_Ca_Mg)*Mgi*kon_TroponinC_Ca_Mg_Mg
    dCa_Calmodulin = (Bmax_Calmodulin - Ca_Calmodulin)*Cai*kon_Calmodulin -\
        Ca_Calmodulin*koff_Calmodulin
    dCa_Myosin = -Ca_Myosin*koff_Myosin_Ca + (-Mg_Myosin - Ca_Myosin +\
        Bmax_Myosin_Ca)*Cai*kon_Myosin_Ca
    dMg_Myosin = -Mg_Myosin*koff_Myosin_Mg + (-Mg_Myosin - Ca_Myosin +\
        Bmax_Myosin_Mg)*Mgi*kon_Myosin_Mg
    dCa_SRB = (Bmax_SRB - Ca_SRB)*Cai*kon_SRB - Ca_SRB*koff_SRB
    dCa_cytosol_tot_bound = dCa_SRB + dCa_TroponinC_Ca_Mg + dMg_Myosin +\
        dMg_TroponinC_Ca_Mg + dCa_TroponinC + dCa_Calmodulin + dCa_Myosin

    # Ikr
    G_IKr = 0.0129099444873581*math.sqrt(Ko)
    i_Kr = (-E_K + V)*G_IKr*Rr*Xr

    # Ik1
    G_K1 = 0.387298334620742*math.sqrt(Ko)
    i_K1 = (-E_K + V)*G_K1*K1_infinity

    # Ical
    Q_CaL = math.pow(Q10_CaL, -31.0 + T/10.0)
    temp = 0.45*(F*F)*Q_CaL*V*d*f/(Rgas*T)
    i_CaL_Ca_jct = 4.0*(Ca_jct*gamma_Cai*math.exp(2.0*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*Fx_ICaL_jct*PCa*fCa_jct*temp/(-1.0 +\
        math.exp(2.0*F*V/(Rgas*T)))
    i_CaL_Na_jct = (Na_jct*gamma_Nai*math.exp(F*V/(Rgas*T)) -\
        Nao*gamma_Nao)*Fx_ICaL_jct*PNa*fCa_jct*temp/(-1.0 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_Ca_SL = 4.0*(Ca_SL*gamma_Cai*math.exp(2.0*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*Fx_ICaL_SL*PCa*fCa_SL*temp/(-1.0 +\
        math.exp(2.0*F*V/(Rgas*T)))
    i_CaL_Na_SL = (-Nao*gamma_Nao +\
        Na_SL*gamma_Nai*math.exp(F*V/(Rgas*T)))*Fx_ICaL_SL*PNa*fCa_SL*temp/(-1.0 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_K = (Fx_ICaL_SL*fCa_SL +\
        Fx_ICaL_jct*fCa_jct)*(Ki*gamma_Ki*math.exp(F*V/(Rgas*T)) -\
        Ko*gamma_Ko)*PK*temp/(-1.0 + math.exp(F*V/(Rgas*T)))
    i_CaL = i_CaL_Ca_jct + i_CaL_K + i_CaL_Ca_SL + i_CaL_Na_jct + i_CaL_Na_SL

    # Na buffer
    dNa_jct_buf = -Na_jct_buf*koff + (Bmax_jct - Na_jct_buf)*Na_jct*kon
    dNa_SL_buf = (Bmax_SL - Na_SL_buf)*Na_SL*kon - Na_SL_buf*koff

    # Ca buffer
    dCalsequestrin = -Ca_Calsequestrin*koff_Calsequestrin +\
        (Bmax_Calsequestrin*Vol_myo/Vol_SR -\
        Ca_Calsequestrin)*Ca_SR*kon_Calsequestrin
    dCa_SLB_SL = (Bmax_SLB_SL*Vol_myo/Vol_SL - Ca_SLB_SL)*Ca_SL*kon_SL -\
        Ca_SLB_SL*koff_SLB
    dCa_SLB_jct = -Ca_SLB_jct*koff_SLB + (0.1*Bmax_SLB_jct*Vol_myo/Vol_jct -\
        Ca_SLB_jct)*Ca_jct*kon_SL
    dCa_SLHigh_SL = -Ca_SLHigh_SL*koff_SLHigh + (-Ca_SLHigh_SL +\
        Bmax_SLHigh_SL*Vol_myo/Vol_SL)*Ca_SL*kon_SL
    dCa_SLHigh_jct = -Ca_SLHigh_jct*koff_SLHigh +\
        (0.1*Bmax_SLHigh_jct*Vol_myo/Vol_jct - Ca_SLHigh_jct)*Ca_jct*kon_SL
    dCa_jct_tot_bound = dCa_SLHigh_jct + dCa_SLB_jct
    dCa_SL_tot_bound = dCa_SLB_SL + dCa_SLHigh_SL
    i_Ca_jct_tot = i_Cab_jct + i_CaL_Ca_jct + i_Cap_jct - 2.0*i_NaCa_jct
    i_Ca_SL_tot = i_Cap_SL - 2.0*i_NaCa_SL + i_Cab_SL + i_CaL_Ca_SL

    # Cell
    i_Stim = (-stim_amplitude if (-stim_period*math.floor(time/stim_period) +\
        time <= stim_start + stim_duration) and\
        (-stim_period*math.floor(time/stim_period) + time >= stim_start) else\
        0.0)

    # The ODE system: 39 states

    # Init dy
    if dy is None:
        dy = np.zeros_like(states)
    dy[0] = (1.0 - h)*alpha_h - beta_h*h
    dy[1] = (1.0 - j)*alpha_j - beta_j*j
    dy[2] = (1.0 - m)*alpha_m - beta_m*m
    dy[3] = (-Xr + Xr_infinity)/tau_Xr
    dy[4] = (Xs_infinity - Xs)/tau_Xs
    dy[5] = (-X_tos + X_tos_infinity)/tau_X_tos
    dy[6] = (-Y_tos + Y_tos_infinity)/tau_Y_tos
    dy[7] = (R_tos_infinity - R_tos)/tau_R_tos
    dy[8] = (-X_tof + X_tof_infinity)/tau_X_tof
    dy[9] = (-Y_tof + Y_tof_infinity)/tau_Y_tof
    dy[10] = (-d + d_infinity)/tau_d
    dy[11] = (-f + f_infinity)/tau_f
    dy[12] = 1.7*(1.0 - fCaB_SL)*Ca_SL - 0.0119*fCaB_SL
    dy[13] = 1.7*(1.0 - fCaB_jct)*Ca_jct - 0.0119*fCaB_jct
    dy[14] = Ca_jct*O*kiSRCa + (Ca_jct*Ca_jct)*RI*koSRCa - I*kom - I*kim
    dy[15] = -Ca_jct*O*kiSRCa + I*kim + (Ca_jct*Ca_jct)*R*koSRCa - O*kom
    dy[16] = RI*kim + O*kom - Ca_jct*R*kiSRCa - (Ca_jct*Ca_jct)*R*koSRCa
    dy[17] = dCa_Calmodulin
    dy[18] = dCa_Myosin
    dy[19] = dCa_SRB
    dy[20] = dCa_TroponinC
    dy[21] = dCa_TroponinC
    dy[22] = dMg_Myosin
    dy[23] = dMg_TroponinC_Ca_Mg
    dy[24] = (-J_Na_SL_myo + J_Na_jct_SL)/Vol_SL - dNa_SL_buf -\
        (3.0*i_NaCa_SL + 3.0*i_NaK_SL + i_Nab_SL + i_Na_SL +\
        i_CaL_Na_SL)*Cm/(F*Vol_SL)
    dy[25] = dNa_SL_buf
    dy[26] = -J_Na_jct_SL/Vol_jct - (3.0*i_NaCa_jct + i_Nab_jct +\
        3.0*i_NaK_jct + i_Na_jct + i_CaL_Na_jct)*Cm/(F*Vol_jct) - dNa_jct_buf
    dy[27] = dNa_jct_buf
    dy[28] = J_Na_SL_myo/Vol_myo
    dy[29] = dCalsequestrin
    dy[30] = (-J_Ca_SL_myo + J_Ca_jct_SL)/Vol_SL - dCa_SL_tot_bound -\
        Cm*i_Ca_SL_tot/(2.0*F*Vol_SL)
    dy[31] = dCa_SLB_SL
    dy[32] = dCa_SLB_jct
    dy[33] = dCa_SLHigh_SL
    dy[34] = dCa_SLHigh_jct
    dy[35] = -j_rel_SR - Vol_myo*j_leak_SR/Vol_SR - dCalsequestrin + j_pump_SR
    dy[36] = Vol_SR*j_rel_SR/Vol_jct - dCa_jct_tot_bound -\
        J_Ca_jct_SL/Vol_jct - Cm*i_Ca_jct_tot/(2.0*F*Vol_jct) +\
        Vol_myo*j_leak_SR/Vol_jct
    dy[37] = -Vol_SR*j_pump_SR/Vol_myo - dCa_cytosol_tot_bound +\
        J_Ca_SL_myo/Vol_myo
    dy[38] = -i_NaCa - i_Cl_Ca - i_Kp - i_Stim - i_tof - i_Na - i_Cap - i_K1 -\
        i_tos - i_Cab - i_CaL - i_Ks - i_Kr - i_Clb - i_NaK - i_Nab

    # Return dy
    return dy

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict(h=0, j=1, m=2, Xr=3, Xs=4, X_tos=5, Y_tos=6, R_tos=7,\
        X_tof=8, Y_tof=9, d=10, f=11, fCaB_SL=12, fCaB_jct=13, I=14, O=15,\
        R=16, Ca_Calmodulin=17, Ca_Myosin=18, Ca_SRB=19, Ca_TroponinC=20,\
        Ca_TroponinC_Ca_Mg=21, Mg_Myosin=22, Mg_TroponinC_Ca_Mg=23, Na_SL=24,\
        Na_SL_buf=25, Na_jct=26, Na_jct_buf=27, Nai=28, Ca_Calsequestrin=29,\
        Ca_SL=30, Ca_SLB_SL=31, Ca_SLB_jct=32, Ca_SLHigh_SL=33,\
        Ca_SLHigh_jct=34, Ca_SR=35, Ca_jct=36, Cai=37, V=38)

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    return indices if len(indices)>1 else indices[0]

def param_indices(*params):
    """
    Param indices
    """
    param_inds = dict(Cao=0, Cli=1, Clo=2, Cm=3, F=4, Ki=5, Ko=6, Mgi=7,\
        Nao=8, Rgas=9, T=10, cell_length=11, cell_radius=12, Fx_Na_SL=13,\
        Fx_Na_jct=14, G_INa=15, Fx_NaBk_SL=16, Fx_NaBk_jct=17, G_NaBk=18,\
        Fx_NaK_SL=19, Fx_NaK_jct=20, H_NaK=21, I_NaK_max=22, Km_Ko=23,\
        Km_Nai=24, Fx_Ks_SL=25, Fx_Ks_jct=26, pKNa=27, g_Kp=28, G_tos=29,\
        G_tof=30, Fx_Cl_SL=31, Fx_Cl_jct=32, G_Cl=33, Kd_ClCa=34, G_ClBk=35,\
        Fx_NCX_SL=36, Fx_NCX_jct=37, HNa=38, K_mCai=39, K_mCao=40, K_mNai=41,\
        K_mNao=42, Kd_act=43, Q10_NCX=44, V_max=45, eta=46, ksat=47,\
        Fx_SLCaP_SL=48, Fx_SLCaP_jct=49, H=50, Km=51, Q10_SLCaP=52,\
        V_maxAF=53, Fx_CaBk_SL=54, Fx_CaBk_jct=55, G_CaBk=56, EC50_SR=57,\
        HSR=58, Max_SR=59, Min_SR=60, kiCa=61, kim=62, koCa=63, kom=64,\
        ks=65, KSRleak=66, H_Jpump=67, Kmf=68, Kmr=69, Q10_SRCaP=70,\
        V_max_Jpump=71, Bmax_Calmodulin=72, Bmax_Myosin_Ca=73,\
        Bmax_Myosin_Mg=74, Bmax_SRB=75, Bmax_TroponinC=76,\
        Bmax_TroponinC_Ca_Mg_Ca=77, Bmax_TroponinC_Ca_Mg_Mg=78,\
        koff_Calmodulin=79, koff_Myosin_Ca=80, koff_Myosin_Mg=81,\
        koff_SRB=82, koff_TroponinC=83, koff_TroponinC_Ca_Mg_Ca=84,\
        koff_TroponinC_Ca_Mg_Mg=85, kon_Calmodulin=86, kon_Myosin_Ca=87,\
        kon_Myosin_Mg=88, kon_SRB=89, kon_TroponinC=90,\
        kon_TroponinC_Ca_Mg_Ca=91, kon_TroponinC_Ca_Mg_Mg=92, Fx_ICaL_SL=93,\
        Fx_ICaL_jct=94, PCa=95, PK=96, PNa=97, Q10_CaL=98, gamma_Cai=99,\
        gamma_Cao=100, gamma_Ki=101, gamma_Ko=102, gamma_Nai=103,\
        gamma_Nao=104, Bmax_SL=105, Bmax_jct=106, koff=107, kon=108,\
        Bmax_Calsequestrin=109, Bmax_SLB_SL=110, Bmax_SLB_jct=111,\
        Bmax_SLHigh_SL=112, Bmax_SLHigh_jct=113, koff_Calsequestrin=114,\
        koff_SLB=115, koff_SLHigh=116, kon_Calsequestrin=117, kon_SL=118,\
        stim_amplitude=119, stim_duration=120, stim_period=121,\
        stim_start=122)

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    return indices if len(indices)>1 else indices[0]

def monitor(states, time, parameters, monitored=None,debug=False):
    """
    Compute monitored intermediates
    """
    # Imports
    import numpy as np
    import math
    from math import pow, sqrt, log

    # Assign states
    assert(len(states) == 39)
    fCaB_SL, fCaB_jct, O, Ca_TroponinC, Na_SL, Na_jct, Ca_SL, Ca_SR, Ca_jct,\
        Cai, V = states[12], states[13], states[15], states[20], states[24],\
        states[26], states[30], states[35], states[36], states[37],\
        states[38]

    # Assign parameters
    assert(len(parameters) == 123)
    Cao, F, Nao, Rgas, T, cell_length, cell_radius, Fx_NCX_SL, Fx_NCX_jct,\
        HNa, K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX, V_max, eta,\
        ksat, ks, H_Jpump, Kmf, Kmr, Q10_SRCaP, V_max_Jpump, Bmax_TroponinC,\
        koff_TroponinC, kon_TroponinC, stim_amplitude, stim_duration,\
        stim_period, stim_start = parameters[0], parameters[4],\
        parameters[8], parameters[9], parameters[10], parameters[11],\
        parameters[12], parameters[36], parameters[37], parameters[38],\
        parameters[39], parameters[40], parameters[41], parameters[42],\
        parameters[43], parameters[44], parameters[45], parameters[46],\
        parameters[47], parameters[65], parameters[67], parameters[68],\
        parameters[69], parameters[70], parameters[71], parameters[76],\
        parameters[83], parameters[90], parameters[119], parameters[120],\
        parameters[121], parameters[122]

    # Common Sub Expressions for monitored intermediates
    cse_monitored_0 = math.pow(K_mNao, HNa)
    cse_monitored_1 = (Kd_act*Kd_act*Kd_act)
    cse_monitored_2 = math.pow(Na_SL, HNa)
    cse_monitored_3 = math.pow(Na_jct, HNa)
    cse_monitored_4 = math.pow(Nao, HNa)
    cse_monitored_5 = -31.0 + T/10.0
    cse_monitored_6 = math.pow(Ca_SR/Kmr, H_Jpump)
    cse_monitored_7 = math.pow(Cai/Kmf, H_Jpump)
    cse_monitored_8 = -stim_period*math.floor(time/stim_period) + time
    cse_monitored_9 = math.exp(F*V*eta/(Rgas*T))
    cse_monitored_10 = math.exp((-1.0 + eta)*F*V/(Rgas*T))
    cse_monitored_11 = -Ca_TroponinC
    cse_monitored_12 = -Ca_jct
    cse_monitored_13 = cse_monitored_4*Ca_SL
    cse_monitored_14 = cse_monitored_2*Cao
    cse_monitored_15 = cse_monitored_3*Cao
    cse_monitored_16 = cse_monitored_4*K_mCai
    cse_monitored_17 = (cell_radius*cell_radius)*cell_length
    cse_monitored_18 = math.pow(Q10_NCX, cse_monitored_5)*V_max/(1.0 +\
        cse_monitored_10*ksat)

    # Init monitored
    if monitored is None:
        monitored = np.zeros(10, dtype=np.float_)

    # Monitored intermediates
    monitored[0] = 1.0995574289e-16*cse_monitored_17
    monitored[1] = 2.0420352251e-15*cse_monitored_17
    monitored[2] = 1.0 - fCaB_SL
    monitored[3] = 1.0 - fCaB_jct
    monitored[4] = cse_monitored_18*(-cse_monitored_10*cse_monitored_13 +\
        cse_monitored_14*cse_monitored_9)*Fx_NCX_SL/((1.0 +\
        cse_monitored_1/(Ca_SL*Ca_SL*Ca_SL))*(cse_monitored_13 +\
        cse_monitored_14 + cse_monitored_2*K_mCao + cse_monitored_0*(1.0 +\
        Ca_SL/K_mCai)*Ca_SL + cse_monitored_16*(1.0 + math.pow(Na_SL/K_mNai,\
        HNa)))) +\
        cse_monitored_18*(cse_monitored_10*cse_monitored_12*cse_monitored_4 +\
        cse_monitored_15*cse_monitored_9)*Fx_NCX_jct/((1.0 +\
        cse_monitored_1/(Ca_jct*Ca_jct*Ca_jct))*(cse_monitored_15 +\
        cse_monitored_16*(1.0 + math.pow(Na_jct/K_mNai, HNa)) +\
        cse_monitored_4*Ca_jct + cse_monitored_0*(1.0 + Ca_jct/K_mCai)*Ca_jct\
        + cse_monitored_3*K_mCao))
    monitored[5] = (Ca_SR + cse_monitored_12)*O*ks
    ## HACK
    # NOTE: can't simply replace Ca/Kmf term with Markov, since Ca/Kmf-->inf
    # while Markov normalizes to 1.  
    # for now will just ignore reverse flux of pump 
    # orig
    temp1 = math.pow(Q10_SRCaP, cse_monitored_5)*(-cse_monitored_6 +\
        cse_monitored_7)*V_max_Jpump/(1.0 + cse_monitored_7 +\
        cse_monitored_6)
    if debug:
      monitored[10] = temp1 # only works when passing in monitored[10] 

    # modified to use 'original' SB equation for reverse flux and Markov
    # for fwd 
    monitored[6] = math.pow(Q10_SRCaP, cse_monitored_5)*(-1*cse_monitored_6 +\
        0*cse_monitored_7)*V_max_Jpump/(1.0 + cse_monitored_7 +\
        cse_monitored_6)
    cm7n = srMarkov(Cai)
    monitored[6] += math.pow(Q10_SRCaP, cse_monitored_5)*V_max_Jpump*cm7n 

    monitored[7] = -3.7243e-12*Cai + 3.7243e-12*Ca_SL
    monitored[8] = cse_monitored_11*koff_TroponinC + (cse_monitored_11 +\
        Bmax_TroponinC)*Cai*kon_TroponinC
    monitored[9] = (-stim_amplitude if (cse_monitored_8 <= stim_start +\
        stim_duration) and (cse_monitored_8 >= stim_start) else 0.0)

    # Return monitored
    return (monitored)            

