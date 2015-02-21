# Gotran generated code for the  "shannon_splitcleft" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # h=0.9867005, j=0.991562, m=0.001405627, Xr=0.008641386, Xs=0.005412034,
    # X_tos=0.004051574, Y_tos=0.9945511, R_tos=0.9946,
    # X_tof=0.004051574, Y_tof=0.9945511, d=7.175662e-06, f=1.000681,
    # fCaB_SL=0.01452605, fCaB_jct1=0.02421991, fCaB_jct2=0.02421991,
    # R1=0.8884332, R2=0.8884332, O1=8.156628e-07, O2=8.156628e-07,
    # I1=1.024274e-07, I2=1.024274e-07, Ca_TroponinC=0.008773191,
    # Ca_TroponinC_Ca_Mg=0.1078283, Mg_TroponinC_Ca_Mg=0.01524002,
    # Ca_Calmodulin=0.0002911916, Ca_Myosin=0.001298754,
    # Mg_Myosin=0.1381982, Ca_SRB=0.002143165, Na_jct1_buf=3.539892,
    # Na_jct2_buf=3.539892, Na_SL_buf=0.7720854, Na_jct1=8.80329,
    # Na_jct2=8.80329, Na_SL=8.80733, Nai=8.80853,
    # Ca_Calsequestrin=1.242988, Ca_SLB_SL=0.1110363,
    # Ca_SLB_jct1=0.009566355, Ca_SLB_jct2=0.009566355,
    # Ca_SLHigh_SL=0.07297378, Ca_SLHigh_jct1=0.007347888,
    # Ca_SLHigh_jct2=0.007347888, Ca_SR=0.5545201, Ca_jct1=0.0001737475,
    # Ca_jct2=0.0001737475, Ca_SL=0.0001031812, Cai=8.597401e-05,
    # V=-85.56885
    init_values = np.array([0.9867005, 0.991562, 0.001405627, 0.008641386,\
        0.005412034, 0.004051574, 0.9945511, 0.9946, 0.004051574, 0.9945511,\
        7.175662e-06, 1.000681, 0.01452605, 0.02421991, 0.02421991,\
        0.8884332, 0.8884332, 8.156628e-07, 8.156628e-07, 1.024274e-07,\
        1.024274e-07, 0.008773191, 0.1078283, 0.01524002, 0.0002911916,\
        0.001298754, 0.1381982, 0.002143165, 3.539892, 3.539892, 0.7720854,\
        8.80329, 8.80329, 8.80733, 8.80853, 1.242988, 0.1110363, 0.009566355,\
        0.009566355, 0.07297378, 0.007347888, 0.007347888, 0.5545201,\
        0.0001737475, 0.0001737475, 0.0001031812, 8.597401e-05, -85.56885],\
        dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("h",(0, Range())), ("j",(1, Range())), ("m",(2,\
        Range())), ("Xr",(3, Range())), ("Xs",(4, Range())), ("X_tos",(5,\
        Range())), ("Y_tos",(6, Range())), ("R_tos",(7, Range())),\
        ("X_tof",(8, Range())), ("Y_tof",(9, Range())), ("d",(10, Range())),\
        ("f",(11, Range())), ("fCaB_SL",(12, Range())), ("fCaB_jct1",(13,\
        Range())), ("fCaB_jct2",(14, Range())), ("R1",(15, Range())),\
        ("R2",(16, Range())), ("O1",(17, Range())), ("O2",(18, Range())),\
        ("I1",(19, Range())), ("I2",(20, Range())), ("Ca_TroponinC",(21,\
        Range())), ("Ca_TroponinC_Ca_Mg",(22, Range())),\
        ("Mg_TroponinC_Ca_Mg",(23, Range())), ("Ca_Calmodulin",(24,\
        Range())), ("Ca_Myosin",(25, Range())), ("Mg_Myosin",(26, Range())),\
        ("Ca_SRB",(27, Range())), ("Na_jct1_buf",(28, Range())),\
        ("Na_jct2_buf",(29, Range())), ("Na_SL_buf",(30, Range())),\
        ("Na_jct1",(31, Range())), ("Na_jct2",(32, Range())), ("Na_SL",(33,\
        Range())), ("Nai",(34, Range())), ("Ca_Calsequestrin",(35, Range())),\
        ("Ca_SLB_SL",(36, Range())), ("Ca_SLB_jct1",(37, Range())),\
        ("Ca_SLB_jct2",(38, Range())), ("Ca_SLHigh_SL",(39, Range())),\
        ("Ca_SLHigh_jct1",(40, Range())), ("Ca_SLHigh_jct2",(41, Range())),\
        ("Ca_SR",(42, Range())), ("Ca_jct1",(43, Range())), ("Ca_jct2",(44,\
        Range())), ("Ca_SL",(45, Range())), ("Cai",(46, Range())), ("V",(47,\
        Range()))])

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind, range = state_ind[state_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(state_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def init_parameter_values(**values):
    """
    Initialize parameter values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Param values
    # Cao=1.8, Cli=15, Clo=150, Cm=1.381e-10, F=96485, Ki=135, Ko=5.4, Mgi=1,
    # Nao=140, Rgas=8314.3, T=310, cell_length=100, cell_radius=10.25,
    # Fx_Na_SL=0.89, Fx_Na_jct1=0.11, Fx_Na_jct2=0.0, G_INa=16,
    # Fx_NaBk_SL=0.89, Fx_NaBk_jct1=0.11, Fx_NaBk_jct2=0.0,
    # G_NaBk=0.000297, Fx_NaK_SL=0.89, Fx_NaK_jct1=0.11, Fx_NaK_jct2=0.0,
    # H_NaK=4, I_NaK_max=1.90719, Km_Ko=1.5, Km_Nai=11, Fx_Ks_SL=0.89,
    # Fx_Ks_jct1=0.11, Fx_Ks_jct2=0.0, pKNa=0.01833, g_Kp=0.001,
    # G_tos=0.06, G_tof=0.02, Fx_Cl_SL=0.89, Fx_Cl_jct1=0.11,
    # Fx_Cl_jct2=0.0, G_Cl=0.109625, Kd_ClCa=0.1, G_ClBk=0.009,
    # Fx_ICaL_SL=0.1, Fx_ICaL_jct1=0.9, Fx_ICaL_jct2=0.0, PCa=0.00054,
    # PK=2.7e-07, PNa=1.5e-08, Q10_CaL=1.8, gamma_Cai=0.341,
    # gamma_Cao=0.341, gamma_Ki=0.75, gamma_Ko=0.75, gamma_Nai=0.75,
    # gamma_Nao=0.75, Fx_NCX_SL=0.89, Fx_NCX_jct1=0.11, Fx_NCX_jct2=0.0,
    # HNa=3, K_mCai=0.00359, K_mCao=1.3, K_mNai=12.29, K_mNao=87.5,
    # Kd_act=0.000256, Q10_NCX=1.57, V_max_INaCa=9, eta=0.35,
    # ksat=0.27, Fx_SLCaP_SL=0.89, Fx_SLCaP_jct1=0.11, Fx_SLCaP_jct2=0.0,
    # H_ICap=1.6, Km=0.0005, Q10_SLCaP=2.35, V_maxAF=0.0673,
    # Fx_CaBk_SL=0.89, Fx_CaBk_jct1=0.11, Fx_CaBk_jct2=0.0,
    # G_CaBk=0.0002513, EC50_SR=0.45, HSR=2.5, Max_SR=15, Min_SR=1,
    # kiCa=0.5, kim=0.005, koCa=10, kom=0.06, ks1=25, ks2=0,
    # KSRleak1=5.348e-06, KSRleak2=0.0, H=1.787, Kmf=0.000246, Kmr=1.7,
    # Q10_SRCaP=2.6, V_max=0.0053114, Bmax_Calsequestrin=0.14,
    # Bmax_SLB_SL=0.0374, Bmax_SLB_jct1=0.0046, Bmax_SLB_jct2=0.0046,
    # Bmax_SLHigh_SL=0.0134, Bmax_SLHigh_jct1=0.00165,
    # Bmax_SLHigh_jct2=0.00165, koff_Calsequestrin=65, koff_SLB=1.3,
    # koff_SLHigh=0.03, kon_Calsequestrin=100, kon_SL=100,
    # Bmax_Calmodulin=0.024, Bmax_Myosin_Ca=0.14, Bmax_Myosin_Mg=0.14,
    # Bmax_SRB=0.0171, Bmax_TroponinC=0.07, Bmax_TroponinC_Ca_Mg_Ca=0.14,
    # Bmax_TroponinC_Ca_Mg_Mg=0.14, koff_Calmodulin=0.238,
    # koff_Myosin_Ca=0.00046, koff_Myosin_Mg=5.7e-05, koff_SRB=0.06,
    # koff_TroponinC=0.0196, koff_TroponinC_Ca_Mg_Ca=3.2e-05,
    # koff_TroponinC_Ca_Mg_Mg=0.00333, kon_Calmodulin=34,
    # kon_Myosin_Ca=13.8, kon_Myosin_Mg=0.0157, kon_SRB=100,
    # kon_TroponinC=32.7, kon_TroponinC_Ca_Mg_Ca=2.37,
    # kon_TroponinC_Ca_Mg_Mg=0.003, Bmax_SL=1.65, Bmax_jct1=7.561,
    # Bmax_jct2=7.561, koff=0.001, kon=0.0001, stim_amplitude=9.5,
    # stim_duration=5, stim_period=1000, stim_start=100
    init_values = np.array([1.8, 15, 150, 1.381e-10, 96485, 135, 5.4, 1, 140,\
        8314.3, 310, 100, 10.25, 0.89, 0.11, 0.0, 16, 0.89, 0.11, 0.0,\
        0.000297, 0.89, 0.11, 0.0, 4, 1.90719, 1.5, 11, 0.89, 0.11, 0.0,\
        0.01833, 0.001, 0.06, 0.02, 0.89, 0.11, 0.0, 0.109625, 0.1, 0.009,\
        0.1, 0.9, 0.0, 0.00054, 2.7e-07, 1.5e-08, 1.8, 0.341, 0.341, 0.75,\
        0.75, 0.75, 0.75, 0.89, 0.11, 0.0, 3, 0.00359, 1.3, 12.29, 87.5,\
        0.000256, 1.57, 9, 0.35, 0.27, 0.89, 0.11, 0.0, 1.6, 0.0005, 2.35,\
        0.0673, 0.89, 0.11, 0.0, 0.0002513, 0.45, 2.5, 15, 1, 0.5, 0.005, 10,\
        0.06, 25, 0, 5.348e-06, 0.0, 1.787, 0.000246, 1.7, 2.6, 0.0053114,\
        0.14, 0.0374, 0.0046, 0.0046, 0.0134, 0.00165, 0.00165, 65, 1.3,\
        0.03, 100, 100, 0.024, 0.14, 0.14, 0.0171, 0.07, 0.14, 0.14, 0.238,\
        0.00046, 5.7e-05, 0.06, 0.0196, 3.2e-05, 0.00333, 34, 13.8, 0.0157,\
        100, 32.7, 2.37, 0.003, 1.65, 7.561, 7.561, 0.001, 0.0001, 9.5, 5,\
        1000, 100], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("Cao", (0, Range())), ("Cli", (1, Range())), ("Clo",\
        (2, Range())), ("Cm", (3, Range())), ("F", (4, Range())), ("Ki", (5,\
        Range())), ("Ko", (6, Range())), ("Mgi", (7, Range())), ("Nao", (8,\
        Range())), ("Rgas", (9, Range())), ("T", (10, Range())),\
        ("cell_length", (11, Range())), ("cell_radius", (12, Range())),\
        ("Fx_Na_SL", (13, Range())), ("Fx_Na_jct1", (14, Range())),\
        ("Fx_Na_jct2", (15, Range())), ("G_INa", (16, Range())),\
        ("Fx_NaBk_SL", (17, Range())), ("Fx_NaBk_jct1", (18, Range())),\
        ("Fx_NaBk_jct2", (19, Range())), ("G_NaBk", (20, Range())),\
        ("Fx_NaK_SL", (21, Range())), ("Fx_NaK_jct1", (22, Range())),\
        ("Fx_NaK_jct2", (23, Range())), ("H_NaK", (24, Range())),\
        ("I_NaK_max", (25, Range())), ("Km_Ko", (26, Range())), ("Km_Nai",\
        (27, Range())), ("Fx_Ks_SL", (28, Range())), ("Fx_Ks_jct1", (29,\
        Range())), ("Fx_Ks_jct2", (30, Range())), ("pKNa", (31, Range())),\
        ("g_Kp", (32, Range())), ("G_tos", (33, Range())), ("G_tof", (34,\
        Range())), ("Fx_Cl_SL", (35, Range())), ("Fx_Cl_jct1", (36,\
        Range())), ("Fx_Cl_jct2", (37, Range())), ("G_Cl", (38, Range())),\
        ("Kd_ClCa", (39, Range())), ("G_ClBk", (40, Range())), ("Fx_ICaL_SL",\
        (41, Range())), ("Fx_ICaL_jct1", (42, Range())), ("Fx_ICaL_jct2",\
        (43, Range())), ("PCa", (44, Range())), ("PK", (45, Range())),\
        ("PNa", (46, Range())), ("Q10_CaL", (47, Range())), ("gamma_Cai",\
        (48, Range())), ("gamma_Cao", (49, Range())), ("gamma_Ki", (50,\
        Range())), ("gamma_Ko", (51, Range())), ("gamma_Nai", (52, Range())),\
        ("gamma_Nao", (53, Range())), ("Fx_NCX_SL", (54, Range())),\
        ("Fx_NCX_jct1", (55, Range())), ("Fx_NCX_jct2", (56, Range())),\
        ("HNa", (57, Range())), ("K_mCai", (58, Range())), ("K_mCao", (59,\
        Range())), ("K_mNai", (60, Range())), ("K_mNao", (61, Range())),\
        ("Kd_act", (62, Range())), ("Q10_NCX", (63, Range())),\
        ("V_max_INaCa", (64, Range())), ("eta", (65, Range())), ("ksat", (66,\
        Range())), ("Fx_SLCaP_SL", (67, Range())), ("Fx_SLCaP_jct1", (68,\
        Range())), ("Fx_SLCaP_jct2", (69, Range())), ("H_ICap", (70,\
        Range())), ("Km", (71, Range())), ("Q10_SLCaP", (72, Range())),\
        ("V_maxAF", (73, Range())), ("Fx_CaBk_SL", (74, Range())),\
        ("Fx_CaBk_jct1", (75, Range())), ("Fx_CaBk_jct2", (76, Range())),\
        ("G_CaBk", (77, Range())), ("EC50_SR", (78, Range())), ("HSR", (79,\
        Range())), ("Max_SR", (80, Range())), ("Min_SR", (81, Range())),\
        ("kiCa", (82, Range())), ("kim", (83, Range())), ("koCa", (84,\
        Range())), ("kom", (85, Range())), ("ks1", (86, Range())), ("ks2",\
        (87, Range())), ("KSRleak1", (88, Range())), ("KSRleak2", (89,\
        Range())), ("H", (90, Range())), ("Kmf", (91, Range())), ("Kmr", (92,\
        Range())), ("Q10_SRCaP", (93, Range())), ("V_max", (94, Range())),\
        ("Bmax_Calsequestrin", (95, Range())), ("Bmax_SLB_SL", (96,\
        Range())), ("Bmax_SLB_jct1", (97, Range())), ("Bmax_SLB_jct2", (98,\
        Range())), ("Bmax_SLHigh_SL", (99, Range())), ("Bmax_SLHigh_jct1",\
        (100, Range())), ("Bmax_SLHigh_jct2", (101, Range())),\
        ("koff_Calsequestrin", (102, Range())), ("koff_SLB", (103, Range())),\
        ("koff_SLHigh", (104, Range())), ("kon_Calsequestrin", (105,\
        Range())), ("kon_SL", (106, Range())), ("Bmax_Calmodulin", (107,\
        Range())), ("Bmax_Myosin_Ca", (108, Range())), ("Bmax_Myosin_Mg",\
        (109, Range())), ("Bmax_SRB", (110, Range())), ("Bmax_TroponinC",\
        (111, Range())), ("Bmax_TroponinC_Ca_Mg_Ca", (112, Range())),\
        ("Bmax_TroponinC_Ca_Mg_Mg", (113, Range())), ("koff_Calmodulin",\
        (114, Range())), ("koff_Myosin_Ca", (115, Range())),\
        ("koff_Myosin_Mg", (116, Range())), ("koff_SRB", (117, Range())),\
        ("koff_TroponinC", (118, Range())), ("koff_TroponinC_Ca_Mg_Ca", (119,\
        Range())), ("koff_TroponinC_Ca_Mg_Mg", (120, Range())),\
        ("kon_Calmodulin", (121, Range())), ("kon_Myosin_Ca", (122,\
        Range())), ("kon_Myosin_Mg", (123, Range())), ("kon_SRB", (124,\
        Range())), ("kon_TroponinC", (125, Range())),\
        ("kon_TroponinC_Ca_Mg_Ca", (126, Range())),\
        ("kon_TroponinC_Ca_Mg_Mg", (127, Range())), ("Bmax_SL", (128,\
        Range())), ("Bmax_jct1", (129, Range())), ("Bmax_jct2", (130,\
        Range())), ("koff", (131, Range())), ("kon", (132, Range())),\
        ("stim_amplitude", (133, Range())), ("stim_duration", (134,\
        Range())), ("stim_period", (135, Range())), ("stim_start", (136,\
        Range()))])

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind, range = param_ind[param_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(param_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict([("h", 0), ("j", 1), ("m", 2), ("Xr", 3), ("Xs", 4),\
        ("X_tos", 5), ("Y_tos", 6), ("R_tos", 7), ("X_tof", 8), ("Y_tof", 9),\
        ("d", 10), ("f", 11), ("fCaB_SL", 12), ("fCaB_jct1", 13),\
        ("fCaB_jct2", 14), ("R1", 15), ("R2", 16), ("O1", 17), ("O2", 18),\
        ("I1", 19), ("I2", 20), ("Ca_TroponinC", 21), ("Ca_TroponinC_Ca_Mg",\
        22), ("Mg_TroponinC_Ca_Mg", 23), ("Ca_Calmodulin", 24), ("Ca_Myosin",\
        25), ("Mg_Myosin", 26), ("Ca_SRB", 27), ("Na_jct1_buf", 28),\
        ("Na_jct2_buf", 29), ("Na_SL_buf", 30), ("Na_jct1", 31), ("Na_jct2",\
        32), ("Na_SL", 33), ("Nai", 34), ("Ca_Calsequestrin", 35),\
        ("Ca_SLB_SL", 36), ("Ca_SLB_jct1", 37), ("Ca_SLB_jct2", 38),\
        ("Ca_SLHigh_SL", 39), ("Ca_SLHigh_jct1", 40), ("Ca_SLHigh_jct2", 41),\
        ("Ca_SR", 42), ("Ca_jct1", 43), ("Ca_jct2", 44), ("Ca_SL", 45),\
        ("Cai", 46), ("V", 47)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict([("Cao", 0), ("Cli", 1), ("Clo", 2), ("Cm", 3), ("F",\
        4), ("Ki", 5), ("Ko", 6), ("Mgi", 7), ("Nao", 8), ("Rgas", 9), ("T",\
        10), ("cell_length", 11), ("cell_radius", 12), ("Fx_Na_SL", 13),\
        ("Fx_Na_jct1", 14), ("Fx_Na_jct2", 15), ("G_INa", 16), ("Fx_NaBk_SL",\
        17), ("Fx_NaBk_jct1", 18), ("Fx_NaBk_jct2", 19), ("G_NaBk", 20),\
        ("Fx_NaK_SL", 21), ("Fx_NaK_jct1", 22), ("Fx_NaK_jct2", 23),\
        ("H_NaK", 24), ("I_NaK_max", 25), ("Km_Ko", 26), ("Km_Nai", 27),\
        ("Fx_Ks_SL", 28), ("Fx_Ks_jct1", 29), ("Fx_Ks_jct2", 30), ("pKNa",\
        31), ("g_Kp", 32), ("G_tos", 33), ("G_tof", 34), ("Fx_Cl_SL", 35),\
        ("Fx_Cl_jct1", 36), ("Fx_Cl_jct2", 37), ("G_Cl", 38), ("Kd_ClCa",\
        39), ("G_ClBk", 40), ("Fx_ICaL_SL", 41), ("Fx_ICaL_jct1", 42),\
        ("Fx_ICaL_jct2", 43), ("PCa", 44), ("PK", 45), ("PNa", 46),\
        ("Q10_CaL", 47), ("gamma_Cai", 48), ("gamma_Cao", 49), ("gamma_Ki",\
        50), ("gamma_Ko", 51), ("gamma_Nai", 52), ("gamma_Nao", 53),\
        ("Fx_NCX_SL", 54), ("Fx_NCX_jct1", 55), ("Fx_NCX_jct2", 56), ("HNa",\
        57), ("K_mCai", 58), ("K_mCao", 59), ("K_mNai", 60), ("K_mNao", 61),\
        ("Kd_act", 62), ("Q10_NCX", 63), ("V_max_INaCa", 64), ("eta", 65),\
        ("ksat", 66), ("Fx_SLCaP_SL", 67), ("Fx_SLCaP_jct1", 68),\
        ("Fx_SLCaP_jct2", 69), ("H_ICap", 70), ("Km", 71), ("Q10_SLCaP", 72),\
        ("V_maxAF", 73), ("Fx_CaBk_SL", 74), ("Fx_CaBk_jct1", 75),\
        ("Fx_CaBk_jct2", 76), ("G_CaBk", 77), ("EC50_SR", 78), ("HSR", 79),\
        ("Max_SR", 80), ("Min_SR", 81), ("kiCa", 82), ("kim", 83), ("koCa",\
        84), ("kom", 85), ("ks1", 86), ("ks2", 87), ("KSRleak1", 88),\
        ("KSRleak2", 89), ("H", 90), ("Kmf", 91), ("Kmr", 92), ("Q10_SRCaP",\
        93), ("V_max", 94), ("Bmax_Calsequestrin", 95), ("Bmax_SLB_SL", 96),\
        ("Bmax_SLB_jct1", 97), ("Bmax_SLB_jct2", 98), ("Bmax_SLHigh_SL", 99),\
        ("Bmax_SLHigh_jct1", 100), ("Bmax_SLHigh_jct2", 101),\
        ("koff_Calsequestrin", 102), ("koff_SLB", 103), ("koff_SLHigh", 104),\
        ("kon_Calsequestrin", 105), ("kon_SL", 106), ("Bmax_Calmodulin",\
        107), ("Bmax_Myosin_Ca", 108), ("Bmax_Myosin_Mg", 109), ("Bmax_SRB",\
        110), ("Bmax_TroponinC", 111), ("Bmax_TroponinC_Ca_Mg_Ca", 112),\
        ("Bmax_TroponinC_Ca_Mg_Mg", 113), ("koff_Calmodulin", 114),\
        ("koff_Myosin_Ca", 115), ("koff_Myosin_Mg", 116), ("koff_SRB", 117),\
        ("koff_TroponinC", 118), ("koff_TroponinC_Ca_Mg_Ca", 119),\
        ("koff_TroponinC_Ca_Mg_Mg", 120), ("kon_Calmodulin", 121),\
        ("kon_Myosin_Ca", 122), ("kon_Myosin_Mg", 123), ("kon_SRB", 124),\
        ("kon_TroponinC", 125), ("kon_TroponinC_Ca_Mg_Ca", 126),\
        ("kon_TroponinC_Ca_Mg_Mg", 127), ("Bmax_SL", 128), ("Bmax_jct1",\
        129), ("Bmax_jct2", 130), ("koff", 131), ("kon", 132),\
        ("stim_amplitude", 133), ("stim_duration", 134), ("stim_period",\
        135), ("stim_start", 136)])

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict([("Vol_Cell", 0), ("Vol_SR", 1), ("Vol_SL", 2),\
        ("Vol_jct1", 3), ("Vol_jct2", 4), ("Vol_myo", 5), ("openProb", 6),\
        ("i_Na_jct1", 7), ("i_Na_jct2", 8), ("i_Na_SL", 9), ("i_Na", 10),\
        ("alpha_h", 11), ("beta_h", 12), ("alpha_j", 13), ("beta_j", 14),\
        ("alpha_m", 15), ("beta_m", 16), ("i_Nab_jct1", 17), ("i_Nab_jct2",\
        18), ("i_Nab_SL", 19), ("i_Nab", 20), ("sigma", 21), ("f_NaK", 22),\
        ("i_NaK_jct1", 23), ("i_NaK_jct2", 24), ("i_NaK_SL", 25), ("i_NaK",\
        26), ("G_IKr", 27), ("i_Kr", 28), ("Xr_infinity", 29), ("tau_Xr",\
        30), ("Rr", 31), ("pCa_jct1", 32), ("pCa_jct2", 33), ("pCa_SL", 34),\
        ("G_Ks_jct1", 35), ("G_Ks_jct2", 36), ("G_Ks_SL", 37), ("E_Ks", 38),\
        ("i_Ks_jct1", 39), ("i_Ks_jct2", 40), ("i_Ks_SL", 41), ("i_Ks", 42),\
        ("Xs_infinity", 43), ("tau_Xs", 44), ("i_Kp", 45), ("i_tos", 46),\
        ("X_tos_infinity", 47), ("tau_X_tos", 48), ("Y_tos_infinity", 49),\
        ("tau_Y_tos", 50), ("R_tos_infinity", 51), ("tau_R_tos", 52),\
        ("i_tof", 53), ("X_tof_infinity", 54), ("tau_X_tof", 55),\
        ("Y_tof_infinity", 56), ("tau_Y_tof", 57), ("i_Cl_Ca", 58), ("i_Clb",\
        59), ("Q_CaL", 60), ("temp", 61), ("i_CaL_Ca_jct1", 62),\
        ("i_CaL_Ca_jct2", 63), ("i_CaL_Na_jct1", 64), ("i_CaL_Na_jct2", 65),\
        ("i_CaL_Ca_SL", 66), ("i_CaL_Na_SL", 67), ("i_CaL_K", 68), ("i_CaL",\
        69), ("d_infinity", 70), ("tau_d", 71), ("f_infinity", 72), ("tau_f",\
        73), ("fCa_SL", 74), ("fCa_jct1", 75), ("fCa_jct2", 76),\
        ("temp_jct1", 77), ("temp_jct2", 78), ("temp_SL", 79), ("Q_NCX", 80),\
        ("Ka_SL", 81), ("Ka_jct1", 82), ("Ka_jct2", 83), ("i_NaCa_jct1", 84),\
        ("i_NaCa_jct2", 85), ("i_NaCa_SL", 86), ("i_NaCa", 87), ("Q_SLCaP",\
        88), ("i_Cap_jct1", 89), ("i_Cap_jct2", 90), ("i_Cap_SL", 91),\
        ("i_Cap", 92), ("i_Cab_jct1", 93), ("i_Cab_jct2", 94), ("i_Cab_SL",\
        95), ("i_Cab", 96), ("kCaSR", 97), ("koSRCa", 98), ("kiSRCa", 99),\
        ("RI1", 100), ("RI2", 101), ("j_rel_SR1", 102), ("j_rel_SR2", 103),\
        ("j_leak_SR1", 104), ("j_leak_SR2", 105), ("Q_SRCaP", 106),\
        ("j_pump_SR", 107), ("dCalsequestrin", 108), ("dCa_SLB_SL", 109),\
        ("dCa_SLB_jct1", 110), ("dCa_SLB_jct2", 111), ("dCa_SLHigh_SL", 112),\
        ("dCa_SLHigh_jct1", 113), ("dCa_SLHigh_jct2", 114),\
        ("dCa_jct1_tot_bound", 115), ("dCa_jct2_tot_bound", 116),\
        ("dCa_SL_tot_bound", 117), ("i_Ca_jct1_tot", 118), ("i_Ca_jct2_tot",\
        119), ("i_Ca_SL_tot", 120), ("dCa_TroponinC", 121),\
        ("dCa_TroponinC_Ca_Mg", 122), ("dMg_TroponinC_Ca_Mg", 123),\
        ("dCa_Calmodulin", 124), ("dCa_Myosin", 125), ("dMg_Myosin", 126),\
        ("dCa_SRB", 127), ("dCa_cytosol_tot_bound", 128), ("dNa_jct1_buf",\
        129), ("dNa_jct2_buf", 130), ("dNa_SL_buf", 131), ("i_Stim", 132),\
        ("E_Na_jct1", 133), ("E_Na_jct2", 134), ("E_Na_SL", 135),\
        ("E_Ca_jct1", 136), ("E_Ca_jct2", 137), ("E_Ca_SL", 138), ("E_K",\
        139), ("E_Cl", 140), ("G_K1", 141), ("i_K1", 142), ("alpha_K1", 143),\
        ("beta_K1", 144), ("K1_infinity", 145), ("J_Na_jct1_SL", 146),\
        ("J_Na_jct2_SL", 147), ("J_Na_SL_myo", 148), ("J_Ca_jct1_SL", 149),\
        ("J_Ca_jct2_SL", 150), ("J_Ca_SL_myo", 151), ("dh_dt", 152),\
        ("dj_dt", 153), ("dm_dt", 154), ("dXr_dt", 155), ("dXs_dt", 156),\
        ("dX_tos_dt", 157), ("dY_tos_dt", 158), ("dR_tos_dt", 159),\
        ("dX_tof_dt", 160), ("dY_tof_dt", 161), ("dd_dt", 162), ("df_dt",\
        163), ("dfCaB_SL_dt", 164), ("dfCaB_jct1_dt", 165), ("dfCaB_jct2_dt",\
        166), ("dR1_dt", 167), ("dR2_dt", 168), ("dO1_dt", 169), ("dO2_dt",\
        170), ("dI1_dt", 171), ("dI2_dt", 172), ("dCa_TroponinC_dt", 173),\
        ("dCa_TroponinC_Ca_Mg_dt", 174), ("dMg_TroponinC_Ca_Mg_dt", 175),\
        ("dCa_Calmodulin_dt", 176), ("dCa_Myosin_dt", 177), ("dMg_Myosin_dt",\
        178), ("dCa_SRB_dt", 179), ("dNa_jct1_buf_dt", 180),\
        ("dNa_jct2_buf_dt", 181), ("dNa_SL_buf_dt", 182), ("dNa_jct1_dt",\
        183), ("dNa_jct2_dt", 184), ("dNa_SL_dt", 185), ("dNai_dt", 186),\
        ("dCa_Calsequestrin_dt", 187), ("dCa_SLB_SL_dt", 188),\
        ("dCa_SLB_jct1_dt", 189), ("dCa_SLB_jct2_dt", 190),\
        ("dCa_SLHigh_SL_dt", 191), ("dCa_SLHigh_jct1_dt", 192),\
        ("dCa_SLHigh_jct2_dt", 193), ("dCa_SR_dt", 194), ("dCa_jct1_dt",\
        195), ("dCa_jct2_dt", 196), ("dCa_SL_dt", 197), ("dCai_dt", 198),\
        ("dV_dt", 199)])

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the shannon_splitcleft ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 48)
    h, j, m, Xr, Xs, X_tos, Y_tos, R_tos, X_tof, Y_tof, d, f, fCaB_SL,\
        fCaB_jct1, fCaB_jct2, R1, R2, O1, O2, I1, I2, Ca_TroponinC,\
        Ca_TroponinC_Ca_Mg, Mg_TroponinC_Ca_Mg, Ca_Calmodulin, Ca_Myosin,\
        Mg_Myosin, Ca_SRB, Na_jct1_buf, Na_jct2_buf, Na_SL_buf, Na_jct1,\
        Na_jct2, Na_SL, Nai, Ca_Calsequestrin, Ca_SLB_SL, Ca_SLB_jct1,\
        Ca_SLB_jct2, Ca_SLHigh_SL, Ca_SLHigh_jct1, Ca_SLHigh_jct2, Ca_SR,\
        Ca_jct1, Ca_jct2, Ca_SL, Cai, V = states

    # Assign parameters
    assert(len(parameters) == 137)
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length,\
        cell_radius, Fx_Na_SL, Fx_Na_jct1, Fx_Na_jct2, G_INa, Fx_NaBk_SL,\
        Fx_NaBk_jct1, Fx_NaBk_jct2, G_NaBk, Fx_NaK_SL, Fx_NaK_jct1,\
        Fx_NaK_jct2, H_NaK, I_NaK_max, Km_Ko, Km_Nai, Fx_Ks_SL, Fx_Ks_jct1,\
        Fx_Ks_jct2, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1,\
        Fx_Cl_jct2, G_Cl, Kd_ClCa, G_ClBk, Fx_NCX_SL, Fx_NCX_jct1,\
        Fx_NCX_jct2, HNa, K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX,\
        V_max_INaCa, eta, ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, Fx_SLCaP_jct2,\
        H_ICap, Km, Q10_SLCaP, V_maxAF, Fx_CaBk_SL, Fx_CaBk_jct1,\
        Fx_CaBk_jct2, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa,\
        kom, ks1, ks2, KSRleak1, KSRleak2, H, Kmf, Kmr, Q10_SRCaP, V_max,\
        Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg, Bmax_SRB,\
        Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca, Bmax_TroponinC_Ca_Mg_Mg,\
        koff_Calmodulin, koff_Myosin_Ca, koff_Myosin_Mg, koff_SRB,\
        koff_TroponinC, koff_TroponinC_Ca_Mg_Ca, koff_TroponinC_Ca_Mg_Mg,\
        kon_Calmodulin, kon_Myosin_Ca, kon_Myosin_Mg, kon_SRB, kon_TroponinC,\
        kon_TroponinC_Ca_Mg_Ca, kon_TroponinC_Ca_Mg_Mg, Fx_ICaL_SL,\
        Fx_ICaL_jct1, Fx_ICaL_jct2, PCa, PK, PNa, Q10_CaL, gamma_Cai,\
        gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai, gamma_Nao, Bmax_SL,\
        Bmax_jct1, Bmax_jct2, koff, kon, Bmax_Calsequestrin, Bmax_SLB_SL,\
        Bmax_SLB_jct1, Bmax_SLB_jct2, Bmax_SLHigh_SL, Bmax_SLHigh_jct1,\
        Bmax_SLHigh_jct2, koff_Calsequestrin, koff_SLB, koff_SLHigh,\
        kon_Calsequestrin, kon_SL, stim_amplitude, stim_duration,\
        stim_period, stim_start = parameters
    Cao, Cli, Clo,Cm, F, Ki,Ko, Mgi, Nao,Rgas, T,cell_length, cell_radius,Fx_Na_SL, Fx_Na_jct1,Fx_Na_jct2, G_INa,Fx_NaBk_SL, Fx_NaBk_jct1,Fx_NaBk_jct2, G_NaBk,Fx_NaK_SL, Fx_NaK_jct1,Fx_NaK_jct2, H_NaK,I_NaK_max, Km_Ko, Km_Nai,Fx_Ks_SL, Fx_Ks_jct1,Fx_Ks_jct2, pKNa,g_Kp, G_tos, G_tof,Fx_Cl_SL, Fx_Cl_jct1,Fx_Cl_jct2, G_Cl,Kd_ClCa, G_ClBk, Fx_ICaL_SL,Fx_ICaL_jct1, Fx_ICaL_jct2,PCa, PK,PNa, Q10_CaL, gamma_Cai,gamma_Cao, gamma_Ki,gamma_Ko, gamma_Nai,gamma_Nao, Fx_NCX_SL,Fx_NCX_jct1, Fx_NCX_jct2,HNa, K_mCai, K_mCao,K_mNai, K_mNao,Kd_act, Q10_NCX,V_max_INaCa, eta, ksat,Fx_SLCaP_SL, Fx_SLCaP_jct1,Fx_SLCaP_jct2, H_ICap,Km, Q10_SLCaP,V_maxAF, Fx_CaBk_SL,Fx_CaBk_jct1, Fx_CaBk_jct2,G_CaBk, EC50_SR, HSR,Max_SR, Min_SR,kiCa, kim, koCa,kom, ks1, ks2,KSRleak1, KSRleak2,H, Kmf, Kmr,Q10_SRCaP, V_max,Bmax_Calsequestrin, Bmax_SLB_SL,Bmax_SLB_jct1, Bmax_SLB_jct2,Bmax_SLHigh_SL, Bmax_SLHigh_jct1,Bmax_SLHigh_jct2,koff_Calsequestrin, koff_SLB,koff_SLHigh, kon_Calsequestrin,kon_SL, Bmax_Calmodulin,Bmax_Myosin_Ca, Bmax_Myosin_Mg,Bmax_SRB, Bmax_TroponinC,Bmax_TroponinC_Ca_Mg_Ca,Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin,koff_Myosin_Ca,koff_Myosin_Mg, koff_SRB,koff_TroponinC, koff_TroponinC_Ca_Mg_Ca,koff_TroponinC_Ca_Mg_Mg,kon_Calmodulin, kon_Myosin_Ca,kon_Myosin_Mg, kon_SRB,kon_TroponinC,kon_TroponinC_Ca_Mg_Ca,kon_TroponinC_Ca_Mg_Mg, Bmax_SL,Bmax_jct1, Bmax_jct2,koff, kon,stim_amplitude, stim_duration,stim_period, stim_start = parameters

    # Init return args
    if values is None:
        values = np.zeros((48,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (48,)

    # Expressions for the Model parameters component
    Vol_Cell = 3.141592654e-15*cell_length*(cell_radius*cell_radius)
    Vol_SR = 0.035*Vol_Cell
    Vol_SL = 0.02*Vol_Cell
    Vol_jct1 = 0.000539*Vol_Cell
    Vol_jct2 = 1
    Vol_myo = 0.65*Vol_Cell

    # Expressions for the Reversal potentials component
    E_Na_jct1 = Rgas*T*math.log(Nao/Na_jct1)/F
    E_Na_jct2 = Rgas*T*math.log(Nao/Na_jct2)/F
    E_Na_SL = Rgas*T*math.log(Nao/Na_SL)/F
    E_Ca_jct1 = Rgas*T*math.log(Cao/Ca_jct1)/(2*F)
    E_Ca_jct2 = Rgas*T*math.log(Cao/Ca_jct2)/(2*F)
    E_Ca_SL = Rgas*T*math.log(Cao/Ca_SL)/(2*F)
    E_K = Rgas*T*math.log(Ko/Ki)/F
    E_Cl = Rgas*T*math.log(Cli/Clo)/F

    # Expressions for the INa component
    openProb = (m*m*m)*h*j
    i_Na_jct1 = Fx_Na_jct1*G_INa*(V - E_Na_jct1)*openProb
    i_Na_jct2 = Fx_Na_jct2*G_INa*(-E_Na_jct2 + V)*openProb
    i_Na_SL = Fx_Na_SL*G_INa*(V - E_Na_SL)*openProb
    i_Na = i_Na_jct1 + i_Na_SL + i_Na_jct2

    # Expressions for the h gate component
    alpha_h = (1.04951082543e-06*math.exp(-0.147058823529*V) if V < -40 else 0)
    beta_h = (3.56*math.exp(0.079*V) + 310000.0*math.exp(0.35*V) if V < -40 else\
        1.0/(0.13 + 0.0497581410839*math.exp(-0.0900900900901*V)))
    values[0] = -beta_h*h + (1 - h)*alpha_h

    # Expressions for the j gate component
    alpha_j = (1.0*(37.78 + V)*(-127140.0*math.exp(0.2444*V) -\
        3.474e-05*math.exp(-0.04391*V))/(1 + 50262745826.0*math.exp(0.311*V))\
        if V < -40 else 0)
    beta_j = (0.1212*math.exp(-0.01052*V)/(1 +\
        0.0039608683399*math.exp(-0.1378*V)) if V < -40 else\
        0.3*math.exp(-2.535e-07*V)/(1 + 0.0407622039784*math.exp(-0.1*V)))
    values[1] = -beta_j*j + (1 - j)*alpha_j

    # Expressions for the m gate component
    alpha_m = (15.0816 + 0.32*V)/(1 - 0.0089778037307*math.exp(-0.1*V))
    beta_m = 0.08*math.exp(-0.0909090909091*V)
    values[2] = (1 - m)*alpha_m - beta_m*m

    # Expressions for the INab component
    i_Nab_jct1 = Fx_NaBk_jct1*G_NaBk*(V - E_Na_jct1)
    i_Nab_jct2 = Fx_NaBk_jct2*G_NaBk*(-E_Na_jct2 + V)
    i_Nab_SL = Fx_NaBk_SL*G_NaBk*(V - E_Na_SL)
    i_Nab = i_Nab_jct2 + i_Nab_jct1 + i_Nab_SL

    # Expressions for the INaK component
    sigma = -1/7. + math.exp(0.0148588410104*Nao)/7.
    f_NaK = 1.0/(1 + 0.0365*math.exp(-F*V/(Rgas*T))*sigma +\
        0.1245*math.exp(-0.1*F*V/(Rgas*T)))
    i_NaK_jct1 = Fx_NaK_jct1*I_NaK_max*Ko*f_NaK/((1 +\
        math.pow(Km_Nai/Na_jct1, H_NaK))*(Km_Ko + Ko))
    i_NaK_jct2 = Fx_NaK_jct2*I_NaK_max*Ko*f_NaK/((1 +\
        math.pow(Km_Nai/Na_jct2, H_NaK))*(Km_Ko + Ko))
    i_NaK_SL = Fx_NaK_SL*I_NaK_max*Ko*f_NaK/((1 + math.pow(Km_Nai/Na_SL,\
        H_NaK))*(Km_Ko + Ko))
    i_NaK = i_NaK_SL + i_NaK_jct1 + i_NaK_jct2

    # Expressions for the Xr gate component
    Xr_infinity = 1.0/(1 + 0.00127263380134*math.exp(-0.133333333333*V))
    tau_Xr = 1.0/((0.0061 + 0.00061*V)/(-1 + 4.26311451517*math.exp(0.145*V))\
        + (0.00966 + 0.00138*V)/(1 - 0.422739131746*math.exp(-0.123*V)))
    values[3] = (Xr_infinity - Xr)/tau_Xr

    # Expressions for the Rr gate component
    Rr = 1.0/(1 + 4.36323731689*math.exp(0.0446428571429*V))

    # Expressions for the IKs component
    pCa_jct1 = 3 - math.log(1.0*Ca_jct1)
    pCa_jct2 = 3 - math.log(1.0*Ca_jct2)
    pCa_SL = 3 - math.log(1.0*Ca_SL)
    G_Ks_jct1 = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*pCa_jct1))
    G_Ks_jct2 = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*pCa_jct2))
    G_Ks_SL = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*pCa_SL))
    E_Ks = Rgas*T*math.log((Nao*pKNa + Ko)/(pKNa*Nai + Ki))/F
    i_Ks_jct1 = Fx_Ks_jct1*(Xs*Xs)*(-E_Ks + V)*G_Ks_jct1
    i_Ks_jct2 = Fx_Ks_jct2*(Xs*Xs)*(-E_Ks + V)*G_Ks_jct2
    i_Ks_SL = Fx_Ks_SL*(Xs*Xs)*(-E_Ks + V)*G_Ks_SL
    i_Ks = i_Ks_jct1 + i_Ks_SL + i_Ks_jct2

    # Expressions for the Xs gate component
    Xs_infinity = 1.0/(1 + 1.0939777431*math.exp(-0.059880239521*V))
    tau_Xs = 1.0/((0.002157 + 7.19e-05*V)/(1 -\
        0.0117959385198*math.exp(-0.148*V)) + (0.00393 + 0.000131*V)/(-1 +\
        7.85381970442*math.exp(0.0687*V)))
    values[4] = (-Xs + Xs_infinity)/tau_Xs

    # Expressions for the IKp component
    i_Kp = g_Kp*(V - E_K)/(1 + 1786.47556538*math.exp(-0.167224080268*V))

    # Expressions for the Itos component
    i_tos = G_tos*(0.5*R_tos + Y_tos)*(V - E_K)*X_tos

    # Expressions for the X_gate component
    #X_tos_infinity = 1.0/(1 + 0.818730753078*math.exp(-0.0666666666667*V))
    #tau_X_tos = 0.5 + 9/(1 + 1.22140275816*math.exp(0.0666666666667*V))
    X_tos_infinity = 1.0/(1 + math.exp(-1/5. - V/15.))
    tau_X_tos = 0.5 + 9/(1 + math.exp(1/5. + V/15.))
    values[5] = (-X_tos + X_tos_infinity)/tau_X_tos

    # Expressions for the Y_gate component
    Y_tos_infinity = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    tau_Y_tos = 30 + 3000/(1 + 403.428793493*math.exp(0.1*V))
    values[6] = (Y_tos_infinity - Y_tos)/tau_Y_tos

    # Expressions for the R_gate component
    R_tos_infinity = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    tau_R_tos = 220 + 2800.0/(1 + 403.428793493*math.exp(0.1*V))
    values[7] = (-R_tos + R_tos_infinity)/tau_R_tos

    # Expressions for the Itof component
    i_tof = G_tof*(V - E_K)*X_tof*Y_tof

    # Expressions for the Itof X gate component
    #X_tof_infinity = 1.0/(1 + 0.818730753078*math.exp(-0.0666666666667*V))
    #tau_X_tof = 1.5 + 3.5*math.exp(-0.00111111111111*(V*V))
    X_tof_infinity = 1.0/(1 + math.exp(-1/5. - V/15.))
    tau_X_tof = 1.5 + 3.5*math.exp(-(V*V)/900.)
    values[8] = (-X_tof + X_tof_infinity)/tau_X_tof

    # Expressions for the Itof Y gate component
    Y_tof_infinity = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    tau_Y_tof = 20 + 20/(1 + 28.5027336438*math.exp(0.1*V))
    values[9] = (Y_tof_infinity - Y_tof)/tau_Y_tof

    # Expressions for the K1 gate component
    alpha_K1 = 1.02/(1 + 7.35454251046e-07*math.exp(-0.2385*E_K + 0.2385*V))
    beta_K1 = (0.762624006506*math.exp(-0.08032*E_K + 0.08032*V) +\
        1.15340563519e-16*math.exp(-0.06175*E_K + 0.06175*V))/(1 +\
        0.0867722941577*math.exp(0.5143*E_K - 0.5143*V))
    K1_infinity = alpha_K1/(beta_K1 + alpha_K1)

    # Expressions for the ICl Ca component
    i_Cl_Ca = G_Cl*(V - E_Cl)*(Fx_Cl_SL/(1 + Kd_ClCa/Ca_SL) + Fx_Cl_jct2/(1 +\
        Kd_ClCa/Ca_jct2) + Fx_Cl_jct1/(1 + Kd_ClCa/Ca_jct1))

    # Expressions for the IClb component
    i_Clb = G_ClBk*(V - E_Cl)

    # Expressions for the d gate component
    #d_infinity = 1.0/(1 + 0.0892185174093*math.exp(-0.166666666667*V))
    #tau_d = (1 -\
    #    0.0892185174093*math.exp(-0.166666666667*V))*d_infinity/(0.5075 +\
    #    0.035*V)
    d_infinity = 1.0/(1 + 0.0892185174093*math.exp(-V/6.))
    tau_d = (1 - 0.0892185174093*math.exp(-V/6.))*d_infinity/(0.5075 + 0.035*V)
    values[10] = (d_infinity - d)/tau_d

    # Expressions for the f gate component
    f_infinity = 1.0/(1 + 16964.681259*math.exp(0.277777777778*V)) + 0.6/(1 +\
        12.1824939607*math.exp(-0.05*V))
    tau_f = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 + 0.0337*V)*(0.48865 +\
        0.0337*V))))
    values[11] = (f_infinity - f)/tau_f

    # Expressions for the FCa gate component
    fCa_SL = 1 - fCaB_SL
    fCa_jct1 = 1 - fCaB_jct1
    fCa_jct2 = 1 - fCaB_jct2
    values[12] = -0.0119*fCaB_SL + 1.7*(1 - fCaB_SL)*Ca_SL
    values[13] = -0.0119*fCaB_jct1 + 1.7*(1 - fCaB_jct1)*Ca_jct1
    values[14] = 1.7*(1 - fCaB_jct2)*Ca_jct2 - 0.0119*fCaB_jct2

    # Expressions for the INaCa component
    temp_jct1 = (-math.pow(Nao, HNa)*Ca_jct1*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct1,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    temp_jct2 = (-math.pow(Nao, HNa)*Ca_jct2*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct2,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    temp_SL = (-math.pow(Nao, HNa)*Ca_SL*math.exp(F*(-1 + eta)*V/(Rgas*T)) +\
        Cao*math.pow(Na_SL, HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 +\
        ksat*math.exp(F*(-1 + eta)*V/(Rgas*T)))
    Q_NCX = math.pow(Q10_NCX, -31.0 + 0.1*T)
    Ka_SL = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_SL*Ca_SL*Ca_SL))
    Ka_jct1 = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct1*Ca_jct1*Ca_jct1))
    Ka_jct2 = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct2*Ca_jct2*Ca_jct2))
    i_NaCa_jct1 =\
        Fx_NCX_jct1*V_max_INaCa*Ka_jct1*Q_NCX*temp_jct1/(K_mCai*math.pow(Nao,\
        HNa)*(1 + math.pow(Na_jct1/K_mNai, HNa)) + K_mCao*math.pow(Na_jct1,\
        HNa) + Cao*math.pow(Na_jct1, HNa) + math.pow(Nao, HNa)*Ca_jct1 +\
        math.pow(K_mNao, HNa)*(1 + Ca_jct1/K_mCai)*Ca_jct1)
    i_NaCa_jct2 =\
        Fx_NCX_jct2*V_max_INaCa*Ka_jct2*Q_NCX*temp_jct2/(math.pow(K_mNao,\
        HNa)*(1 + Ca_jct2/K_mCai)*Ca_jct2 + K_mCai*math.pow(Nao, HNa)*(1 +\
        math.pow(Na_jct2/K_mNai, HNa)) + Cao*math.pow(Na_jct2, HNa) +\
        math.pow(Nao, HNa)*Ca_jct2 + K_mCao*math.pow(Na_jct2, HNa))
    i_NaCa_SL = Fx_NCX_SL*V_max_INaCa*Ka_SL*Q_NCX*temp_SL/(math.pow(K_mNao,\
        HNa)*(1 + Ca_SL/K_mCai)*Ca_SL + math.pow(Nao, HNa)*Ca_SL +\
        Cao*math.pow(Na_SL, HNa) + K_mCao*math.pow(Na_SL, HNa) +\
        K_mCai*math.pow(Nao, HNa)*(1 + math.pow(Na_SL/K_mNai, HNa)))
    i_NaCa = i_NaCa_jct2 + i_NaCa_SL + i_NaCa_jct1

    # Expressions for the ICap component
    Q_SLCaP = math.pow(Q10_SLCaP, -31.0 + 0.1*T)
    i_Cap_jct1 = Fx_SLCaP_jct1*V_maxAF*Q_SLCaP/(1 + math.pow(Km/Ca_jct1,\
        H_ICap))
    i_Cap_jct2 = Fx_SLCaP_jct2*V_maxAF*Q_SLCaP/(1 + math.pow(Km/Ca_jct2,\
        H_ICap))
    i_Cap_SL = Fx_SLCaP_SL*V_maxAF*Q_SLCaP/(1 + math.pow(Km/Ca_SL, H_ICap))
    i_Cap = i_Cap_SL + i_Cap_jct2 + i_Cap_jct1

    # Expressions for the ICab component
    i_Cab_jct1 = Fx_CaBk_jct1*G_CaBk*(-E_Ca_jct1 + V)
    i_Cab_jct2 = Fx_CaBk_jct2*G_CaBk*(V - E_Ca_jct2)
    i_Cab_SL = Fx_CaBk_SL*G_CaBk*(V - E_Ca_SL)
    i_Cab = i_Cab_jct1 + i_Cab_jct2 + i_Cab_SL

    # Expressions for the Jrel SR component
    kCaSR = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR, HSR))
    koSRCa = koCa/kCaSR
    kiSRCa = kiCa*kCaSR
    RI1 = 1 - I1 - O1 - R1
    RI2 = 1 - R2 - O2 - I2
    values[15] = kom*O1 + kim*RI1 - (Ca_jct1*Ca_jct1)*R1*koSRCa -\
        Ca_jct1*R1*kiSRCa
    values[16] = -(Ca_jct2*Ca_jct2)*R2*koSRCa - Ca_jct2*R2*kiSRCa + kim*RI2 +\
        kom*O2
    values[17] = kim*I1 + (Ca_jct1*Ca_jct1)*R1*koSRCa - Ca_jct1*O1*kiSRCa -\
        kom*O1
    values[18] = -Ca_jct2*O2*kiSRCa + kim*I2 - kom*O2 +\
        (Ca_jct2*Ca_jct2)*R2*koSRCa
    values[19] = Ca_jct1*O1*kiSRCa - kom*I1 - kim*I1 +\
        (Ca_jct1*Ca_jct1)*RI1*koSRCa
    values[20] = Ca_jct2*O2*kiSRCa - kim*I2 + (Ca_jct2*Ca_jct2)*RI2*koSRCa -\
        kom*I2
    j_rel_SR1 = ks1*(Ca_SR - Ca_jct1)*O1
    j_rel_SR2 = 0*ks2*(Ca_SR - Ca_jct2)*O2

    # Expressions for the Jleak SR component
    j_leak_SR1 = KSRleak1*(Ca_SR - Ca_jct1)
    j_leak_SR2 = 0*KSRleak2*(Ca_SR - Ca_jct2)

    # Expressions for the Jpump SR component
    Q_SRCaP = math.pow(Q10_SRCaP, -31.0 + 0.1*T)
    j_pump_SR = V_max*(math.pow(Cai/Kmf, H) - math.pow(Ca_SR/Kmr,\
        H))*Q_SRCaP/(1 + math.pow(Cai/Kmf, H) + math.pow(Ca_SR/Kmr, H))

    # Expressions for the Ion diffusion component
    J_Na_jct1_SL = -1.8313e-14*Na_SL + 1.8313e-14*Na_jct1
    J_Na_jct2_SL = 0
    J_Na_SL_myo = 1.6386e-12*Na_SL - 1.6386e-12*Nai
    J_Ca_jct1_SL = 8.2413e-13*Ca_jct1 - 8.2413e-13*Ca_SL
    J_Ca_jct2_SL = 0
    J_Ca_SL_myo = -3.7243e-12*Cai + 3.7243e-12*Ca_SL

    # Expressions for the Cytosolic component
    dCa_TroponinC = kon_TroponinC*(-Ca_TroponinC + Bmax_TroponinC)*Cai -\
        koff_TroponinC*Ca_TroponinC
    dCa_TroponinC_Ca_Mg = -koff_TroponinC_Ca_Mg_Ca*Ca_TroponinC_Ca_Mg +\
        kon_TroponinC_Ca_Mg_Ca*(-Ca_TroponinC_Ca_Mg - Mg_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Ca)*Cai
    dMg_TroponinC_Ca_Mg = Mgi*kon_TroponinC_Ca_Mg_Mg*(-Ca_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Mg - Mg_TroponinC_Ca_Mg) -\
        koff_TroponinC_Ca_Mg_Mg*Mg_TroponinC_Ca_Mg
    dCa_Calmodulin = -koff_Calmodulin*Ca_Calmodulin +\
        kon_Calmodulin*(Bmax_Calmodulin - Ca_Calmodulin)*Cai
    dCa_Myosin = kon_Myosin_Ca*(Bmax_Myosin_Ca - Mg_Myosin - Ca_Myosin)*Cai -\
        koff_Myosin_Ca*Ca_Myosin
    dMg_Myosin = Mgi*kon_Myosin_Mg*(Bmax_Myosin_Mg - Mg_Myosin - Ca_Myosin) -\
        koff_Myosin_Mg*Mg_Myosin
    dCa_SRB = kon_SRB*(-Ca_SRB + Bmax_SRB)*Cai - koff_SRB*Ca_SRB
    dCa_cytosol_tot_bound = dCa_SRB + dCa_Calmodulin + dMg_Myosin +\
        dMg_TroponinC_Ca_Mg + dCa_TroponinC + dCa_TroponinC_Ca_Mg +\
        dCa_Myosin
    values[21] = dCa_TroponinC
    values[22] = dCa_TroponinC_Ca_Mg
    values[23] = dMg_TroponinC_Ca_Mg
    values[24] = dCa_Calmodulin
    values[25] = dCa_Myosin
    values[26] = dMg_Myosin
    values[27] = dCa_SRB

    # Expressions for the IKr component
    G_IKr = 0.0129099444874*math.sqrt(Ko)
    i_Kr = (V - E_K)*G_IKr*Rr*Xr

    # Expressions for the IK1 component
    G_K1 = 0.387298334621*math.sqrt(Ko)
    i_K1 = (V - E_K)*G_K1*K1_infinity

    # Expressions for the ICaL component
    Q_CaL = math.pow(Q10_CaL, -31.0 + 0.1*T)
    temp = 0.45*(F*F)*Q_CaL*V*d*f/(Rgas*T)
    i_CaL_Ca_jct1 =\
        4*Fx_ICaL_jct1*PCa*(gamma_Cai*Ca_jct1*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*fCa_jct1*temp/(-1 + math.exp(2*F*V/(Rgas*T)))
    i_CaL_Ca_jct2 =\
        4*Fx_ICaL_jct2*PCa*(gamma_Cai*Ca_jct2*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*fCa_jct2*temp/(-1 + math.exp(2*F*V/(Rgas*T)))
    i_CaL_Na_jct1 = Fx_ICaL_jct1*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct1*math.exp(F*V/(Rgas*T)))*fCa_jct1*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_Na_jct2 = Fx_ICaL_jct2*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct2*math.exp(F*V/(Rgas*T)))*fCa_jct2*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_Ca_SL = 4*Fx_ICaL_SL*PCa*(gamma_Cai*Ca_SL*math.exp(2*F*V/(Rgas*T))\
        - Cao*gamma_Cao)*fCa_SL*temp/(-1 + math.exp(2*F*V/(Rgas*T)))
    i_CaL_Na_SL = Fx_ICaL_SL*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_SL*math.exp(F*V/(Rgas*T)))*fCa_SL*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_K = PK*(-Ko*gamma_Ko +\
        Ki*gamma_Ki*math.exp(F*V/(Rgas*T)))*(Fx_ICaL_SL*fCa_SL +\
        Fx_ICaL_jct1*fCa_jct1 + Fx_ICaL_jct2*fCa_jct2)*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL = i_CaL_Ca_jct1 + i_CaL_Ca_SL + i_CaL_Na_jct1 + i_CaL_K +\
        i_CaL_Na_jct2 + i_CaL_Na_SL + i_CaL_Ca_jct2

    # Expressions for the Na buffer component
    dNa_jct1_buf = kon*(Bmax_jct1 - Na_jct1_buf)*Na_jct1 - koff*Na_jct1_buf
    dNa_jct2_buf = kon*(Bmax_jct2 - Na_jct2_buf)*Na_jct2 - koff*Na_jct2_buf
    dNa_SL_buf = -koff*Na_SL_buf + kon*(-Na_SL_buf + Bmax_SL)*Na_SL
    values[28] = dNa_jct1_buf
    values[29] = dNa_jct2_buf
    values[30] = dNa_SL_buf
    values[31] = -dNa_jct1_buf - Cm*(3*i_NaCa_jct1 + 3*i_NaK_jct1 + i_Na_jct1 +\
        i_CaL_Na_jct1 + i_Nab_jct1)/(F*Vol_jct1) - J_Na_jct1_SL/Vol_jct1
    values[32] = -J_Na_jct2_SL/Vol_jct2 - dNa_jct2_buf - Cm*(3*i_NaCa_jct2 +\
        i_Nab_jct2 + 3*i_NaK_jct2 + i_CaL_Na_jct2 + i_Na_jct2)/(F*Vol_jct2)
    values[33] = -dNa_SL_buf - Cm*(3*i_NaCa_SL + i_Na_SL + i_CaL_Na_SL +\
        i_Nab_SL + 3*i_NaK_SL)/(F*Vol_SL) + (J_Na_jct2_SL - J_Na_SL_myo +\
        J_Na_jct1_SL)/Vol_SL
    values[34] = J_Na_SL_myo/Vol_myo

    # Expressions for the Ca buffer component
    dCalsequestrin = -koff_Calsequestrin*Ca_Calsequestrin +\
        kon_Calsequestrin*(Bmax_Calsequestrin*Vol_myo/Vol_SR -\
        Ca_Calsequestrin)*Ca_SR
    values[35] = dCalsequestrin
    dCa_SLB_SL = kon_SL*(-Ca_SLB_SL + Bmax_SLB_SL*Vol_myo/Vol_SL)*Ca_SL -\
        koff_SLB*Ca_SLB_SL
    dCa_SLB_jct1 = -koff_SLB*Ca_SLB_jct1 +\
        kon_SL*(0.1*Bmax_SLB_jct1*Vol_myo/Vol_jct1 - Ca_SLB_jct1)*Ca_jct1
    dCa_SLB_jct2 = -koff_SLB*Ca_SLB_jct2 +\
        kon_SL*(0.1*Bmax_SLB_jct2*Vol_myo/Vol_jct2 - Ca_SLB_jct2)*Ca_jct2
    dCa_SLHigh_SL = -koff_SLHigh*Ca_SLHigh_SL +\
        kon_SL*(Bmax_SLHigh_SL*Vol_myo/Vol_SL - Ca_SLHigh_SL)*Ca_SL
    dCa_SLHigh_jct1 = -koff_SLHigh*Ca_SLHigh_jct1 +\
        kon_SL*(0.1*Bmax_SLHigh_jct1*Vol_myo/Vol_jct1 -\
        Ca_SLHigh_jct1)*Ca_jct1
    dCa_SLHigh_jct2 = kon_SL*(-Ca_SLHigh_jct2 +\
        0.1*Bmax_SLHigh_jct2*Vol_myo/Vol_jct2)*Ca_jct2 -\
        koff_SLHigh*Ca_SLHigh_jct2
    values[36] = dCa_SLB_SL
    values[37] = dCa_SLB_jct1
    values[38] = dCa_SLB_jct2
    values[39] = dCa_SLHigh_SL
    values[40] = dCa_SLHigh_jct1
    values[41] = dCa_SLHigh_jct2
    dCa_jct1_tot_bound = dCa_SLHigh_jct1 + dCa_SLB_jct1
    dCa_jct2_tot_bound = dCa_SLB_jct2 + dCa_SLHigh_jct2
    dCa_SL_tot_bound = dCa_SLB_SL + dCa_SLHigh_SL
    i_Ca_jct1_tot = i_CaL_Ca_jct1 - 2*i_NaCa_jct1 + i_Cap_jct1 + i_Cab_jct1
    i_Ca_jct2_tot = i_Cap_jct2 - 2*i_NaCa_jct2 + i_Cab_jct2 + i_CaL_Ca_jct2
    i_Ca_SL_tot = i_Cap_SL + i_CaL_Ca_SL - 2*i_NaCa_SL + i_Cab_SL
    values[42] = -dCalsequestrin - (j_leak_SR1 + j_leak_SR2)*Vol_myo/Vol_SR +\
        j_pump_SR - j_rel_SR1 - j_rel_SR2
    values[43] = -J_Ca_jct1_SL/Vol_jct1 - Cm*i_Ca_jct1_tot/(2*F*Vol_jct1) -\
        dCa_jct1_tot_bound + Vol_myo*j_leak_SR1/Vol_jct1 +\
        Vol_SR*j_rel_SR1/Vol_jct1
    values[44] = Vol_myo*j_leak_SR2/Vol_jct2 -\
        Cm*i_Ca_jct2_tot/(2*F*Vol_jct2) + Vol_SR*j_rel_SR2/Vol_jct2 -\
        dCa_jct2_tot_bound - J_Ca_jct2_SL/Vol_jct2
    values[45] = (J_Ca_jct1_SL - J_Ca_SL_myo + J_Ca_jct2_SL)/Vol_SL -\
        dCa_SL_tot_bound - Cm*i_Ca_SL_tot/(2*F*Vol_SL)
    values[46] = -dCa_cytosol_tot_bound - Vol_SR*j_pump_SR/Vol_myo +\
        J_Ca_SL_myo/Vol_myo

    # Expressions for the Cell component
    i_Stim = (-stim_amplitude if -stim_period*math.floor(t/stim_period) + t\
        <= stim_start + stim_duration and\
        -stim_period*math.floor(t/stim_period) + t >= stim_start else 0)
    values[47] = -i_NaCa - i_Cl_Ca - i_Cab - i_Clb - i_tos - i_CaL - i_tof -\
        i_Stim - i_Na - i_Cap - i_Nab - i_NaK - i_Kr - i_K1 - i_Kp - i_Ks

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the shannon_splitcleft ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 48)
    h, j, m, Xr, Xs, X_tos, Y_tos, R_tos, X_tof, Y_tof, d, f, fCaB_SL,\
        fCaB_jct1, fCaB_jct2, R1, R2, O1, O2, I1, I2, Ca_TroponinC,\
        Ca_TroponinC_Ca_Mg, Mg_TroponinC_Ca_Mg, Ca_Calmodulin, Ca_Myosin,\
        Mg_Myosin, Ca_SRB, Na_jct1_buf, Na_jct2_buf, Na_SL_buf, Na_jct1,\
        Na_jct2, Na_SL, Nai, Ca_Calsequestrin, Ca_SLB_SL, Ca_SLB_jct1,\
        Ca_SLB_jct2, Ca_SLHigh_SL, Ca_SLHigh_jct1, Ca_SLHigh_jct2, Ca_SR,\
        Ca_jct1, Ca_jct2, Ca_SL, Cai, V = states

    # Assign parameters
    assert(len(parameters) == 137)
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length,\
        cell_radius, Fx_Na_SL, Fx_Na_jct1, Fx_Na_jct2, G_INa, Fx_NaBk_SL,\
        Fx_NaBk_jct1, Fx_NaBk_jct2, G_NaBk, Fx_NaK_SL, Fx_NaK_jct1,\
        Fx_NaK_jct2, H_NaK, I_NaK_max, Km_Ko, Km_Nai, Fx_Ks_SL, Fx_Ks_jct1,\
        Fx_Ks_jct2, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1,\
        Fx_Cl_jct2, G_Cl, Kd_ClCa, G_ClBk, Fx_NCX_SL, Fx_NCX_jct1,\
        Fx_NCX_jct2, HNa, K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX,\
        V_max_INaCa, eta, ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, Fx_SLCaP_jct2,\
        H_ICap, Km, Q10_SLCaP, V_maxAF, Fx_CaBk_SL, Fx_CaBk_jct1,\
        Fx_CaBk_jct2, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa,\
        kom, ks1, ks2, KSRleak1, KSRleak2, H, Kmf, Kmr, Q10_SRCaP, V_max,\
        Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg, Bmax_SRB,\
        Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca, Bmax_TroponinC_Ca_Mg_Mg,\
        koff_Calmodulin, koff_Myosin_Ca, koff_Myosin_Mg, koff_SRB,\
        koff_TroponinC, koff_TroponinC_Ca_Mg_Ca, koff_TroponinC_Ca_Mg_Mg,\
        kon_Calmodulin, kon_Myosin_Ca, kon_Myosin_Mg, kon_SRB, kon_TroponinC,\
        kon_TroponinC_Ca_Mg_Ca, kon_TroponinC_Ca_Mg_Mg, Fx_ICaL_SL,\
        Fx_ICaL_jct1, Fx_ICaL_jct2, PCa, PK, PNa, Q10_CaL, gamma_Cai,\
        gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai, gamma_Nao, Bmax_SL,\
        Bmax_jct1, Bmax_jct2, koff, kon, Bmax_Calsequestrin, Bmax_SLB_SL,\
        Bmax_SLB_jct1, Bmax_SLB_jct2, Bmax_SLHigh_SL, Bmax_SLHigh_jct1,\
        Bmax_SLHigh_jct2, koff_Calsequestrin, koff_SLB, koff_SLHigh,\
        kon_Calsequestrin, kon_SL, stim_amplitude, stim_duration,\
        stim_period, stim_start = parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((200,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (200,)

    # Expressions for the Model parameters component
    monitored[0] = 3.141592654e-15*cell_length*(cell_radius*cell_radius)
    monitored[1] = 0.035*monitored[0]
    monitored[2] = 0.02*monitored[0]
    monitored[3] = 0.000539*monitored[0]
    monitored[4] = 0
    monitored[5] = 0.65*monitored[0]

    # Expressions for the Reversal potentials component
    monitored[133] = Rgas*T*math.log(Nao/Na_jct1)/F
    monitored[134] = Rgas*T*math.log(Nao/Na_jct2)/F
    monitored[135] = Rgas*T*math.log(Nao/Na_SL)/F
    monitored[136] = Rgas*T*math.log(Cao/Ca_jct1)/(2*F)
    monitored[137] = Rgas*T*math.log(Cao/Ca_jct2)/(2*F)
    monitored[138] = Rgas*T*math.log(Cao/Ca_SL)/(2*F)
    monitored[139] = Rgas*T*math.log(Ko/Ki)/F
    monitored[140] = Rgas*T*math.log(Cli/Clo)/F

    # Expressions for the INa component
    monitored[6] = (m*m*m)*h*j
    monitored[7] = Fx_Na_jct1*G_INa*(V - monitored[133])*monitored[6]
    monitored[8] = Fx_Na_jct2*G_INa*(-monitored[134] + V)*monitored[6]
    monitored[9] = Fx_Na_SL*G_INa*(V - monitored[135])*monitored[6]
    monitored[10] = monitored[9] + monitored[8] + monitored[7]

    # Expressions for the h gate component
    monitored[11] = (1.04951082543e-06*math.exp(-0.147058823529*V) if V < -40 else\
        0)
    monitored[12] = (3.56*math.exp(0.079*V) + 310000.0*math.exp(0.35*V) if V\
        < -40 else 1.0/(0.13 + 0.0497581410839*math.exp(-0.0900900900901*V)))
    monitored[152] = (1 - h)*monitored[11] - h*monitored[12]

    # Expressions for the j gate component
    monitored[13] = (1.0*(37.78 + V)*(-127140.0*math.exp(0.2444*V) -\
        3.474e-05*math.exp(-0.04391*V))/(1 + 50262745826.0*math.exp(0.311*V))\
        if V < -40 else 0)
    monitored[14] = (0.1212*math.exp(-0.01052*V)/(1 +\
        0.0039608683399*math.exp(-0.1378*V)) if V < -40 else\
        0.3*math.exp(-2.535e-07*V)/(1 + 0.0407622039784*math.exp(-0.1*V)))
    monitored[153] = -j*monitored[14] + (1 - j)*monitored[13]

    # Expressions for the m gate component
    monitored[15] = (15.0816 + 0.32*V)/(1 - 0.0089778037307*math.exp(-0.1*V))
    monitored[16] = 0.08*math.exp(-0.0909090909091*V)
    monitored[154] = -m*monitored[16] + (1 - m)*monitored[15]

    # Expressions for the INab component
    monitored[17] = Fx_NaBk_jct1*G_NaBk*(V - monitored[133])
    monitored[18] = Fx_NaBk_jct2*G_NaBk*(-monitored[134] + V)
    monitored[19] = Fx_NaBk_SL*G_NaBk*(V - monitored[135])
    monitored[20] = monitored[19] + monitored[18] + monitored[17]

    # Expressions for the INaK component
    monitored[21] = -0.142857142857 +\
        0.142857142857*math.exp(0.0148588410104*Nao)
    monitored[22] = 1.0/(1 + 0.0365*math.exp(-F*V/(Rgas*T))*monitored[21] +\
        0.1245*math.exp(-0.1*F*V/(Rgas*T)))
    monitored[23] = Fx_NaK_jct1*I_NaK_max*Ko*monitored[22]/((1 +\
        math.pow(Km_Nai/Na_jct1, H_NaK))*(Km_Ko + Ko))
    monitored[24] = Fx_NaK_jct2*I_NaK_max*Ko*monitored[22]/((1 +\
        math.pow(Km_Nai/Na_jct2, H_NaK))*(Km_Ko + Ko))
    monitored[25] = Fx_NaK_SL*I_NaK_max*Ko*monitored[22]/((1 +\
        math.pow(Km_Nai/Na_SL, H_NaK))*(Km_Ko + Ko))
    monitored[26] = monitored[24] + monitored[23] + monitored[25]

    # Expressions for the Xr gate component
    monitored[29] = 1.0/(1 + 0.00127263380134*math.exp(-0.133333333333*V))
    monitored[30] = 1.0/((0.0061 + 0.00061*V)/(-1 +\
        4.26311451517*math.exp(0.145*V)) + (0.00966 + 0.00138*V)/(1 -\
        0.422739131746*math.exp(-0.123*V)))
    monitored[155] = (monitored[29] - Xr)/monitored[30]

    # Expressions for the Rr gate component
    monitored[31] = 1.0/(1 + 4.36323731689*math.exp(0.0446428571429*V))

    # Expressions for the IKs component
    monitored[32] = 3 - math.log(1.0*Ca_jct1)
    monitored[33] = 3 - math.log(1.0*Ca_jct2)
    monitored[34] = 3 - math.log(1.0*Ca_SL)
    monitored[35] = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*monitored[32]))
    monitored[36] = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*monitored[33]))
    monitored[37] = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*monitored[34]))
    monitored[38] = Rgas*T*math.log((Nao*pKNa + Ko)/(pKNa*Nai + Ki))/F
    monitored[39] = Fx_Ks_jct1*(Xs*Xs)*(V - monitored[38])*monitored[35]
    monitored[40] = Fx_Ks_jct2*(Xs*Xs)*(V - monitored[38])*monitored[36]
    monitored[41] = Fx_Ks_SL*(Xs*Xs)*(V - monitored[38])*monitored[37]
    monitored[42] = monitored[39] + monitored[41] + monitored[40]

    # Expressions for the Xs gate component
    monitored[43] = 1.0/(1 + 1.0939777431*math.exp(-0.059880239521*V))
    monitored[44] = 1.0/((0.002157 + 7.19e-05*V)/(1 -\
        0.0117959385198*math.exp(-0.148*V)) + (0.00393 + 0.000131*V)/(-1 +\
        7.85381970442*math.exp(0.0687*V)))
    monitored[156] = (-Xs + monitored[43])/monitored[44]

    # Expressions for the IKp component
    monitored[45] = g_Kp*(V - monitored[139])/(1 +\
        1786.47556538*math.exp(-0.167224080268*V))

    # Expressions for the Itos component
    monitored[46] = G_tos*(0.5*R_tos + Y_tos)*(V - monitored[139])*X_tos

    # Expressions for the X_gate component
    monitored[47] = 1.0/(1 + 0.818730753078*math.exp(-0.0666666666667*V))
    monitored[48] = 0.5 + 9/(1 + 1.22140275816*math.exp(0.0666666666667*V))
    monitored[157] = (-X_tos + monitored[47])/monitored[48]

    # Expressions for the Y_gate component
    monitored[49] = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    monitored[50] = 30 + 3000/(1 + 403.428793493*math.exp(0.1*V))
    monitored[158] = (monitored[49] - Y_tos)/monitored[50]

    # Expressions for the R_gate component
    monitored[51] = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    monitored[52] = 220 + 2800.0/(1 + 403.428793493*math.exp(0.1*V))
    monitored[159] = (-R_tos + monitored[51])/monitored[52]

    # Expressions for the Itof component
    monitored[53] = G_tof*(V - monitored[139])*X_tof*Y_tof

    # Expressions for the Itof X gate component
    monitored[54] = 1.0/(1 + 0.818730753078*math.exp(-0.0666666666667*V))
    monitored[55] = 1.5 + 3.5*math.exp(-0.00111111111111*(V*V))
    monitored[160] = (-X_tof + monitored[54])/monitored[55]

    # Expressions for the Itof Y gate component
    monitored[56] = 1.0/(1 + 28.5027336438*math.exp(0.1*V))
    monitored[57] = 20 + 20/(1 + 28.5027336438*math.exp(0.1*V))
    monitored[161] = (-Y_tof + monitored[56])/monitored[57]

    # Expressions for the K1 gate component
    monitored[143] = 1.02/(1 +\
        7.35454251046e-07*math.exp(-0.2385*monitored[139] + 0.2385*V))
    monitored[144] = (1.15340563519e-16*math.exp(-0.06175*monitored[139] +\
        0.06175*V) + 0.762624006506*math.exp(-0.08032*monitored[139] +\
        0.08032*V))/(1 + 0.0867722941577*math.exp(0.5143*monitored[139] -\
        0.5143*V))
    monitored[145] = monitored[143]/(monitored[144] + monitored[143])

    # Expressions for the ICl Ca component
    monitored[58] = G_Cl*(V - monitored[140])*(Fx_Cl_SL/(1 + Kd_ClCa/Ca_SL) +\
        Fx_Cl_jct2/(1 + Kd_ClCa/Ca_jct2) + Fx_Cl_jct1/(1 + Kd_ClCa/Ca_jct1))

    # Expressions for the IClb component
    monitored[59] = G_ClBk*(V - monitored[140])

    # Expressions for the d gate component
    monitored[70] = 1.0/(1 + 0.0892185174093*math.exp(-0.166666666667*V))
    monitored[71] = (1 -\
        0.0892185174093*math.exp(-0.166666666667*V))*monitored[70]/(0.5075 +\
        0.035*V)
    monitored[162] = (monitored[70] - d)/monitored[71]

    # Expressions for the f gate component
    monitored[72] = 1.0/(1 + 16964.681259*math.exp(0.277777777778*V)) +\
        0.6/(1 + 12.1824939607*math.exp(-0.05*V))
    monitored[73] = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 +\
        0.0337*V)*(0.48865 + 0.0337*V))))
    monitored[163] = (monitored[72] - f)/monitored[73]

    # Expressions for the FCa gate component
    monitored[74] = 1 - fCaB_SL
    monitored[75] = 1 - fCaB_jct1
    monitored[76] = 1 - fCaB_jct2
    monitored[164] = -0.0119*fCaB_SL + 1.7*(1 - fCaB_SL)*Ca_SL
    monitored[165] = -0.0119*fCaB_jct1 + 1.7*(1 - fCaB_jct1)*Ca_jct1
    monitored[166] = 1.7*(1 - fCaB_jct2)*Ca_jct2 - 0.0119*fCaB_jct2

    # Expressions for the INaCa component
    monitored[77] = (-math.pow(Nao, HNa)*Ca_jct1*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct1,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    monitored[78] = (-math.pow(Nao, HNa)*Ca_jct2*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct2,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    monitored[79] = (-math.pow(Nao, HNa)*Ca_SL*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_SL,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    monitored[80] = math.pow(Q10_NCX, -31.0 + 0.1*T)
    monitored[81] = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_SL*Ca_SL*Ca_SL))
    monitored[82] = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct1*Ca_jct1*Ca_jct1))
    monitored[83] = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct2*Ca_jct2*Ca_jct2))
    monitored[84] =\
        Fx_NCX_jct1*V_max_INaCa*monitored[77]*monitored[80]*monitored[82]/(K_mCai*math.pow(Nao,\
        HNa)*(1 + math.pow(Na_jct1/K_mNai, HNa)) + K_mCao*math.pow(Na_jct1,\
        HNa) + Cao*math.pow(Na_jct1, HNa) + math.pow(Nao, HNa)*Ca_jct1 +\
        math.pow(K_mNao, HNa)*(1 + Ca_jct1/K_mCai)*Ca_jct1)
    monitored[85] =\
        Fx_NCX_jct2*V_max_INaCa*monitored[78]*monitored[80]*monitored[83]/(math.pow(K_mNao,\
        HNa)*(1 + Ca_jct2/K_mCai)*Ca_jct2 + K_mCai*math.pow(Nao, HNa)*(1 +\
        math.pow(Na_jct2/K_mNai, HNa)) + Cao*math.pow(Na_jct2, HNa) +\
        math.pow(Nao, HNa)*Ca_jct2 + K_mCao*math.pow(Na_jct2, HNa))
    monitored[86] =\
        Fx_NCX_SL*V_max_INaCa*monitored[79]*monitored[80]*monitored[81]/(math.pow(K_mNao,\
        HNa)*(1 + Ca_SL/K_mCai)*Ca_SL + math.pow(Nao, HNa)*Ca_SL +\
        Cao*math.pow(Na_SL, HNa) + K_mCao*math.pow(Na_SL, HNa) +\
        K_mCai*math.pow(Nao, HNa)*(1 + math.pow(Na_SL/K_mNai, HNa)))
    monitored[87] = monitored[86] + monitored[85] + monitored[84]

    # Expressions for the ICap component
    monitored[88] = math.pow(Q10_SLCaP, -31.0 + 0.1*T)
    monitored[89] = Fx_SLCaP_jct1*V_maxAF*monitored[88]/(1 +\
        math.pow(Km/Ca_jct1, H_ICap))
    monitored[90] = Fx_SLCaP_jct2*V_maxAF*monitored[88]/(1 +\
        math.pow(Km/Ca_jct2, H_ICap))
    monitored[91] = Fx_SLCaP_SL*V_maxAF*monitored[88]/(1 + math.pow(Km/Ca_SL,\
        H_ICap))
    monitored[92] = monitored[89] + monitored[91] + monitored[90]

    # Expressions for the ICab component
    monitored[93] = Fx_CaBk_jct1*G_CaBk*(V - monitored[136])
    monitored[94] = Fx_CaBk_jct2*G_CaBk*(-monitored[137] + V)
    monitored[95] = Fx_CaBk_SL*G_CaBk*(V - monitored[138])
    monitored[96] = monitored[93] + monitored[94] + monitored[95]

    # Expressions for the Jrel SR component
    monitored[97] = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR,\
        HSR))
    monitored[98] = koCa/monitored[97]
    monitored[99] = kiCa*monitored[97]
    monitored[100] = 1 - I1 - O1 - R1
    monitored[101] = 1 - R2 - O2 - I2
    monitored[167] = kom*O1 + kim*monitored[100] - Ca_jct1*R1*monitored[99] -\
        (Ca_jct1*Ca_jct1)*R1*monitored[98]
    monitored[168] = -(Ca_jct2*Ca_jct2)*R2*monitored[98] -\
        Ca_jct2*R2*monitored[99] + kim*monitored[101] + kom*O2
    monitored[169] = kim*I1 - kom*O1 + (Ca_jct1*Ca_jct1)*R1*monitored[98] -\
        Ca_jct1*O1*monitored[99]
    monitored[170] = (Ca_jct2*Ca_jct2)*R2*monitored[98] + kim*I2 -\
        Ca_jct2*O2*monitored[99] - kom*O2
    monitored[171] = (Ca_jct1*Ca_jct1)*monitored[100]*monitored[98] +\
        Ca_jct1*O1*monitored[99] - kom*I1 - kim*I1
    monitored[172] = (Ca_jct2*Ca_jct2)*monitored[101]*monitored[98] +\
        Ca_jct2*O2*monitored[99] - kim*I2 - kom*I2
    monitored[102] = ks1*(Ca_SR - Ca_jct1)*O1
    monitored[103] = ks2*(Ca_SR - Ca_jct2)*O2

    # Expressions for the Jleak SR component
    monitored[104] = KSRleak1*(Ca_SR - Ca_jct1)
    monitored[105] = KSRleak2*(Ca_SR - Ca_jct2)

    # Expressions for the Jpump SR component
    monitored[106] = math.pow(Q10_SRCaP, -31.0 + 0.1*T)
    monitored[107] = V_max*(math.pow(Cai/Kmf, H) - math.pow(Ca_SR/Kmr,\
        H))*monitored[106]/(1 + math.pow(Cai/Kmf, H) + math.pow(Ca_SR/Kmr,\
        H))

    # Expressions for the Ion diffusion component
    monitored[146] = -1.8313e-14*Na_SL + 1.8313e-14*Na_jct1
    monitored[147] = 0
    monitored[148] = 1.6386e-12*Na_SL - 1.6386e-12*Nai
    monitored[149] = 8.2413e-13*Ca_jct1 - 8.2413e-13*Ca_SL
    monitored[150] = 0
    monitored[151] = -3.7243e-12*Cai + 3.7243e-12*Ca_SL

    # Expressions for the Cytosolic component
    monitored[121] = kon_TroponinC*(-Ca_TroponinC + Bmax_TroponinC)*Cai -\
        koff_TroponinC*Ca_TroponinC
    monitored[122] = -koff_TroponinC_Ca_Mg_Ca*Ca_TroponinC_Ca_Mg +\
        kon_TroponinC_Ca_Mg_Ca*(-Ca_TroponinC_Ca_Mg - Mg_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Ca)*Cai
    monitored[123] = Mgi*kon_TroponinC_Ca_Mg_Mg*(-Ca_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Mg - Mg_TroponinC_Ca_Mg) -\
        koff_TroponinC_Ca_Mg_Mg*Mg_TroponinC_Ca_Mg
    monitored[124] = -koff_Calmodulin*Ca_Calmodulin +\
        kon_Calmodulin*(Bmax_Calmodulin - Ca_Calmodulin)*Cai
    monitored[125] = kon_Myosin_Ca*(Bmax_Myosin_Ca - Mg_Myosin -\
        Ca_Myosin)*Cai - koff_Myosin_Ca*Ca_Myosin
    monitored[126] = Mgi*kon_Myosin_Mg*(Bmax_Myosin_Mg - Mg_Myosin -\
        Ca_Myosin) - koff_Myosin_Mg*Mg_Myosin
    monitored[127] = kon_SRB*(-Ca_SRB + Bmax_SRB)*Cai - koff_SRB*Ca_SRB
    monitored[128] = monitored[124] + monitored[127] + monitored[125] +\
        monitored[126] + monitored[122] + monitored[123] + monitored[121]
    monitored[173] = monitored[121]
    monitored[174] = monitored[122]
    monitored[175] = monitored[123]
    monitored[176] = monitored[124]
    monitored[177] = monitored[125]
    monitored[178] = monitored[126]
    monitored[179] = monitored[127]

    # Expressions for the IKr component
    monitored[27] = 0.0129099444874*math.sqrt(Ko)
    monitored[28] = (V - monitored[139])*Xr*monitored[27]*monitored[31]

    # Expressions for the IK1 component
    monitored[141] = 0.387298334621*math.sqrt(Ko)
    monitored[142] = (V - monitored[139])*monitored[141]*monitored[145]

    # Expressions for the ICaL component
    monitored[60] = math.pow(Q10_CaL, -31.0 + 0.1*T)
    monitored[61] = 0.45*(F*F)*V*d*f*monitored[60]/(Rgas*T)
    monitored[62] =\
        4*Fx_ICaL_jct1*PCa*(gamma_Cai*Ca_jct1*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*monitored[61]*monitored[75]/(-1 +\
        math.exp(2*F*V/(Rgas*T)))
    monitored[63] =\
        4*Fx_ICaL_jct2*PCa*(gamma_Cai*Ca_jct2*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*monitored[61]*monitored[76]/(-1 +\
        math.exp(2*F*V/(Rgas*T)))
    monitored[64] = Fx_ICaL_jct1*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct1*math.exp(F*V/(Rgas*T)))*monitored[61]*monitored[75]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[65] = Fx_ICaL_jct2*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct2*math.exp(F*V/(Rgas*T)))*monitored[61]*monitored[76]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[66] =\
        4*Fx_ICaL_SL*PCa*(gamma_Cai*Ca_SL*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*monitored[61]*monitored[74]/(-1 +\
        math.exp(2*F*V/(Rgas*T)))
    monitored[67] = Fx_ICaL_SL*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_SL*math.exp(F*V/(Rgas*T)))*monitored[61]*monitored[74]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[68] = PK*(-Ko*gamma_Ko +\
        Ki*gamma_Ki*math.exp(F*V/(Rgas*T)))*(Fx_ICaL_SL*monitored[74] +\
        Fx_ICaL_jct1*monitored[75] +\
        Fx_ICaL_jct2*monitored[76])*monitored[61]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[69] = monitored[66] + monitored[63] + monitored[67] +\
        monitored[65] + monitored[62] + monitored[64] + monitored[68]

    # Expressions for the Na buffer component
    monitored[129] = kon*(Bmax_jct1 - Na_jct1_buf)*Na_jct1 - koff*Na_jct1_buf
    monitored[130] = kon*(Bmax_jct2 - Na_jct2_buf)*Na_jct2 - koff*Na_jct2_buf
    monitored[131] = -koff*Na_SL_buf + kon*(-Na_SL_buf + Bmax_SL)*Na_SL
    monitored[180] = monitored[129]
    monitored[181] = monitored[130]
    monitored[182] = monitored[131]
    monitored[183] = -Cm*(3*monitored[23] + monitored[17] + monitored[7] +\
        3*monitored[84] + monitored[64])/(F*monitored[3]) -\
        monitored[146]/monitored[3] - monitored[129]
    monitored[184] = -monitored[147]/monitored[4] - Cm*(monitored[18] +\
        3*monitored[24] + 3*monitored[85] + monitored[65] +\
        monitored[8])/(F*monitored[4]) - monitored[130]
    monitored[185] = -Cm*(monitored[9] + monitored[67] + monitored[19] +\
        3*monitored[86] + 3*monitored[25])/(F*monitored[2]) - monitored[131]\
        + (monitored[146] - monitored[148] + monitored[147])/monitored[2]
    monitored[186] = monitored[148]/monitored[5]

    # Expressions for the Ca buffer component
    monitored[108] = kon_Calsequestrin*(-Ca_Calsequestrin +\
        Bmax_Calsequestrin*monitored[5]/monitored[1])*Ca_SR -\
        koff_Calsequestrin*Ca_Calsequestrin
    monitored[187] = monitored[108]
    monitored[109] = -koff_SLB*Ca_SLB_SL + kon_SL*(-Ca_SLB_SL +\
        Bmax_SLB_SL*monitored[5]/monitored[2])*Ca_SL
    monitored[110] = kon_SL*(0.1*Bmax_SLB_jct1*monitored[5]/monitored[3] -\
        Ca_SLB_jct1)*Ca_jct1 - koff_SLB*Ca_SLB_jct1
    monitored[111] = -koff_SLB*Ca_SLB_jct2 + kon_SL*(-Ca_SLB_jct2 +\
        0.1*Bmax_SLB_jct2*monitored[5]/monitored[4])*Ca_jct2
    monitored[112] = -koff_SLHigh*Ca_SLHigh_SL +\
        kon_SL*(Bmax_SLHigh_SL*monitored[5]/monitored[2] -\
        Ca_SLHigh_SL)*Ca_SL
    monitored[113] = kon_SL*(0.1*Bmax_SLHigh_jct1*monitored[5]/monitored[3] -\
        Ca_SLHigh_jct1)*Ca_jct1 - koff_SLHigh*Ca_SLHigh_jct1
    monitored[114] = kon_SL*(0.1*Bmax_SLHigh_jct2*monitored[5]/monitored[4] -\
        Ca_SLHigh_jct2)*Ca_jct2 - koff_SLHigh*Ca_SLHigh_jct2
    monitored[188] = monitored[109]
    monitored[189] = monitored[110]
    monitored[190] = monitored[111]
    monitored[191] = monitored[112]
    monitored[192] = monitored[113]
    monitored[193] = monitored[114]
    monitored[115] = monitored[113] + monitored[110]
    monitored[116] = monitored[111] + monitored[114]
    monitored[117] = monitored[109] + monitored[112]
    monitored[118] = monitored[93] - 2*monitored[84] + monitored[89] +\
        monitored[62]
    monitored[119] = monitored[63] + monitored[94] + monitored[90] -\
        2*monitored[85]
    monitored[120] = monitored[66] - 2*monitored[86] + monitored[91] +\
        monitored[95]
    monitored[194] = -(monitored[105] +\
        monitored[104])*monitored[5]/monitored[1] - monitored[103] +\
        monitored[107] - monitored[108] - monitored[102]
    monitored[195] = monitored[104]*monitored[5]/monitored[3] +\
        monitored[102]*monitored[1]/monitored[3] - monitored[115] -\
        monitored[149]/monitored[3] - Cm*monitored[118]/(2*F*monitored[3])
    monitored[196] = -monitored[150]/monitored[4] -\
        Cm*monitored[119]/(2*F*monitored[4]) +\
        monitored[105]*monitored[5]/monitored[4] - monitored[116] +\
        monitored[103]*monitored[1]/monitored[4]
    monitored[197] = (-monitored[151] + monitored[149] +\
        monitored[150])/monitored[2] - monitored[117] -\
        Cm*monitored[120]/(2*F*monitored[2])
    monitored[198] = monitored[151]/monitored[5] - monitored[128] -\
        monitored[107]*monitored[1]/monitored[5]

    # Expressions for the Cell component
    monitored[132] = (-stim_amplitude if\
        -stim_period*math.floor(t/stim_period) + t <= stim_start +\
        stim_duration and -stim_period*math.floor(t/stim_period) + t >=\
        stim_start else 0)
    monitored[199] = -monitored[87] - monitored[53] - monitored[69] -\
        monitored[58] - monitored[96] - monitored[59] - monitored[92] -\
        monitored[142] - monitored[10] - monitored[20] - monitored[42] -\
        monitored[26] - monitored[132] - monitored[46] - monitored[28] -\
        monitored[45]

    # Return results
    return monitored
