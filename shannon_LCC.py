# Gotran generated code for: shannon_LCC

from __future__ import division

def init_values(**values):
    """
    Init values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # d=7.175662e-06, f=1.000681, fCaB_SL=0.01452605, fCaB_jct=0.02421991
    init_values = np.array([7.175662e-06, 1.000681, 0.01452605, 0.02421991],\
        dtype=np.float_)

    # State indices and limit checker
    state_ind = dict(d=(0, Range()), f=(1, Range()), fCaB_SL=(2, Range()),\
        fCaB_jct=(3, Range()))

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
    # Ca_SL=0.0001031812, Ca_jct=0.0001737475, Na_SL=8.80733, Na_jct=8.80329,
    # V=-85.0, dV=14.5, Cao=1.8, Cli=15, Clo=150, Cm=1.381e-10,
    # F=96485, Ki=135, Ko=5.4, Mgi=1, Nao=140, Rgas=8314.3, T=310,
    # cell_length=100, cell_radius=10.25, Fx_ICaL_SL=0.1,
    # Fx_ICaL_jct=0.9, PCa=0.00054, PK=2.7e-07, PNa=1.5e-08,
    # Q10_CaL=1.8, gamma_Cai=0.341, gamma_Cao=0.341, gamma_Ki=0.75,
    # gamma_Ko=0.75, gamma_Nai=0.75, gamma_Nao=0.75
    param_values = np.array([0.0001031812, 0.0001737475, 8.80733, 8.80329,\
        -85.0, 14.5, 1.8, 15, 150, 1.381e-10, 96485, 135, 5.4, 1, 140,\
        8314.3, 310, 100, 10.25, 0.1, 0.9, 0.00054, 2.7e-07, 1.5e-08, 1.8,\
        0.341, 0.341, 0.75, 0.75, 0.75, 0.75], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict(Ca_SL=(0, Range()), Ca_jct=(1, Range()), Na_SL=(2,\
        Range()), Na_jct=(3, Range()), V=(4, Range()), dV=(5, Range()),\
        Cao=(6, Range()), Cli=(7, Range()), Clo=(8, Range()), Cm=(9,\
        Range()), F=(10, Range()), Ki=(11, Range()), Ko=(12, Range()),\
        Mgi=(13, Range()), Nao=(14, Range()), Rgas=(15, Range()), T=(16,\
        Range()), cell_length=(17, Range()), cell_radius=(18, Range()),\
        Fx_ICaL_SL=(19, Range()), Fx_ICaL_jct=(20, Range()), PCa=(21,\
        Range()), PK=(22, Range()), PNa=(23, Range()), Q10_CaL=(24, Range()),\
        gamma_Cai=(25, Range()), gamma_Cao=(26, Range()), gamma_Ki=(27,\
        Range()), gamma_Ko=(28, Range()), gamma_Nai=(29, Range()),\
        gamma_Nao=(30, Range()))

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
    assert(len(states) == 4)
    d, f, fCaB_SL, fCaB_jct = states

    # Assign parameters
    assert(len(parameters) == 31)
    Ca_SL, Ca_jct, Na_SL, Na_jct, V, dV, Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi,\
        Nao, Rgas, T, cell_length, cell_radius, Fx_ICaL_SL, Fx_ICaL_jct, PCa,\
        PK, PNa, Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki, gamma_Ko,\
        gamma_Nai, gamma_Nao = parameters

    # Ical d gate
    d_infinity = 1.0/(1.0 + math.exp(-dV/6.0 - V/6.0))
    tau_d = (1.0 - math.exp(-dV/6.0 - V/6.0))*d_infinity/(0.035*dV + 0.035*V)

    # Ical f gate
    f_infinity = 0.6/(1.0 + math.exp(5/2 - V/20.0)) + 1.0/(1.0 +\
        16964.6812589854*math.exp(0.277777777777778*V))
    tau_f = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 + 0.0337*V)*(0.48865 +\
        0.0337*V))))

    # Ical fca gate
    fCa_SL = 1.0 - fCaB_SL
    fCa_jct = 1.0 - fCaB_jct

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
    i_CaL_Na_SL = (Na_SL*gamma_Nai*math.exp(F*V/(Rgas*T)) -\
        Nao*gamma_Nao)*Fx_ICaL_SL*PNa*fCa_SL*temp/(-1.0 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_K = (Fx_ICaL_jct*fCa_jct +\
        Fx_ICaL_SL*fCa_SL)*(Ki*gamma_Ki*math.exp(F*V/(Rgas*T)) -\
        Ko*gamma_Ko)*PK*temp/(-1.0 + math.exp(F*V/(Rgas*T)))
    i_CaL = i_CaL_Na_jct + i_CaL_Na_SL + i_CaL_Ca_jct + i_CaL_K + i_CaL_Ca_SL

    # The ODE system: 39 states

    # Init dy
    if dy is None:
        dy = np.zeros_like(states)
    dy[0] = (-d + d_infinity)/tau_d
    dy[1] = (f_infinity - f)/tau_f
    dy[2] = -0.0119*fCaB_SL + 1.7*(1.0 - fCaB_SL)*Ca_SL
    dy[3] = 1.7*(1.0 - fCaB_jct)*Ca_jct - 0.0119*fCaB_jct

    # Return dy
    return dy

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict(d=0, f=1, fCaB_SL=2, fCaB_jct=3)

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
    param_inds = dict(Ca_SL=0, Ca_jct=1, Na_SL=2, Na_jct=3, V=4, dV=5, Cao=6,\
        Cli=7, Clo=8, Cm=9, F=10, Ki=11, Ko=12, Mgi=13, Nao=14, Rgas=15,\
        T=16, cell_length=17, cell_radius=18, Fx_ICaL_SL=19, Fx_ICaL_jct=20,\
        PCa=21, PK=22, PNa=23, Q10_CaL=24, gamma_Cai=25, gamma_Cao=26,\
        gamma_Ki=27, gamma_Ko=28, gamma_Nai=29, gamma_Nao=30)

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    return indices if len(indices)>1 else indices[0]

def monitor(states, time, parameters, monitored=None):
    """
    Compute monitored intermediates
    """
    # Imports
    import numpy as np
    import math
    from math import pow, sqrt, log

    # Assign states
    assert(len(states) == 4)
    d, f, fCaB_SL, fCaB_jct = states[0], states[1], states[2], states[3]

    # Assign parameters
    assert(len(parameters) == 31)
    Ca_SL, Ca_jct, Na_SL, Na_jct, V, Cao, F, Ki, Ko, Nao, Rgas, T,\
        Fx_ICaL_SL, Fx_ICaL_jct, PCa, PK, PNa, Q10_CaL, gamma_Cai, gamma_Cao,\
        gamma_Ki, gamma_Ko, gamma_Nai, gamma_Nao = parameters[0],\
        parameters[1], parameters[2], parameters[3], parameters[4],\
        parameters[6], parameters[10], parameters[11], parameters[12],\
        parameters[14], parameters[15], parameters[16], parameters[19],\
        parameters[20], parameters[21], parameters[22], parameters[23],\
        parameters[24], parameters[25], parameters[26], parameters[27],\
        parameters[28], parameters[29], parameters[30]

    # Common Sub Expressions for monitored intermediates
    cse_monitored_0 = (F*F)
    cse_monitored_1 = math.pow(Q10_CaL, -31.0 + T/10.0)
    cse_monitored_2 = math.exp(F*V/(Rgas*T))
    cse_monitored_3 = (cse_monitored_2*cse_monitored_2)
    cse_monitored_4 = -1.0 + cse_monitored_2
    cse_monitored_5 = -Cao*gamma_Cao
    cse_monitored_6 = -Nao*gamma_Nao
    cse_monitored_7 = (1.0 - fCaB_SL)*Fx_ICaL_SL
    cse_monitored_8 = (1.0 - fCaB_jct)*Fx_ICaL_jct
    cse_monitored_9 = cse_monitored_2*gamma_Nai
    cse_monitored_10 = cse_monitored_3*gamma_Cai
    cse_monitored_11 =\
        0.45*cse_monitored_0*cse_monitored_1*PNa*V*d*f/(cse_monitored_4*Rgas*T)
    cse_monitored_12 = 1.8*cse_monitored_0*cse_monitored_1*PCa*V*d*f/((-1.0 +\
        cse_monitored_3)*Rgas*T)

    # Init monitored
    if monitored is None:
        monitored = np.zeros(1, dtype=np.float_)

    # Monitored intermediates
    monitored[0] = cse_monitored_11*cse_monitored_8*(cse_monitored_6 +\
        cse_monitored_9*Na_jct) +\
        cse_monitored_12*cse_monitored_7*(cse_monitored_10*Ca_SL +\
        cse_monitored_5) +\
        cse_monitored_11*cse_monitored_7*(cse_monitored_9*Na_SL +\
        cse_monitored_6) +\
        0.45*cse_monitored_0*cse_monitored_1*(cse_monitored_7 +\
        cse_monitored_8)*(cse_monitored_2*Ki*gamma_Ki -\
        Ko*gamma_Ko)*PK*V*d*f/(cse_monitored_4*Rgas*T) +\
        cse_monitored_12*cse_monitored_8*(cse_monitored_5 +\
        cse_monitored_10*Ca_jct)

    # Return monitored
    return monitored

