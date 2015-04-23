# Gotran generated code for the  "simple" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # R=0.8884332, O=8.156628e-07, I=1.024274e-07
    init_values = np.array([0.8884332, 8.156628e-07, 1.024274e-07],\
        dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("R",(0, Range())), ("O",(1, Range())), ("I",(2,\
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
    # EC50_SR=0.45, HSR=2.5, Max_SR=15, Min_SR=1, kiCa=0.5, kim=0.005,
    # koCa=10, kom=0.06, ks=25, Ca_SR=0.5545201, Ca_jct=0.0001737475
    init_values = np.array([0.45, 2.5, 15, 1, 0.5, 0.005, 10, 0.06, 25,\
        0.5545201, 0.0001737475], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("EC50_SR", (0, Range())), ("HSR", (1, Range())),\
        ("Max_SR", (2, Range())), ("Min_SR", (3, Range())), ("kiCa", (4,\
        Range())), ("kim", (5, Range())), ("koCa", (6, Range())), ("kom", (7,\
        Range())), ("ks", (8, Range())), ("Ca_SR", (9, Range())), ("Ca_jct",\
        (10, Range()))])

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
    state_inds = dict([("R", 0), ("O", 1), ("I", 2)])

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
    param_inds = dict([("EC50_SR", 0), ("HSR", 1), ("Max_SR", 2), ("Min_SR",\
        3), ("kiCa", 4), ("kim", 5), ("koCa", 6), ("kom", 7), ("ks", 8),\
        ("Ca_SR", 9), ("Ca_jct", 10)])

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
    monitor_inds = dict([("kCaSR", 0), ("koSRCa", 1), ("kiSRCa", 2), ("RI",\
        3), ("j_rel_SR", 4), ("dR_dt", 5), ("dO_dt", 6), ("dI_dt", 7)])

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
    Compute the right hand side of the simple ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 3)
    R, O, I = states

    # Assign parameters
    assert(len(parameters) == 11)
    EC50_SR=parameters[0]; HSR=parameters[1]; Max_SR=parameters[2];\
        Min_SR=parameters[3]; kiCa=parameters[4]; kim=parameters[5];\
        koCa=parameters[6]; kom=parameters[7]; Ca_SR=parameters[9];\
        Ca_jct=parameters[10]

    # Init return args
    if values is None:
        values = np.zeros((3,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (3,)

    # Expressions for the Jrel sr component
    kCaSR = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR, HSR))
    koSRCa = koCa/kCaSR
    kiSRCa = kiCa*kCaSR
    RI = 1 - I - O - R

    # The ODE system: 3 states
    values[0] = -(Ca_jct*Ca_jct)*R*koSRCa + kom*O + kim*RI - Ca_jct*R*kiSRCa
    values[1] = kim*I + (Ca_jct*Ca_jct)*R*koSRCa - Ca_jct*O*kiSRCa - kom*O
    values[2] = -kim*I + Ca_jct*O*kiSRCa - kom*I + (Ca_jct*Ca_jct)*RI*koSRCa

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the simple ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 3)
    R, O, I = states

    # Assign parameters
    assert(len(parameters) == 11)
    EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa, kom, ks, Ca_SR, Ca_jct =\
        parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((8,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (8,)

    # Expressions for the Jrel sr component
    monitored[0] = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR,\
        HSR))
    monitored[1] = koCa/monitored[0]
    monitored[2] = kiCa*monitored[0]
    monitored[3] = 1 - I - O - R
    monitored[4] = ks*(-Ca_jct + Ca_SR)*O

    # The ODE system: 3 states
    monitored[5] = kom*O - (Ca_jct*Ca_jct)*R*monitored[1] -\
        Ca_jct*R*monitored[2] + kim*monitored[3]
    monitored[6] = (Ca_jct*Ca_jct)*R*monitored[1] + kim*I -\
        Ca_jct*O*monitored[2] - kom*O
    monitored[7] = (Ca_jct*Ca_jct)*monitored[1]*monitored[3] - kim*I - kom*I\
        + Ca_jct*O*monitored[2]

    # Return results
    return monitored
