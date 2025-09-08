import numpy as np

def solve_Q_dot_conduction(delta_T, h_tank, k_w, A_xsec):

    L_w_cond = 0.5*h_tank
    Q_dot_conduction = k_w *(delta_T)*A_xsec/L_w_cond
    
    return Q_dot_conduction