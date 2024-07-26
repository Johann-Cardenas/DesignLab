# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:10:46 2022

@author: aravind4
"""

import numpy as np
from copy import deepcopy

def Nf_cracking(Structure, Single, Tandem, InputMat,CT): 
    if Structure == ['NSS','NTS']:
        Eo = 32592.26 * 145.038     # Instantaneous Modulus converted to psi
        H = 4.0 # AC thickness in inches    
        gi = np.array([0.0583, 0.0572, 0.0970, 0.1281, 0.1569, 0.1697, 0.1428, 0.0927, 0.0544, 0.0187, 0.0155])
        tt = np.array([1.00e-05, 1.00e-04, 0.001, 0.01, 0.1, 1, 10., 100., 1000., 10000., 100000.])
    if Structure == ['NSW','NTW']:
        Eo = 35668.37 * 145.038     # Instantaneous Modulus converted to psi
        H = 4.0 # AC thickness in inches    
        gi = np.array([0.2768, 0.0534, 0.1806, 0.1386, 0.1309, 0.0938, 0.0601, 0.0331, 0.0174, 0.0076, 0.0042])
        tt = np.array([1.00e-05, 1.00e-04, 0.001, 0.01, 0.1, 1, 10., 100., 1000., 10000., 100000.])
    if Structure == ['KSS','KTS']:
        Eo = 29602.67 * 145.038     # Instantaneous Modulus converted to psi
        H = 11.0 # AC thickness in inches    
        gi = np.array([0.1564, 0.0388, 0.1499, 0.1390, 0.1639, 0.1468, 0.1019, 0.0551, 0.0274, 0.0104, 0.0054])
        tt = np.array([1.00e-05, 1.00e-04, 0.001, 0.01, 0.1, 1, 10., 100., 1000., 10000., 100000.])
    if Structure == ['KSW','KTW']:
        Eo = 40317.79 * 145.038     # Instantaneous Modulus converted to psi
        H = 11.0 # AC thickness in inches    
        gi = np.array([0.0035, 0.1003, 0.0476, 0.2531, 0.0000, 0.1202, 0.1167, 0.1108, 0.0850, 0.0651, 0.0437, 0.0259, 0.0139, 0.0073, 0.0040])
        tt = np.array([1.00e-10, 1.00e-09, 1.00e-08, 1.00e-07, 1.00e-06, 1.00e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10., 100., 1000., 10000.])
    
    k_f1, k_f2, k_f3 = 0.007566, 3.9492, 1.281
    bb_f1, bb_f2, bb_f3 = 1.0, 1.0, 1.0
    V_be, V_a = 12.36, 3.34
    M = 4.84 * (V_be / (V_a + V_be) - 0.69)    
    C = 10**M
    E = 35000
    if CT == 'BU':
        resp = '_BotAC_Z'
        C_h = 1 / (0.000398 + 0.003602 / (1 + np.exp(11.02-3.49*H)))
    elif CT == 'TD':
        resp ='_Surf_Z'
        C_h = 1 / (0.01 + 12.00 / (1 + np.exp(15.676 - 2.8186*H)))
    E11_Single = deepcopy(Single.loc[:,'E11'+resp].values)
    E11_Tandem = deepcopy(Tandem.loc[:,'E11'+resp].values)
    E33_Single = deepcopy(Single.loc[:,'E33'+resp].values)
    E33_Tandem = deepcopy(Tandem.loc[:,'E33'+resp].values)
    E11_Single[E11_Single < 1e-8] = 1e-8
    E11_Tandem[E11_Tandem < 1e-8] = 1e-8
    E33_Single[E33_Single < 1e-8] = 1e-8
    E33_Tandem[E33_Tandem < 1e-8] = 1e-8        
    ee_s = np.array([max(ii, jj) for ii, jj in zip(E11_Single, E33_Single)])
    ee_t = np.array([max(ii, jj) for ii, jj in zip(E11_Tandem, E33_Tandem)])
    
    Nf_s = k_f1 * C * C_h * bb_f1 * (ee_s)**(-k_f2*bb_f2) * (E)**(-k_f3*bb_f3)
    Nf_t = k_f1 * C * C_h * bb_f1 * (ee_t)**(-k_f2*bb_f2) * (E)**(-k_f3*bb_f3)
    return Nf_s, Nf_t

def FC_BU_(DI_bottom, struct):
    if struct == ['KSS','KTS'] or struct == ['KSW','KTW']:
        H  = 11
    elif struct== ['NSS','NTS'] or struct ==['NSW','NTW']:
        H = 4
    C1, C2, C4 = 1, 1, 6000
    C2s = -2.40874 - 39.748 * (1 + H)**(-2.856)
    C1s = -2 * C2s
    return 1/60 * C4 / (1 + np.exp(C1*C1s + C2*C2s*np.log10(DI_bottom * 100)))

def FC_TD_(DI_tdc):
    C1, C2, C4 = 7.00, 3.50, 1000
    return 10.56 * (C4 / (1 + np.exp(C1 - C2*np.log10(DI_tdc * 100))))