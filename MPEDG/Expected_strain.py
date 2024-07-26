# -*- coding: utf-8 -*-
"""
Created on Thu May 26 08:38:17 2022

@author: aravind4
"""


import numpy as np
import pandas as pd
from copy import deepcopy

def Expected(freq_npl,bin_npl,freq_pl,bin_pl,PS,Single,Tandem,RPFactor):
    FER = deepcopy(Single) #Final Expected Response
    S = deepcopy(Single.iloc[:,1:].values)
    T = deepcopy(Tandem.iloc[:,1:].values)
    ER_npl = np.zeros_like(S,dtype=float)
    ER_pl = np.zeros_like(S,dtype=float)
    ER = np.zeros_like(S,dtype=float)
    # =============================================================================
    #     Non Platooning Trucks
    # =============================================================================
    for resp in range(len(S[0])):
        temp = np.zeros((2*len(S)-1,))
        for f in range(len(freq_npl)):
            temp[f:f+len(S)] = (temp[f:f+len(S)] + 1*freq_npl[f]*S[:,resp] +
                                2*2*freq_npl[f]*T[:,resp])
        ER_npl[:,resp] = temp[int(len(S)/2):int(len(S)/2)+len(S)]
    # =============================================================================
    #     Platooning Trucks
    # =============================================================================
    for resp in range(len(S[0])):
        temp = np.zeros((2*len(S)-1,))
        for pos in range(len(freq_pl)):
            locs = (2*(np.linspace(bin_pl[pos][0], bin_pl[pos][1],5)+24.5))
            for f in range(len(locs)):
                temp[int(locs[f]):int(locs[f])+len(S)] = (temp[int(locs[f]):int(locs[f])+len(S)] + 
                                                1*((PS-1)*RPFactor+1)*freq_pl[pos][f]*S[:,resp] +
                                    ((PS-1)*RPFactor+1)*2*2*freq_pl[pos][f]*T[:,resp])
        ER_pl[:,resp] = temp[int(len(S)/2):int(len(S)/2)+len(S)]
        ER = (ER_npl + ER_pl)/5/(np.sum(freq_npl)+PS*np.sum(freq_pl))
        FER.iloc[:,1:] = ER
    return FER
    
def Expected_Nf(freq_npl,bin_npl,freq_pl,bin_pl,PS,Single,Tandem):
    
    S = np.array(Single).T
    T = np.array(Tandem).T
    ER_S = np.zeros_like(S,dtype=float)
    ER_T = np.zeros_like(T,dtype=float)
    # =============================================================================
    #     Non Platooning Trucks
    # =============================================================================
    for resp in range(len(S[0])):
        temp1 = np.zeros((2*len(S)-1,))
        temp2 = np.zeros((2*len(T)-1,))
        for f in range(len(freq_npl)):
            temp1[f:f+len(S)] = (temp1[f:f+len(S)] + 1*freq_npl[f]*S[:,resp])
            temp2[f:f+len(S)] = (temp2[f:f+len(S)] + 2*2*freq_npl[f]*T[:,resp])
        ER_S[:,resp] = temp1[int(len(S)/2):int(len(S)/2)+len(S)]
        ER_T[:,resp] = temp2[int(len(S)/2):int(len(S)/2)+len(S)]
    # =============================================================================
    #     Platooning Trucks
    # =============================================================================
    for resp in range(len(S[0])):
        temp1 = np.zeros((2*len(S)-1,))
        temp2 = np.zeros((2*len(T)-1,))
        for pos in range(len(freq_pl)):
            locs = (2*(np.linspace(bin_pl[pos][0], bin_pl[pos][1],5)+24.5))
            for f in range(len(locs)):
                temp1[int(locs[f]):int(locs[f])+len(S)] = (temp1[int(locs[f]):int(locs[f])+len(S)] + 
                                                1*(PS)*freq_pl[pos][f]*S[:,resp])
                temp2[int(locs[f]):int(locs[f])+len(S)] = (temp2[int(locs[f]):int(locs[f])+len(S)] +
                                    (PS)*2*2*freq_pl[pos][f]*T[:,resp])
        ER_S[:,resp] = ER_S[:,resp] + temp1[int(len(S)/2):int(len(S)/2)+len(S)]
        ER_T[:,resp] = ER_T[:,resp] + temp2[int(len(S)/2):int(len(S)/2)+len(S)]
    DI_BU = (ER_S[:,0] + 4*ER_T[:,0])
    DI_TD = (ER_S[:,1] + 4*ER_T[:,1])
    return DI_BU, DI_TD   
    

if __name__ == "__main__":
    Resp = pd.ExcelFile('FE_INPUTS - Final.xlsx').parse(['KSS','KTS'])
    Single = Resp['KSS']
    Tandem = Resp['KTS']
    Npl_Trucks = 18000
    Pl_groups = 4000
    PS = 3
    s = np.random.normal(0,0.001,Npl_Trucks)
    b = np.histogram(s,bins=np.linspace(-24.75, 24.75,100))
    freq_npl = b[0]
    # difference_npl = Npl_Trucks-np.sum(freq_npl)
    # if difference_npl>0:
    #     freq_npl[0] = freq_npl[0] + int(difference_npl/2)
    #     freq_npl[-1] = freq_npl[-1] + int(difference_npl/2)
    bin_npl = np.linspace(-24.5,24.5,99) 
    bin_pl = np.array([[-15,-13],[-1,1],[13,15]])
    prob_pl = np.array([1/3,1/3,1/3])
    freq_pl = []
    for i in range(len(bin_pl)):
        s = np.random.uniform(bin_pl[i][0],bin_pl[i][1],int(prob_pl[i]*Pl_groups))
        b = np.histogram(s,bins=np.linspace(bin_pl[i][0], bin_pl[i][1],6))
        freq_pl.append(b[0])
    freq_pl = np.array(freq_pl)
    FER = Expected(freq_npl,bin_npl,freq_pl,bin_pl,PS,Single,Tandem,2)
    
    
    
    
    
    
    
    
    