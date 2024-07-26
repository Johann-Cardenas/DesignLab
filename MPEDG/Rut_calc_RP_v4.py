# -*- coding: utf-8 -*-
"""
Created on Sat May 28 17:50:33 2022

@author: aravind4
"""
import numpy as np
import pandas as pd
from copy import deepcopy
import math
import Temperature_Parameters as TP
import Expected_strain as ES
###############################################################################
########################### ASPHALT RUTTING ###################################
###############################################################################
def RutAccum(ev,hac,depth,n,Rutprev,T,ii):#equivalent number of cycles calcultions
    totthick=sum(hac)
    k1=-3.3512
    k2=0.4791
    k3=1.5606
    C1=-0.1039*totthick**2 +2.4868*totthick -17.342
    C2=0.0172*totthick**2 -1.7331*totthick +27.428
    kz=(C1+C2*depth)*0.328196**depth
    if kz<0:
        kz=0
    cons=hac[ii]*kz*(10**k1)*(T**k3)
    Rut=np.power(np.power(Rutprev,1/k2)+n*np.power(cons*ev,1/k2),k2)
    return Rut

def ACrutting(struct, temperature, FER, Single, Tandem, InputMat, RPFactor, rutAC, rutAC_ref
              ,IP_npl,IP_pl):
    if struct == ['KSS','KTS'] or struct == ['KSW','KTW']:
        sublayerAC = [0.46, 0.48, 1.06, 0.84, 0.99, 4.17, 3]
    elif struct== ['NSS','NTS'] or struct ==['NSW','NTW']:
        sublayerAC = [0.55, 0.59, 0.97, 1.09, 0.80]
    depth = 0
    for ii in range(len(sublayerAC)):
        response = 'E22_AC_'+str(ii+1)+'_Z'
        E22AC_S = -deepcopy(Single.loc[:,response].values)
        E22AC_T = -deepcopy(Tandem.loc[:,response].values)
        E22AC_ref = -deepcopy(FER.loc[:,response].values)
        E22AC_S[E22AC_S< 1e-6] = 1e-6
        E22AC_T[E22AC_T< 1e-6] = 1e-6
        E22AC_ref[E22AC_ref< 1e-6] = 1e-6
        maxloc = np.argmax(E22AC_ref)
        #print(maxloc)
        if ii == 0:
            depth = depth + sublayerAC[ii]/2
        else:
            depth = depth + sublayerAC[ii-1]/2+sublayerAC[ii]/2
        Temp = TP.TProfile(depth, temperature)
        RutprevAC = rutAC[ii,:]
        RutprevAC_ref = rutAC_ref[ii,:]
        TempRut = RutprevAC_ref
        RDAC = np.zeros_like(RutprevAC)
        NPL = InputMat[InputMat == 1]
        PL = InputMat[InputMat != 1]
        
        if len(NPL)>0:
            Coeff = np.zeros((len(InputMat),))
            for i in range(len(NPL)):
                temploc = int(maxloc-2*IP_npl[i])
                if temploc>len(E22AC_ref)-1:
                    temploc = len(E22AC_ref)-1
                elif temploc<0:
                    temploc = 0
                RDAC = RutAccum(E22AC_S[temploc], sublayerAC, depth, 1, RutprevAC[maxloc], Temp,ii)
                RDAC = RutAccum(E22AC_T[temploc], sublayerAC, depth, 4, RDAC, Temp,ii)
                RutprevAC[maxloc] = RDAC
                RutprevAC_ref[maxloc] = RutAccum(max(E22AC_ref), sublayerAC, depth, 5,
                                             RutprevAC_ref[maxloc],Temp,ii)
                Coeff[i] = (RutprevAC[maxloc])/(RutprevAC_ref[maxloc])
                if math.isnan(Coeff[i]):
                    Coeff[i] = 0
                if i>= 400 and abs(np.mean(Coeff[i-50:i])-np.mean(Coeff[i-25:i]))<1e-3: #CONVERGENCE CRITERIA
                    break
            RutprevAC_ref = RutAccum((E22AC_ref), sublayerAC, depth, 5*len(NPL), TempRut,Temp,ii)
            RutprevAC = RutprevAC_ref * Coeff[i]
            TempRut  = RutprevAC_ref
        if len(PL) >0:
            Coeff = np.zeros((len(InputMat),))
            for i in range(len(PL)):
                temploc = int(maxloc-2*IP_pl[i])
                if temploc>len(E22AC_ref)-1:
                    temploc = len(E22AC_ref)-1
                elif temploc<0:
                    temploc = 0
                for j in range(int(PL[i])):
                    if j==0:
                        RDAC = RutAccum(E22AC_S[temploc], sublayerAC, depth, 1, 
                                        RutprevAC[maxloc], Temp,ii)
                        RDAC = RutAccum(E22AC_T[temploc], sublayerAC, depth, 4, RDAC, Temp,ii)
                        RutprevAC[maxloc] = RDAC
                    else:
                        RDAC = RutAccum(E22AC_S[temploc], sublayerAC, depth, 1*RPFactor, 
                                        RutprevAC[maxloc], Temp,ii)
                        RDAC = RutAccum(E22AC_T[temploc], sublayerAC, depth, 4*RPFactor, RDAC, Temp,ii)
                        RutprevAC[maxloc] = RDAC
                RutprevAC_ref[maxloc] = RutAccum(max(E22AC_ref), sublayerAC, depth,  1*(RPFactor*(PL[0]-1)+1)
                                         +4*(RPFactor*(PL[0]-1)+1),
                                             RutprevAC_ref[maxloc],Temp,ii)
                Coeff[i] = (RutprevAC[maxloc])/(RutprevAC_ref[maxloc])
                if math.isnan(Coeff[i]):
                    Coeff[i] = 0
                if i>= 400 and abs(np.mean(Coeff[i-50:i])-np.mean(Coeff[i-25:i]))<1e-3: #CONVERGENCE CRITERIA
                    break
            RutprevAC_ref = RutAccum(E22AC_ref, sublayerAC, depth, 1*(RPFactor*(PL[0]-1)+1)*len(PL)+
                                     4*(RPFactor*(PL[0]-1)+1)*len(PL), 
                                     TempRut,Temp,ii)
            RutprevAC = RutprevAC_ref * Coeff[i]
            TempRut = RutprevAC_ref
        rutAC[ii,:] = RutprevAC
        rutAC_ref[ii,:] = RutprevAC_ref
    return rutAC, rutAC_ref


###############################################################################
########################### UNBOUND RUTTING ###################################
###############################################################################

def RutUB(hB,ev, E_r, GWT, k_r1,n):#unbound material
    a_1 = 0.15
    b_1 = 0
    a_9 = 20
    b_9 = 0
    W_c = 51.712*(((E_r/2555)**(1/0.64))**(-0.3586*GWT**(0.1192)))
    C_0 = np.log((a_1*E_r**b_1)/(a_9*E_r**b_9))
    beta =10**( -0.61119 - 0.017638*W_c)
    ro = 10**9*((C_0)/(1-(10**9)**(beta)))**(1/beta)
    eober=((np.exp(ro**beta)*a_1*E_r**b_1)+(np.exp((ro/10**9)**beta)*a_9*E_r**b_9))/2
    Rut=hB*k_r1*eober*np.exp(-np.power((ro/n),beta))*ev
    return Rut

def NeqUB(Rut,hB,ev, E_r, GWT, k_r1):#E in psi, hac in inches
    a_1 = 0.15
    b_1 = 0
    a_9 = 20
    b_9 = 0
    W_c = 51.712*(((E_r/2555)**(1/0.64))**(-0.3586*GWT**(0.1192)))
    C_0 = np.log((a_1*E_r**b_1)/(a_9*E_r**b_9))
    beta =10**( -0.61119 - 0.017638*W_c)
    ro = 10**9*((C_0)/(1-(10**9)**(beta)))**(1/beta)
    eober=((np.exp(ro**beta)*a_1*E_r**b_1)+(np.exp((ro/10**9)**beta)*a_9*E_r**b_9))/2
    a=-np.log(Rut/hB/ev/eober/k_r1)
    a[a<1/1000] = 1/1000
    b=np.power(a,-1/beta)
    n=ro*np.abs(b)   
    n[a==1/1000]=-0.001
    #print(a)
    return n

def RutAccumUB(hB,ev, E_r, GWT, k_r1, n, Rutprev):
    if np.sum(Rutprev)==0:
        Ne=ev*0
    else:
        Ne = NeqUB(Rutprev,hB,ev, E_r, GWT, k_r1)
    Rd = RutUB(hB,ev, E_r, GWT, k_r1, Ne+n)
    Rd[np.where(Ne==-0.001)]=Rutprev[np.where(Ne==-0.001)]
    return Rd     

def Unboundrutting(struct, temperature, FER, Single, Tandem, InputMat, rutB, rutS,IP_npl,IP_pl):
    if struct == ['KSS','KTS'] or struct == ['KSW','KTW']:
        sublayerB = [1.71, 4.47, 5.82]
        sublayerS = [12.97, 18.80]
        E_br = 60000
    elif struct== ['NSS','NTS'] or struct ==['NSW','NTW']:
        sublayerB = [2.52, 3.83, 5.65]
        sublayerS = [10.95, 17.67]
        E_br = 40000
    for ii in range(len(sublayerB)):
        response = 'E22_Base_'+str(ii+1)+'_Z'
        E22B_ref = -deepcopy(FER.loc[:,response].values)
        E22B_ref[E22B_ref<1e-6] = 1e-6
        RutprevB = rutB[ii,:]
        RutprevB = RutAccumUB(sublayerB[ii],E22B_ref, E_br, 20, 2.03,5*sum(InputMat),RutprevB)        
        rutB[ii,:] = RutprevB
        
    if temperature == 'T2':
        E_r = 10000
    elif temperature == 'T3':
        E_r = 10000
    for ii in range(len(sublayerS)):
        response = 'E22_Subgrade_' + str(ii+1)+'_Z'
        E22SG_ref = -deepcopy(FER.loc[:,response].values)
        E22SG_ref[E22SG_ref< 1e-6] = 1e-6
        RutprevS = rutS[ii,:]
        RutprevS = RutAccumUB(sublayerS[ii],E22SG_ref,E_r, 20, 2.03,5*sum(InputMat), RutprevS)
        rutS[ii,:] = RutprevS
    return [rutB, rutS]



if __name__ == "__main__":
    Resp = pd.ExcelFile('FE_INPUTS - Final.xlsx').parse(['KSS','KTS'])
    Single = Resp['KSS']
    Tandem = Resp['KTS']
    Npl_Trucks = 18000
    Pl_groups = 4000
    PS = 3
    s = np.random.normal(0,10,Npl_Trucks)
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
    FER = ES.Expected(freq_npl,bin_npl,freq_pl,bin_pl,PS,Single,Tandem)
    
