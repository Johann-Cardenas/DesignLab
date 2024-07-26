# -*- coding: utf-8 -*-
"""
Created on Sat May 28 21:10:59 2022

@author: aravind4
"""

import numpy as np
import pandas as pd
import time
import Temperature_Parameters as TP
import Expected_strain as ES
import Rut_calc_RP_v4 as RC_V4
import Fatigue_calc as FC
import IRI_calc
start = time.time()
# =============================================================================
# INITIALIZING PARAMETERS
# =============================================================================
bin_npl = np.linspace(-24.5,24.5,99) 
bin_pl = np.array([[-24,-22],[-1,1],[22,24]])
prob_pl = np.array([1/3,1/3,1/3])
PS = 3 # Platoon Size
PL = 75 # Penetration Level in %
GR = 0.01                                 # Traffic growth rate
AP = 10                                    # Analysis period
RPFactor = 0.0096 # Need to write a function
# =============================================================================
# TRAFFIC AND STRUCTURE
# =============================================================================
Traffic = pd.ExcelFile('Traffic_info.xlsx')
Sheet_Names = Traffic.sheet_names
AMDTT, AADTT, Axles_per_TruckClass, Normalized_Truck_Dist, Single_Axle, Tandem_Axle = [Traffic.parse(i) for i in Sheet_Names]
P_si = Single_Axle['Bin_Center_kN'].values/1./2.
P_ta = Tandem_Axle['Bin_Center_kN'].values/1./2.
PP = [P_si, P_ta]
Structure = ['NSS','NTS']                      # Type of Structure, ThinS, ThinW, ThickS, ThickW
AADTT = AADTT.values[0][0]
AADTT_yearly = [AADTT*(1+GR)**i for i in range(AP)]
# =============================================================================
# INITIALIZING DISTRESS PARAMETERS
# =============================================================================
rutfinalAC = []
maxrutAC = [0]
rutfinalBase = []
maxrutBase = [0]
rutfinalSubg = []
maxrutSubg = [0]
rutfinal = [0]
rutfinalAC_sd = [0]
rutfinalBase_sd = [0]
rutfinalSubg_sd = [0]
rutfinal_sd = [0]
iteration = 0
FC_total_max = [0]
DI_BU = np.zeros((99,))
DI_TD = np.zeros((99,))
# =============================================================================
# MAIN CODE
# =============================================================================
for yy in range(AP):
    AADTT_y = AADTT_yearly[yy]
    print(yy)
    for mm in range(12):
        N_ycm_si, N_ycm_ta, N_ycm_tr = 0, 0, 0
        Temp = TP.temperature(mm)
        Resp = pd.ExcelFile('FE_INPUTS-'+Temp+'.xlsx').parse(Structure)
        Single = Resp[Structure[0]] # Single Axle Strain Profile for 12 kips axle load
        Tandem = Resp[Structure[1]] # Tandem Axle Strain Profile for 34 kips axle load
        
        for cc in range(10):
            AADTT_yc = AADTT_y * Normalized_Truck_Dist['Percent_AADTT'][cc]/100.
            TT_ycm = 30. * AADTT_yc * AMDTT['Class_' + str(cc+4)][mm]/100.
            Si = TT_ycm * Axles_per_TruckClass['Single'][cc]
            Ta = TT_ycm * Axles_per_TruckClass['Tandem'][cc]
            N_ycm_si += Si * Single_Axle['Class_' + str(cc+4)].values/100.
            N_ycm_ta += Ta * Tandem_Axle['Class_' + str(cc+4)].values/100.
        
        NN = [N_ycm_si, N_ycm_ta]
        iteration += 1        
        # =============================================================================
        # Truck parameters initialization        
        # =============================================================================
        Total_Trucks = np.sum(NN[0]) # Total Number of Trucks
        Pl_Trucks = int(Total_Trucks*PL/100) # Total Number of Trucks in platoon
        Pl_groups = int(Pl_Trucks/PS) # Total number of platoons
        Npl_Trucks = int(Total_Trucks - Pl_groups*PS) # Manual Traffic
        s_normal = np.random.normal(0,10,Npl_Trucks)
        b = np.histogram(s_normal,bins=np.linspace(-24.75, 24.75,100))
        freq_npl = b[0]
        freq_pl = []
        for i in range(len(bin_pl)):
            s = np.random.uniform(bin_pl[i][0],bin_pl[i][1],int(prob_pl[i]*Pl_groups))
            b = np.histogram(s,bins=np.linspace(bin_pl[i][0], bin_pl[i][1],6))
            freq_pl.append(b[0])
        freq_pl = np.array(freq_pl)
        IP_npl = np.concatenate([(bin_npl[i]*np.ones((freq_npl[i],))) 
                                       for i in range(len(freq_npl))])
        np.random.shuffle(IP_npl)
        locs = np.linspace(bin_pl[:,0],bin_pl[:,1],5).T
        IP_pl = np.concatenate([(locs[i][j]*np.ones((freq_pl[i][j],)))
                 for i in range(len(freq_pl)) for j in range(len(freq_pl[0]))])
        np.random.shuffle(IP_pl)
        InputMat = np.append(PS*np.ones(Pl_groups),np.ones(sum(freq_npl)))
        #np.random.shuffle(InputMat)
        # =============================================================================
        # Expected Strain
        # =============================================================================
        FER_rut_AC = ES.Expected(freq_npl, bin_npl, freq_pl, bin_pl, PS, Single, Tandem, RPFactor)
        cols = len(FER_rut_AC) # Transverse Strain Profile bins
        FER_rut_UB = ES.Expected(freq_npl, bin_npl, freq_pl, bin_pl, PS, Single, Tandem, 1)
        # =============================================================================
        # Rutting calculation        
        # =============================================================================
        if iteration == 1 and (Structure == ['NSS','NTS'] or Structure == ['NSW','NTW']):
            rutAC_D, rutAC_ref = np.zeros((5,cols)), np.zeros((5,cols))
            rutB_D = np.zeros((3,cols))
            rutS_D = np.zeros((2,cols))
        elif iteration == 1 and (Structure == ['KSS','KTS'] or Structure == ['KSW','KTW']):
            rutAC_D, rutAC_ref = np.zeros((7,cols)), np.zeros((7,cols))
            rutB_D = np.zeros((3,cols))
            rutS_D = np.zeros((2,cols))
        rutAC_D, rutAC_ref = RC_V4.ACrutting(Structure, Temp, FER_rut_AC, Single, Tandem, 
                                  InputMat, RPFactor, rutAC_D, rutAC_ref,IP_npl,IP_pl)                       
        rutfinalAC.append(np.sum(rutAC_D,0))
        maxrutAC.append(max(np.sum(rutAC_D,0)))
        rutB_D, rutS_D = RC_V4.Unboundrutting(Structure, Temp, FER_rut_UB, Single, Tandem, 
                                              InputMat, rutB_D, rutS_D,IP_npl,IP_pl)
        rutfinalBase.append(np.sum(rutB_D,0))
        maxrutBase.append(max(np.sum(rutB_D,0)))
        rutfinalSubg.append(np.sum(rutS_D,0))
        maxrutSubg.append(max(np.sum(rutS_D,0)))
        # =============================================================================
        # Fatigue Calculation
        # =============================================================================
        BU_Nf_s, BU_Nf_t = FC.Nf_cracking(Structure, Single, Tandem, InputMat, 'BU')
        TD_Nf_s, TD_Nf_t = FC.Nf_cracking(Structure, Single, Tandem, InputMat, 'TD')
        DI_bu, DI_td = ES.Expected_Nf(freq_npl, bin_npl, freq_pl, bin_pl, PS, 
                                    [1/0.44/BU_Nf_s,1/TD_Nf_s], [1/0.44/BU_Nf_t,1/TD_Nf_t])
        DI_BU += DI_bu
        DI_TD += DI_td
        FC_BU = FC.FC_BU_(DI_BU, Structure)
        FC_TD = FC.FC_TD_(DI_TD)
        DI_prev_BU, DI_prev_TD = DI_BU, DI_TD
        FC_total_max.append(max(FC_BU+FC_TD*100/5280.))
        
        
rutfinalAC = np.array(rutfinalAC)
maxrutAC = np.array(maxrutAC)
rutfinalBase = np.array(rutfinalBase)
maxrutBase = np.array(maxrutBase)
rutfinalSubg = np.array(rutfinalSubg)
maxrutSubg = np.array(maxrutSubg)
maxrut = maxrutAC + maxrutBase + maxrutSubg
rutfinal = rutfinalAC + rutfinalBase + rutfinalSubg
FC_total_max = np.array(FC_total_max)*100
IRI_max = IRI_calc.IRI_(FC_total_max, maxrut, np.linspace(0,AP,12*AP+1), 0)
Total_Combined = pd.DataFrame(np.array([maxrutAC, maxrut,FC_total_max,IRI_max]).T,
                              columns=['AC_rut','Total_rut','FC','IRI'])
end = time.time()
elapsed = (end - start)
print('Elapsed Time:', format(elapsed, '.2f'), 'seconds.')