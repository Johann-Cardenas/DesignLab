# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 13:04:43 2022

@author: aravind4
"""
import numpy as np

def IRI_(FC_total, RD, Age, TC):
    PI=1 #%
    FI=567
    p0075=8.7 #%
    p002=2 #%
    prec=42 #in
    SF=np.power(Age,1.5)*(np.log((prec+1)*(FI+1)*p002))+np.log((prec+1)*(PI+1)*p0075)
    IRIo = 1.0 * 1000 * 1.6 / 25.4 # 1 m/km IRI converted to in/mile
    C1, C2, C3, C4 = 40.0, 0.40, 0.008, 0.015
    return IRIo + C1 * RD + C2 * FC_total + C3*TC + C4*SF