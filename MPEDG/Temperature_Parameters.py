# -*- coding: utf-8 -*-
"""
Created on Sun May 29 17:29:52 2022

@author: aravind4
"""

def TProfile(depth, temperature):
    h = depth*25
    if temperature == 'T1':
        T = (0.0002*h**2 -0.1384*h + 50.827)*1.8 + 32
    elif temperature == 'T2':
        T = (0.00007*h**2 -0.0334*h + 19.963)*1.8 + 32
    elif temperature == 'T3':
        T = (-0.0001*h**2 +0.0727*h -7.8002)*1.8 + 32
    return T

def temperature(mm):
    if mm == 10 or mm == 0 or mm == 4 or mm == 8 or mm == 9 or mm == 6 or mm == 7 or mm == 5:
        Temp = 'T2'
    elif mm==3 or mm==1 or mm == 2 or mm == 11:
        Temp = 'T2'    
    return Temp

