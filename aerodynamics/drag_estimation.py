import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
from parameters import UAV
import matplotlib.pyplot as plt
import pandas as pd

#Importing needed stuff
aircraft = UAV('aircraft')
speed_of_sound = 328.387
V_c = aircraft.V_cruise
rho = 0.904637
MAC = aircraft.MAC_length


#Cf coefficients
def Cf(section):
    #calcultaing reynolds number and mach
    mu = 1.802E-5
    k = 0.634E-5
    Re = min(rho*V_c*MAC/mu, 38.21*((MAC/k)**1.053))
    M = V_c/speed_of_sound
    
    Cf_turb = 0.445/(((np.log10(Re))**2.58) *  (1+0.144*(M**2)**0.65))
    Cf_lam = 1.328/np.sqrt(Re)
    if section == 'wing' or 'htail' or 'vtail':
        Cf = 0.1*Cf_lam + 0.9*Cf_turb
    elif section == 'landing' or 'fuselage' or 'strut' or 'boom':
        Cf = Cf_turb
    else: 
        print('wrong section imported')

    return Cf
    

#FF 
def FF(section):    
    #calcultaing reynolds number and mach
    M = V_c/speed_of_sound
    if section == 'wing':
        type = 'lifting'
        airfoil = '4415'
    elif section == 'htail':
        type = 'lifting'
        airfoil = '0012'
    elif section == 'vtail':
        type = 'lifting'
        airfoil = '0009'
    elif section == 'landing':
        type = 'fus'
        l = 0.4
        d = 0.04
    elif section == 'fuselage':
        type = 'fus'
        l = 4.6342
        d = 0.67
    elif section == 'strut':
        type = 'fus'
        l = 2.56
        d = 0.015
    elif section == 'boom':
        type = 'fus'
        l = 2
        d = 0.05
    else: 
        print('wrong section imported')

    if type == 'lifting':
        x_c_max = float(airfoil[1]) / 10
        t_c = float(airfoil[2]+airfoil[3])
        FF = (1 + (0.6/x_c_max) * t_c + 100 * t_c**4) * (1.34 * M**0.18)
    elif type == 'fus':
        FF = 1 + 60/((l/d)**3) + (l/d)/400

    return FF

#IF
def IF(section):
    if section == 'wing':
        IF = 1
    elif section == 'htail' or 'vtail':
        IF = 1.04
    elif section == 'landing':
        IF = 1.1
    elif section == 'fuselage':
        IF = 1
    elif section == 'strut' or 'boom':
        IF = 1.1
    else: 
         print('wrong section imported')

    return IF

#Swet
def Swet(section):
    if section == 'wing':
        Swet = 1.07 * 2 * aircraft.Sw
    elif section == 'htail' or 'vtail':
        Swet = 1.05 * 2 * aircraft.Sh
    elif section == 'fuselage':
        D = 0.67
        L1 = 0.9342
        L2 = 2.9
        L3 = 0.8
        Swet = (np.pi * D / 4) * (1/(3*L1) * ((4*L1**2 + D**2 /4)**1.5 - D**3 /8) - D + 4*L2 + 2 * np.sqrt(L3**2 + D**2 /4))
    elif section == 'boom':
        #l*pi*r**2
        Swet = 2 * np.pi * 0.05**2
    elif section == 'strut':
        #2*l*pi*r**2 twice for both struts
        Swet = 2 * 2.56 * np.pi * 0.015**2
    elif section == 'landing':
        Swet = 3 * 0.4*np.sqrt(2) * 0.04**2 * np.pi
    else: 
         print('wrong section imported')

    return Swet

parts = ['wing', 'htail', 'vtail', 'landing', 'fuselage', 'strut', 'boom']

dataframe = {'Cf': Cf(parts), "FF": FF(parts), "IF": IF(parts), "Swet": Swet(parts)}
drag_info = pd.DataFrame(data=dataframe, index=parts)

print(drag_info)