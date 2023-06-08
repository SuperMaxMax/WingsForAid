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

#intermediate calculation for consistent updates in function
l_fus_1 = aircraft.l_fus_nose_cone
l_fus_2 = aircraft.l_fus_main_cone
l_fus_3 = aircraft.l_fus_tail_cone     
l_fus_total = l_fus_1 + l_fus_2 + l_fus_3
d_fus = aircraft.d_fuselage

variable_list2 = np.arange(0,3,0.5)


for j in range(2):
    nose = [False,True]
    for i in range(len(variable_list2)):
        if nose[j] == True:
            l_fus_1 = variable_list2[i]
            l_fus_3 = 0.8
        else:
            l_fus_1 = 0.9342
            l_fus_3 = variable_list2[i]
        
        l_fus_2 = 2.9
        l_fus_total = l_fus_1 + l_fus_2 + l_fus_3
        d_fus = aircraft.d_fuselage

        #Cf coefficients
        def Cf(section):
            #calcultaing reynolds number and mach
            mu = 1.802E-5
            k = 0.634E-5
            Re = min(rho*V_c*MAC/mu, 38.21*((MAC/k)**1.053))
            M = V_c/speed_of_sound
            
            Cf_turb = 0.445/(((np.log10(Re))**2.58) *  (1+0.144*(M**2)**0.65))
            Cf_lam = 1.328/np.sqrt(Re)

            if section == 'wing' or section == 'htail' or section == 'vtail':
                Cf = 0.1*Cf_lam + 0.9*Cf_turb
            elif section == 'landing' or section == 'fuselage' or section == 'strut' or section == 'boom':
                Cf = 0*Cf_lam + 1*Cf_turb 
            else: 
                print('wrong section imported')

            return Cf
            

        #FF 
        def FF(section):    
            #calcultaing reynolds number and mach
            M = V_c/speed_of_sound
            if section == 'wing':
                type = 'lifting'
                t_c = 0.15
                x_c_max = 0.4
            elif section == 'htail':
                type = 'lifting'
                t_c = 0.12
                x_c_max = 0.3
            elif section == 'vtail':
                type = 'lifting'
                t_c = 0.09
                x_c_max = 0.3
            elif section == 'landing':
                type = 'fus'
                l = 0.4
                d = 0.04
            elif section == 'fuselage':
                type = 'fus'
                l = l_fus_total
                d = d_fus
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
                FF = (1 + (0.6/x_c_max) * t_c + 100 * t_c**4) * (1.34 * M**0.18)
            elif type == 'fus':
                FF = 1 + 60/((l/d)**3) + (l/d)/400

            return FF

        #IF
        def IF(section):
            if section == 'wing':
                IF = 1
            elif section == 'htail' or section == 'vtail':
                IF = 1.04
            elif section == 'landing':
                IF = 1.1
            elif section == 'fuselage':
                IF = 1
            elif section == 'strut' or section == 'boom':
                IF = 1.1
            else: 
                print('wrong section imported')

            return IF

        #Swet
        def Swet(section):
            if section == 'wing':
                Swet = 1.07 * 2 * aircraft.Sw
            elif section == 'htail':
                Swet = 1.05 * 2 * aircraft.Sh_S * aircraft.Sw
            elif section == 'vtail':
                Swet = 1.05 * 2 * aircraft.Sv_S * aircraft.Sw
            elif section == 'fuselage':
                D = d_fus
                L1 = l_fus_1
                L2 = l_fus_2
                L3 = l_fus_3
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

        def Cd_misc(section):
            if section == 'landing':
                Cd = 2 * 0.55 *  aircraft.tire_main_height * aircraft.tire_main_width + 0.55 * aircraft.tire_nose_height * aircraft.tire_nose_width
            else:
                Cd = 0
            return Cd/aircraft.Sw

        #define all sections of the aircraft
        parts = ['wing', 'htail', 'vtail', 'landing', 'fuselage', 'strut', 'boom']

        #run all drag contributions per part
        Cf_list = []
        FF_list = []
        IF_list = []
        Swet_list = []
        CD0_list = []
        Cd_misc_list = []

        for part in parts:
            Cf_list.append(Cf(part))
            FF_list.append(FF(part))
            IF_list.append(IF(part))
            Swet_list.append(Swet(part))
            Cd_misc_list.append(Cd_misc(part))
            CD0_list.append((Cf(part) * FF(part) * IF(part) * Swet(part) / aircraft.Sw + Cd_misc(part)))

        # first_sum = sum(sum1_list)
        # second_sum = sum(Cd_misc_list)
        # leakage = (first_sum + second_sum)*0.05
        # CD_0 = first_sum + second_sum + leakage
        CD_0 = sum(CD0_list)*1.05
        print(CD_0)

        #save data to data frame
        dataframe = {'Cf': Cf_list, "FF": FF_list, "IF": IF_list, "Swet": Swet_list, 'Cd_misc': Cd_misc_list, 'CD_0': CD0_list}
        drag_info = pd.DataFrame(data=dataframe, index=parts)
        #drag_info = drag_info.round(3)
        #print(drag_info)
        print(l_fus_1, l_fus_3)
        print("Wetted area fuselage:", Swet_list[4])
        print("Zero lift drag fuselage:", CD0_list[4])
        print(l_fus_total / d_fus)
        plt.scatter(l_fus_total/d_fus, CD0_list[4], s = 4, color = 'black')

    plt.show()
drag_info.to_csv('drag_estimation.csv', index=True)
