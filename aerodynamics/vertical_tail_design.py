# 24. Select the vertical tail volume coefficient, V V (Table 6.4).
# 25. Assume the vertical tail moment arm (lv) equal to the horizontal tail moment arm (l).
# 26. Calculate vertical tail planform area, Sv (Equation (6.74)).
# 27. Select vertical tail airfoil section (Section 6.8.2.4).
# 28. Select vertical tail aspect ratio, ARv(Section 6.8.2.6).
# 29. Select vertical tail taper ratio, λv(Section 6.8.2.7).
# 30. Determine the vertical tail incidence angle (Section 6.8.2.5).
# 31. Determine the vertical tail sweep angle (Section 6.8.2.8).
# 32. Determine the vertical tail dihedral angle (Section 6.8.2.9).
# 33. Calculate vertical tail span (bv), root chord (Cvroot), and tip chord (Cvtip ), and
# MACv(Equations (6.76)–(6.79)).
# 34. Check the spin recovery.
# 35. Adjust the location of the vertical tail relative to the horizontal tail by changing lv to
# satisfy the spin recovery requirements (Section 6.8.2.2).
# 36. Analyze directional trim (Section 6.8.1).
# 37. Analyze directional stability (Section 6.8.1).
# 38. Modify to meet the design requirements.
# 39. Optimize the tail.

import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from parameters import UAV
object = UAV('aircraft')

def vertical_wing_design(aircraft):
    #airfoil data
    #NACA 0009
    #data from https://digital.library.unt.edu/ark:/67531/metadc65459/m2/1/high_res_d/19930090937.pdf except alpha stall
    list_Cl0 = 0
    list_Cd_min = 0.0064
    list_Cm_0 = 0
    list_alpha_0 = 0
    list_alpha_s = 13.2
    list_Cl_max = 1.39
    list_Cl_alpha = 0.098 * 180/np.pi
    list_tc = 0.09
    
    Sv_Sw = aircraft.AE_Sv_S
    b = aircraft.AE_b
    l_v = aircraft.AE_l_h
    Sw = aircraft.AE_Sw

    V_v = Sv_Sw * l_v / b
    S_v = Sv_Sw * Sw
    print(V_v)


vertical_wing_design(object)

#All vertical tail volume things known so only influence planform sizing

#Function that determines required lift based on counter torque
def required_lift(aircraft):
    V_c = aircraft.V_cruise
    rho_c = aircraft.rho_cruise
    W_TO = aircraft.W_TO * 9.80665


    rpm = 3000 #Change to object variable [1/min]
    omega = rpm / 60 * 2 * np.pi
    power = 100*745.7 #Change to object variable [W]

    M = power / omega

    ver_dist = 2 #Change to object variable, distance between centre of pressure vertical tail and vertical cg location

    L_h = M / ver_dist

    Sw = aircraft.AE_Sw
    Sv_Sw = aircraft.AE_Sv_S
    S_v = Sv_Sw * Sw

    C_L_h = L_h / (0.5 * rho_c * V_c**2 * S_v)
    #C_L_h = 0.4


    C_L_W_c = 2*W_TO/(rho_c * (V_c**2) * Sw) #lift in cruise

    return C_L_h, C_L_W_c



#Sizing function
    #Incidence angle such that alpha is zero in cruise, so calc cruise torque by prop
    #Determine rest based on systems engineering book

def airfoil_select(C_L_h, change):
    #airfoil data
    #NACA0006, 0009, 0012
    #data from https://digital.library.unt.edu/ark:/67531/metadc65459/m2/1/high_res_d/19930090937.pdf except alpha stall
    list_Cl0 = [0, 0, 0]
    list_Cd_min = [0.0054, 0.0064, 0.0069]
    list_Cm_0 =[0, 0, 0]
    list_alpha_0 = [0, 0, 0]
    list_alpha_s = [11.0, 13.2, 16.4]
    list_Cl_max = [0.91, 1.39, 1.66]
    list_Cl_alpha = [0.098 * 180/np.pi, 0.098 * 180/np.pi, 0.099 * 180/np.pi]
    list_tc = [0.06, 0.09, 0.12]

    dataframe = {'C_l_0': list_Cl0, "Cdmin": list_Cd_min, "Cm0": list_Cm_0, "alpha_0": list_alpha_0, "alpha_s": list_alpha_s, "C_l_max": list_Cl_max, "C_l_alpha": list_Cl_alpha, "t/c": list_tc}
    df = pd.DataFrame(data=dataframe, index=["0006", "0009", "0012"]) #NACA

    airfoils = ["0006", "0009", "0012"] #NACA

    if abs(C_L_h) < 0.1:
        if change == 0 or change == -1:
            airfoil = airfoils[0]
        elif change == 1:
            airfoil = airfoils[0+change]
    elif abs(C_L_h) < 0.2:
        airfoil = airfoils[1+change]
    else:
        if abs(C_L_h) > 0.5:
            print("Required lift coefficient of horizontal too high something must be changed in the design to limit it. Currently C_L_h = ", C_L_h)
        else:
            if change == 0 or change == -1:
                airfoil = airfoils[2+change]
            else: 
                airfoil = airfoils[2]
    
    return df.loc[[airfoil]]

def horizontal_tail_planform(aircraft):
    def plot_lift_distr(i_w, full_print = False):
        variable = "Lambda"      #Lambda, AR or Twist
        plot_mode = "Normalize"         #"Normalized" for normalized plots
        if variable == "Lambda":    
            variable_list2 = [0.6]
        elif variable == "AR":  
            variable_list2 = [5.1666666]
        elif variable == "Twist":
            variable_list2 = [0]
        
        airfoildata_temp = airfoil_select(required_lift(aircraft)[0], 0)
        a_stall = airfoildata_temp['alpha_s'].tolist()[0] * np.pi /180
        
        for parameter in variable_list2:
            if abs(i_w/a_stall) > 0.1333: #2degree out of 15 is 0.13333
                change = 1
            elif abs(i_w/a_stall) < 0.06666: #1degree out of 15 is 0.06666
                change = -1
            else: 
                change = 0
            
            airfoildata = airfoil_select(required_lift(aircraft)[0], change)

            segments = 100
            N = segments - 1
            S = aircraft.AE_Sw * aircraft.AE_Sh_S  #tail.S
            if variable == "AR":
                AR = parameter
            else:
                AR = 5.1666666
            if variable == "Lambda":
                Lambda = parameter
            else:
                Lambda = 0.6
            if variable == "Twist":
                alpha_twist = parameter
            else:
                alpha_twist = 0 * np.pi / 180
            
            C_L_h = required_lift(aircraft)[0]
            i_w = i_w[0]
            a_2d = airfoildata['C_l_alpha'].tolist()[0]         #AIRFOIL PARAMETER                                  
            alpha_0 = airfoildata['alpha_0'].tolist()[0]        #AIRFOIL PARAMETER
            
            b = (AR * S)**0.5
            Croot = 2/(1+Lambda) * S/b
            MAC = Croot * 2 / 3 * ((1 + Lambda + Lambda**2)/(1+Lambda))

            theta = np.linspace(np.pi/(2*N), np.pi/2, N, endpoint = True) 
            alpha = np.linspace(i_w+alpha_twist, i_w, N, endpoint = False)
            z = (b/2) * np.cos(theta)
            c = Croot * (1 - (1-Lambda) * np.cos(theta))
            mu = (c * a_2d) / (4 * b)

            #solve Ansin(ntheta)
            LHS = mu * (alpha - alpha_0)
            RHS = np.zeros((N,N))
            for i in range(N):
                for j in range(N):
                    RHS[i,j] = np.sin((2*j + 1) * theta[i]) * (1 + (mu[i] * (2 * j + 1)) / np.sin(theta[i]))
            A = np.linalg.solve(RHS, LHS)
            sum = np.zeros(N)
            for i in range(N):
                for j in range(N):
                    sum[i] += A[j] * np.sin((2 * j + 1) * theta[i])


            CL = 4 * b * sum / c
            CL1 = np.insert(CL, 0, 0)
            y_s = np.insert(z, 0, b / 2)
            if plot_mode == "Normalized":
                CL1 = CL1/max(CL1)
                y_s = y_s / (b/2)
                
        
            label = variable + "= " + str(parameter)
            plt.plot(y_s, CL1, marker = "s", label = label)

            ##Wing Lift Coefficient
            C_L_wing = np.pi * AR * A[0]
            
            ##Wing INDUCED DRAG
            cdi_sum = 0
            for i in range(len(A)):
                cdi_sum += (i+1) * A[i]**2
            CD_induced = np.pi * AR *  cdi_sum #(Wing) 

            #Span efficiency factor
            delta = 0
            for i in range(1,len(A)):
                delta += (i+1) * (A[i] / A[0])**2
            span_eff = 1 / (1 + delta)
            tau = 1/span_eff - 1
            
            CL_a_h = a_2d / (1+(a_2d/(np.pi*AR))*(1+tau))


            #print('=====================================================================')
            #print('current option is: AR = ', AR, 'taper ratio = ', Lambda, 'indidence = ', i_w*180/np.pi)
            #print("Span_eff = ", span_eff, "CL_wing = ", C_L_wing, "CL required for cruis = ", C_L_h, "CD_i = ", CD_induced)
            print(i_w*180/np.pi, C_L_wing - C_L_h, airfoildata.index.tolist())
            #print("C_L", C_L_wing)

        #Find integral current distribution
        area_lift_dist = -integrate.simps(CL1, y_s)

        #Elliptical lift distribution
        y = np.linspace(0, b/2, 50, endpoint = True)
        #Cli_elliptical = (b/2) * np.sqrt(1-(((np.pi * b * y)/(8 * area_lift_dist))**2))
        Cli_elliptical = (8*area_lift_dist)/(np.pi*b) * np.sqrt(1-(2*y/b)**2)
        #plt.plot(y, Cli_elliptical, label = "Elliptical")

        #General plot
        plt.grid()
        plt.xlabel('semi span [m]')
        plt.ylabel('C_L')
        plt.legend()
        #plt.show()
        
        if not full_print:
            return abs(C_L_wing - C_L_h)
        elif full_print:
            return AR, b, Lambda, alpha_twist, S, CL_a_h  #choose whatever

    airfoildata = airfoil_select(required_lift(aircraft)[0], 0)
    C_l_alpha = airfoildata['C_l_alpha'].tolist()[0]
    initial_guess = required_lift(aircraft)[0]/C_l_alpha
    a_h_optimal = optimize.minimize(plot_lift_distr,initial_guess, method = 'Nelder-Mead', tol=1e-06)['x']
    print(a_h_optimal)
    AR, b, Lambda, alpha_twist, S, CL_a_h  = plot_lift_distr(a_h_optimal, full_print = True)


horizontal_tail_planform(object)
