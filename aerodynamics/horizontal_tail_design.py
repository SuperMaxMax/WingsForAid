import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from parameters import UAV


def required_lift(aircraft):
    #load in parameters from current design
    AR = aircraft.A
    sweep_c4 = aircraft.sweep_co4
    taper = aircraft.taper
    C_m_af = aircraft.af_cm0
    alpha_t = aircraft.wing_twist
    V_c = aircraft.V_cruise
    rho_c = aircraft.rho_cruise
    W_TO = aircraft.W_TO * 9.80665
    S = aircraft.Sw
    h0 = aircraft.MAC_ac
    MAC = aircraft.MAC_length
    xcg_fwrd = aircraft.X_cg_fwd
    xcg_aft = aircraft.X_cg_aft
    Sh_SW = aircraft.Sh_S

    #test to verify code
    # AR = 28
    # sweep_c4 = 0
    # taper = 0.8
    # C_m_af = -0.013
    # alpha_t = -1.1 
    # V_c = 48.872222
    # rho_c = 0.905
    # W_TO = 850 * 9.81
    # S = 18
    # h0 = 0.23
    # MAC = 0.8

    #temp values from book
    #V_H = 0.7 #table 6.4 "aircraft design synthesis a systems engineering approach"
    V_H = Sh_SW * aircraft.l_h / MAC

    #inbetween calculutions 
    C_L = 2*W_TO/(rho_c * (V_c**2) * S) #lift in cruise
    sweep_LE = np.arctan(np.tan(sweep_c4) - (4/AR) * (25/100 * (1-taper)/(1+taper))) #leading edge sweep
    #sweep_LE = 8 * np.pi/180 #test to verify code
    C_m_0_wf = C_m_af * (AR * (np.cos(sweep_LE) ** 2)) / (AR + 2 * np.cos(sweep_LE)) + 0.01 * alpha_t

    most_extreme_cg = [xcg_fwrd, xcg_aft]
    C_L_h = 0 
    for xcg_wf in most_extreme_cg:
        h = xcg_wf
        #h = 0.114 #test to verify code
   
        C_L_h_new = (C_m_0_wf + C_L * (h - h0)) / (V_H)
        if abs(C_L_h_new) > abs(C_L_h):
            C_L_h = C_L_h_new
    C_L_W_c = C_L
    
    return C_L_h, C_L_W_c

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
#print(airfoil_select(required_lift()))

def horizontal_tail_planform(aircraft):
    def plot_lift_distr(i_w, full_print = False):
        variable = "Lambda"      #Lambda, AR or Twist
        plot_mode = "Normalize"         #"Normalized" for normalized plots
        if variable == "Lambda":    
            variable_list2 = [1]
        elif variable == "AR":  
            variable_list2 = [5.1666666]
        elif variable == "Twist":
            variable_list2 = [-1 * np.pi / 180 , -2 * np.pi / 180, -3 * np.pi / 180, -4 * np.pi / 180, -5 * np.pi / 180]
        
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
            S = aircraft.Sw * aircraft.Sh_S  #tail.S
            if variable == "AR":
                AR = parameter
            else:
                AR = 4
            if variable == "Lambda":
                Lambda = parameter
            else:
                Lambda = 1
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
            print("CLAH", CL_a_h)


            print('=====================================================================')
            print('current option is: AR = ', AR, 'taper ratio = ', Lambda, 'indidence = ', i_w*180/np.pi)
            print("Span_eff = ", span_eff, "CL_wing = ", C_L_wing, "CL required for cruis = ", C_L_h, "CD_i = ", CD_induced)
            print("aircraft width, container is 2.40 = ", b)
            #print(i_w*180/np.pi, C_L_wing - C_L_h, airfoildata.index.tolist())

        #Find integral current distribution
        area_lift_dist = -integrate.simps(CL1, y_s)

        #Elliptical lift distribution
        y = np.linspace(0, b/2, 50, endpoint = True)
        #Cli_elliptical = (b/2) * np.sqrt(1-(((np.pi * b * y)/(8 * area_lift_dist))**2))
        Cli_elliptical = (8*area_lift_dist)/(np.pi*b) * np.sqrt(1-(2*y/b)**2)
        #plt.plot(y, Cli_elliptical, label = "Elliptical")

        #print('lestgo')
        #General plot
        # plt.grid()
        # plt.xlabel('semi span [m]')
        # plt.ylabel('C_L')
        # plt.legend()
        # plt.show()
        
        if not full_print:
            return abs(C_L_wing - C_L_h)
        elif full_print:
            return AR, b, Lambda, alpha_twist, S, CL_a_h  #choose whatever

    airfoildata = airfoil_select(required_lift(aircraft)[0], 0)
    C_l_alpha = airfoildata['C_l_alpha'].tolist()[0]
    initial_guess = required_lift(aircraft)[0]/C_l_alpha
    a_h_optimal = optimize.minimize(plot_lift_distr,initial_guess, method = 'Nelder-Mead', tol=1e-06)['x']
    #print(a_h_optimal)
    AR, b, Lambda, alpha_twist, S, CL_a_h  = plot_lift_distr(a_h_optimal, full_print = True)

    #Adding effect of downwash
    C_L_W_c = required_lift(aircraft)[1]
    AR_w = aircraft.A
    i_w = aircraft.i_w
    a_f = aircraft.alpha_f
    C_L_a_w = aircraft.CL_a_w
    epsilon_0 = 2 * C_L_W_c / (np.pi * AR_w)
    epsilon_alpha = 2 * C_L_a_w / (np.pi * AR_w)
    epsilon = epsilon_0 + epsilon_alpha * i_w
    i_h = a_h_optimal - a_f + epsilon
    
    # Horizontal tailplane
    aircraft.S_h = S
    aircraft.A_h = AR                        # Aspect ratio horizontal tail. NOTE: This is a guestimate  
    aircraft.lambda_co2_h = 0               # [rad] Half chord sweep of horizontal tailplane [-] NOTE: This is a guestimate  
    aircraft.dEpsilondA = epsilon_alpha              # Downwash [-] TODO: check this value, this is a pure guess
    aircraft.A_h = AR                        
    aircraft.b_h = b                     
    aircraft.i_w_h = i_h       
    aircraft.wing_twist_h = alpha_twist    
    aircraft.sweep_co4_h = 0.0                 # Updated half chord sweep [rad]
    aircraft.sweep_co2_h = np.arctan(np.tan(aircraft.sweep_co4_h) - 4/AR * (25/100*(1-Lambda)/(1+Lambda))) 
    aircraft.sweep_LE_h = np.arctan(np.tan(aircraft.sweep_co4_h) - 4/AR * (-25/100*(1-Lambda)/(1+Lambda)))          
    aircraft.taper_h = Lambda                
    aircraft.rootchord_h = 2 * aircraft.S_h / (aircraft.b_h * (1+Lambda))            
    aircraft.tipchord_h = aircraft.rootchord_h*Lambda        
    aircraft.MAC_length_h = 2/3 * aircraft.rootchord_h * (1 + Lambda + Lambda**2) / (1 + Lambda)        
    aircraft.y_mac_h = 1/3*(aircraft.b_h/2)*(1+2*Lambda)/(1+Lambda)   
    aircraft.x_lemac_h = aircraft.y_mac_h*np.tan(aircraft.sweep_LE_h)
    aircraft.CL_a_h = CL_a_h
    #print("Afsd", CL_a_h)
    #print(b)
    #print(S)
    return 

