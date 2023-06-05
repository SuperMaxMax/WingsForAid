import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
from parameters import UAV
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize





def iw(airfoil):
    ### FIND CL optimal ###
    file_name_clcd = "cl-cd-" + str(airfoil)
    file_path_clcd = "aerodynamics/" + file_name_clcd + ".csv"
    
    cd_list = [] 
    cl_list = []

    with open(file_path_clcd) as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        linecount = 0
        for row in reader:
            if row:    
                #skip first line with explanation of columns
                if linecount == 0:
                    linecount += 1
                else:
                    cd_list.append(float(row[0]))
                    cl_list.append(float(row[1]))
    
    cd_bucket = []
    cl_bucket = []
    for i in range(len(cd_list)-1):
        if abs((cd_list[i] - cd_list[i+1]) / cd_list[i]) < 0.02: 
            cd_bucket.append(cd_list[i])
            cl_bucket.append(cl_list[i])
    cl_optimal = cl_bucket[len(cl_bucket)//2]
    
    ### FIND CL ALPHA ###
    file_name_clalpha = "cl-alpha-" + str(airfoil)
    file_path_clalpha = "aerodynamics/" + file_name_clalpha + ".csv"

    with open(file_path_clalpha) as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        linecount = 0
        for row in reader:
            if row:
                #skip first line with explanation of columns    
                if linecount == 0:
                    linecount += 1
                else:
                    alpha = float(row[0])
                    cl = float(row[1])
                    if alpha == 0:
                        cl_0 = cl
                    elif alpha == -4: 
                        cl_min5 = cl
                    elif alpha == 4:
                        cl_plus5 = cl

    cl_alpha = (cl_plus5 - cl_min5)/8

    ### FIND Iw ###
    Iw = (cl_optimal - cl_0)/cl_alpha #in degrees
    alpha_zero_lift = -cl_0/cl_alpha #in degrees
    
    Iw = Iw * np.pi/180
    cl_alpha = cl_alpha * 180/np.pi
    alpha_zero_lift = alpha_zero_lift * np.pi/180
    
    return Iw, cl_alpha, alpha_zero_lift

#iw(airfoil)

def main_wing_planform(aircraft):
    def plot_lift_distr(i_w, full_print = False):
        i_w = i_w[0]
        variable = "Lambda"      #Lambda, AR or Twist
        plot_mode = "Normalize"         #"Normalized" for normalized plots
        if variable == "Lambda":    
            variable_list2 = [0.4]
        elif variable == "AR":  
            variable_list2 = [7,7.25,7.5,7.75,8]
        elif variable == "Twist":
            variable_list2 = [-1 * np.pi / 180 , -2 * np.pi / 180, -3 * np.pi / 180, -4 * np.pi / 180, -5 * np.pi / 180]

        for parameter in variable_list2:
            segments = 30
            N = segments - 1
            S = aircraft.AE_Sw #aircraft.Sw
            if variable == "AR":
                AR = parameter
            else:
                AR = 7.75
            if variable == "Lambda":
                Lambda = parameter
            else:
                Lambda = 0.4
            if variable == "Twist":
                alpha_twist = parameter
            else:
                alpha_twist = 0 * np.pi / 180


            a_2d = aircraft.AE_cl_alpha       #iw(airfoil)[1]
            alpha_0 = aircraft.AE_alpha0 #iw(airfoil)[2]
            b = (AR * S)**0.5

            Croot = 2/(1+Lambda) * S/b
            MAC = Croot * 2 / 3 * ((1 + Lambda + Lambda**2)/(1+Lambda))                             #Change to iteration between Croot and MAC
        #  MAC = S / b
        #  Croot =  (1.5*(1+Lambda)*MAC)/(1+Lambda+Lambda**2)
            #print(MAC)
            theta = np.linspace(np.pi/(2*N), np.pi/2, N, endpoint = True) #Change to get gooed amount of sections
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

            #print(LHS)
            #print(RHS/mu)
            A = np.linalg.solve(RHS, LHS)
            #print(A)

            sum = np.zeros(N)
            for i in range(N):
                for j in range(N):
                    sum[i] += A[j] * np.sin((2 * j + 1) * theta[i])


            CL = 4 * b * sum / c
            CL1 = np.insert(CL, 0, 0)
            c = np.insert(c,0, Croot*Lambda)
            y_s = np.insert(z, 0, b / 2)

            if plot_mode == "Normalized":
                CL1 = CL1/max(CL1)
                y_s = y_s / (b/2)
                
        
            label = variable + "= " + str(parameter)
        # plt.plot(y_s, CL1, marker = "s", label = label)

            ##Wing Lift Coefficient
            C_L_wing = np.pi * AR * A[0]
            V_c = aircraft.V_cruise
            rho_c = aircraft.rho_cruise
            W_TO = aircraft.W_TO * 9.80665
            C_L_req = 2*W_TO/(rho_c * (V_c**2) * S)
            

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
            oswald = span_eff * 0.971 * 0.804 #from https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf

            tau = 1/span_eff - 1

            
            CL_a_W = a_2d / (1+(a_2d/(np.pi*AR))*(1+tau))

            #print('=====================================================================')
            #print('current option is: AR = ', AR, 'taper ratio = ', Lambda, 'indidence = ', i_w*180/np.pi)
            #print("Span_eff = ", span_eff, "CL_wing = ", C_L_wing, "CL required for cruis = ", C_L_req, "CD_i = ", CD_induced)
            #print(C_L_wing**2 / (AR* np.pi * CD_induced))

            #print("CL_wing", C_L_wing)
            #print("CL required for cruis", C_L_req)
            #print("CDi_wing", CD_induced)
            q = 0.5 * aircraft.rho_cruise * aircraft.V_cruise**2
            #print(CL1, c)
            #print(Croot*Lambda)
            l = CL1*c*q
            plt.plot(y_s, l)
            
            #Getting coefficients for a polynomial fit for the lift distribution
            # lift_coefficients = np.polyfit(y_s, l, 10)
                
            # lift_coefficients = lift_coefficients[::-1]
            # x = np.arange(-0.3, b/2, 0.001)
            # def polynomial(coefficients, x):
            #     y = 0
            #     for i in range(len(coefficients)):
            #         y+= coefficients[i]*x**i
            #         print(coefficients[i], i)

            #     return y
                    
            # y = polynomial(lift_coefficients, x)
            
            # plt.plot(x, y)
            # plt.show()

            # print("asdf", lift_coefficients)
            #print(-2*integrate.simps(l, y_s))
            #print(W_TO)
            #print(l)

        #Find integral current distribution
        area_lift_dist = -integrate.simps(CL1, y_s)
        #print(area_lift_dist)
        

        #Elliptical lift distribution
        y = np.linspace(0, b/2, 50, endpoint = True)
        #Cli_elliptical = (b/2) * np.sqrt(1-(((np.pi * b * y)/(8 * area_lift_dist))**2))
        Cli_elliptical = (8*area_lift_dist)/(np.pi*b) * np.sqrt(1-(2*y/b)**2)
    #   plt.plot(y, Cli_elliptical, label = "Elliptical")

        #General plot
        #plt.grid()
        #plt.xlabel('semi span [m]')
        #plt.ylabel('C_L')
        #plt.legend()
        #plt.show()
        
        if not full_print:
            return abs(C_L_wing-C_L_req)
        elif full_print:
            return AR, Lambda, alpha_twist, span_eff, CD_induced, i_w, tau, CL_a_w
    
    airfoil = aircraft.airfoil
    initial_guess = iw(airfoil)[0]
    i_w_optimal = optimize.minimize(plot_lift_distr,initial_guess, method = 'Nelder-Mead', tol=1e-06)['x']
    AR, Lambda, alpha_twist, span_eff, i_w, tau, CL_a_w = plot_lift_distr(i_w_optimal, full_print=True)

    aircraft.AE_A = AR                        
    aircraft.AE_b = (AR*aircraft.AE_Sw)**0.5                      
    aircraft.AE_span_eff = span_eff                     
    tau = 1/span_eff - 1
    aircraft.AE_CL_a_W = CL_a_w             
    aircraft.AE_tau = tau
    aircraft.AE_i_w = i_w       
    aircraft.AE_wing_twist = alpha_twist    
    aircraft.AE_sweep_co2 = 1 / np.tan(tan(aircraft.AE_sweep_co4) - 4/AR * (25/100*(1-Lambda)/(1+Lambda))) 
    aircraft.AE_sweep_LE = 1 / np.tan(tan(aircraft.AE_sweep_co4) - 4/AR * (-25/100*(1-Lambda)/(1+Lambda)))          
         
    aircraft.AE_taper = Lambda                
    aircraft.AE_rootchord = 2 * aircraft.AE_Sw / (aircraft.AE_b * (1+Lambda))            
    aircraft.AE_tipchord = aircraft.AE_rootchord*Lambda        
    aircraft.AE_MAC_length = 2/3 * aircraft.AE_rootchord * (1 + Lambda + Lambda**2) / (1 + Lambda)        
    aircraft.AE_y_mac = 1/3*(aircraft.AE_b/2)*(1+2*Lambda)/(1+Lambda)   
    aircraft.AE_x_lemac = aircraft.AE_y_mac/np.tan(aircraft.AE_sweep_LE)
               
    return 

aircraft =  UAV('aircraft')
#main_wing_planform(aircraft)
#print(aircraft.AE_A)

def fuel_volume(airfoil, Croot, Lambda, b):
    if len(airfoil) != 4: 
        print('fuel volume calculation is only valid for 4-digit NACA as of now')
    V = 0 
    N = 50 #number of span direction sections
    M = 50 #number of chord direction sections per span location
    for i in range(N):
        S = 0 
        c = Croot * (1 - (1-Lambda) * i/N)
        m = float(airfoil[0]) / 100
        p = float(airfoil[1]) / 10
        X = np.linspace(0, 1, M) #update when location spars known
        dx = 1/M
        for x in X:
            if x > 0.25 and x <= p: 
                y = m/(p**2) * (2*p*x - x**2)
                S += y*dx * c**2
            elif x > p and x < 0.75: 
                y =  m/((1 - p)**2) * ((1 - 2*p) + 2*p*x - x**2) #check wiki
                S += y*dx * c**2
        V += S * (b/2)/N #half span volume in m3
    V = V * 1000 #in liters
    return V 

print(fuel_volume(aircraft.airfoil, aircraft.rootchord, 0.4, aircraft.b))