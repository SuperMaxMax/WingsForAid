import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
from parameters import UAV
import pandas as pd
import csv
import matplotlib.pyplot as plt

aircraft = UAV('aircraft')

airfoil = aircraft.airfoil 

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

iw(airfoil)

def plot_lift_distr():
    segments = 10
    N = segments - 1
    S = 25 #aircraft.Sw
    AR = 8
    Lambda = 0.6
    alpha_twist = -1 * np.pi / 180
    i_w = 2 * np.pi / 180 #iw(airfoil)[0]
    a_2d = 6.3          #iw(airfoil)[1]
    alpha_0 = -1.5 *  np.pi / 180 #iw(airfoil)[2]
    b = (AR * S)**0.5
    MAC = S/b                               #Change to iteration between Croot and MAC
    Croot = (1.5*(1+Lambda)*MAC)/(1+Lambda+Lambda**2)
    theta = np.linspace(0, 90 * np.pi / 180, N, endpoint = False) #Change to get gooed amount of sections
    alpha = np.linspace(i_w, i_w + alpha_twist, N, endpoint = False)
    z = (b/2) * np.cos(theta)
    c = Croot * (1 - (1-Lambda) * np.cos(theta))
    mu = (c * a_2d) / (4 * b)
    #solve Ansin(ntheta)
    LHS = mu * (alpha - alpha_0)
    B = np.array((N,N))
    for i in range(1, N):
        for j in range(1, i):
            B[i,j] = np.sin((2*j - 1) * theta[i]) * (1 + (mu[i] * (2 * j - 1)) / np.sin(theta[i]))

    print(B)
        
    C_L_i = (4 * b)/c * 0#Ansin(ntheta)
    #C_L_wing = 
    plt.plot(z, C_L_i)
    plt.xlabel('semi span [m]')
    plt.ylabel('C_L')
    plt.show()
    
plot_lift_distr()
    
def fuel_volume(airfoil, Croot, Lambda):
    N = 50
    theta = np.linspace(0, 90 * np.pi / 180, N, endpoint = False) #Change to get gooed amount of sections
    chords = Croot * (1 - (1-Lambda) * np.cos(theta))
    m = float(airfoil[0]) / 100
    p = float(airfoil[1]) / 100
    for c in chords:
        x = np.linspace(0.25*c, 0.75*c, 100) #update when location spars known

        
    
