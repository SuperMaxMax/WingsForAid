import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(_file_), os.path.pardir)))

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

def lift(coefficients, theta, b, c): #Gives lift for certain point with position theta and c
    sum = 0
    for i in range(len(coefficients)):
        sum += coefficients[i]*np.sin((i+1) * theta)
    p = 4 * b / c * sum
    return p

def plot_lift_distr():
    segments = 10
    N = segments - 1
    S = 25 #aircraft.Sw
    AR = 8
    Lambda = 1
    alpha_twist = -1 * np.pi / 180
    i_w = 2 * np.pi / 180 #iw(airfoil)[0]
    a_2d = 6.3          #iw(airfoil)[1]
    alpha_0 = -1.5 *  np.pi / 180 #iw(airfoil)[2]
    b = (AR * S)**0.5
    MAC = S/b                               #Change to iteration between Croot and MAC
    
    N = 3
    a_2d = 6
    A = 8
    alpha_0 = -2 * np.pi / 180


    Croot = (1.5*(1+Lambda)*MAC)/(1+Lambda+Lambda**2)
   # if 
    theta = np.linspace(90 * np.pi / 180, 0, N, endpoint = False)[::-1] #Change to get gooed amount of sections
    #theta = np.array([np.pi/6,0, np.pi/3, 0,np.pi/2])
  #  theta = np.array([45 * np.pi / 180, 67.5 * np.pi / 180])
    alpha = np.linspace(i_w + alpha_twist, i_w, N, endpoint = False)[::-1]
 #   z = (b/2) * np.cos(theta)
    c = Croot * (1 - (1-Lambda) * np.cos(theta))
    mu = (c * a_2d) / (4 * b)
    #solve Ansin(ntheta)
  #  LHS = mu * (alpha - alpha_0)
    B = np.zeros((N,N))
    print(theta)
    for i in range(0,N):
        print("i:", i)
        for j in range(0,N):
           # print("j:", j)
            print("theta:", theta)
            print( ((4 * A / a_2d) + (j+1) / np.sin(theta[i])) * np.sin((j+1)* theta[i]))
            B[i,j] = ((4 * A / a_2d) + (j+1) / np.sin(theta[i])) * np.sin((j+1)* theta[i])
          #  B[i,j] = np.sin((2*j - 1) * theta[i]) * (1 + (mu[i] * (2 * j - 1)) / np.sin(theta[i]))
            
    #print(np.linalg.solve(B,np.ones(N)))

    coefficients = np.linalg.solve(B,alpha)
    print(B)

    thetas = np.arange(0.0001, np.pi/2, 0.01)
    z = (b/2) * np.cos(thetas)

    cs = Croot * (1 - (1-Lambda) * np.cos(thetas))
    C_L_i = lift(coefficients, thetas, b, cs)
    print("asdf")
    print(z)
    print(C_L_i)

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

        
    
