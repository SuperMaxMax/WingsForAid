import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from math import tan, atan, cos
from material_dictionary import material, rank_material

def plot_loading_diagram(Ay, Az, By, Bz, load_case):
    plt.figure()
    y = np.linspace(0, aircraft.b/2, 1000)
    # add all forces with hatch pattern
    plt.plot(y, [-spanwise_wing_weight(y)]*1000, label='Wing weight')
    plt.plot(y, [-spanwise_flap_weight()]*1000, label='Flap weight')
    plt.plot(y, [-spanwise_aileron_weight()]*1000, label='Aileron weight')
    plt.plot(y, [-spanwise_fuel_weight(y)*abs(load_case-1)]*1000, label='Fuel weight')
    plt.plot(y, spanwise_wing_loading(y)*load_case, label='Wing loading')
    # set axis limits such that all vector arrows are visible, vector are scaled based on axis limits
    if load_case == 1:
        plt.title('Spanwise force diagram for full wing loading, no fuel')
        global y_range1
        y_range1 = plt.axis()[3]-plt.axis()[2]
        hl1, hw1 = 0.06, 20
        hl2, hw2 = 20, 0.06
        w1, w2 = 5, 0.02
        global fac
        fac = 100
    else:
        plt.title('Spanwise force diagram for no wing loading, full fuel')
        y_range2 = plt.axis()[3]-plt.axis()[2]
        hl1, hw1 = 0.06, 20*y_range2/y_range1
        hl2, hw2 = 20*y_range2/y_range1, 0.06
        w1, w2 = 5*y_range2/y_range1, 0.02
        fac *= y_range2/y_range1
    # add reaction forces
    if np.sign(Ay) == 1:
        plt.arrow(-np.sign(Ay)*0.5, 0, np.sign(Ay)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.arrow(aircraft.wing_strut_location-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
    else:
        plt.arrow(-np.sign(Ay)*0.5, 0, np.sign(Ay)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.arrow(aircraft.wing_strut_location-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
    if np.sign(Az) == 1:
        plt.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    else:
        plt.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    if np.sign(Bz) == 1:
        plt.arrow(aircraft.wing_strut_location, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    else:
        plt.arrow(aircraft.wing_strut_location, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    plt.xlabel('Spanwise location [m]')
    plt.ylabel('Force [N/m]')
    plt.legend()
    plt.grid()
    plt.show()

def spanwise_wing_weight(y):
    # kg / m
    return (aircraft.W_w/aircraft.b)*aircraft.g0

def spanwise_wing_weight_ty(y):
    # kg / m
    return spanwise_wing_weight(y)*y

def spanwise_wing_loading(y):
    # put if for worst case scenarios
    gamma_0 = (aircraft.W_TO-aircraft.W_F)*aircraft.g0*(4/np.pi)*(1/aircraft.rho_cruise)*(1/aircraft.V_cruise)*(1/aircraft.b)
    return aircraft.rho_cruise*aircraft.V_cruise*gamma_0*np.sqrt(1-(2*y/aircraft.b)**2)

def spanwise_wing_loading_ty(y):
    # put if for worst case scenarios
    return spanwise_wing_loading(y)*y

def spanwise_flap_weight():
    # kg / m
    return (aircraft.W_sc/2)*0.6/aircraft.b*aircraft.g0

def spanwise_aileron_weight():
    # kg / m
    return (aircraft.W_sc/2)*0.4/aircraft.b*aircraft.g0

def spanwise_fuel_weight(y):
    # kg / m
    return aircraft.W_F/aircraft.b*aircraft.g0

def spanwise_fuel_weight_ty(y):
    # kg / m
    return spanwise_fuel_weight(y)*y

aircraft = UAV('aircraft')

# Strut_location
aircraft.wing_strut_location = 0.437277492*aircraft.b/2                                # [m] | based on statistical data, value is now based on half span running from middle fuselage
aircraft.strut_angle = atan(1.1/(aircraft.wing_strut_location-aircraft.w_out/2))    # [rad]
aircraft.l_strut = aircraft.wing_strut_location/np.cos(aircraft.strut_angle)        # [m]
aircraft.W_wl = 0.5                                                                 # [kg] for 1 winglet
plot = False

# Material choice for strut weight
prop_weights_strut = [0.1, 0.5, 0.05, 2, 0.05 ,0.05, 0.05]
possible_materials = rank_material(prop_weights_strut)[:5]
ms = pd.DataFrame(index=possible_materials, columns=['Case 1', 'Case 2'])

for i in range(2):
    if i == 0:
        load_case = 1
    else:
        load_case = 0

    # Calculate reaction forces
    Bz = -1/aircraft.wing_strut_location*(quad(spanwise_wing_loading_ty, 0, aircraft.b/2)[0]*load_case-spanwise_flap_weight()*((aircraft.b/2*0.6)**2/2)+spanwise_aileron_weight()*aircraft.b/2*0.4*(0.8*aircraft.b/2)-quad(spanwise_fuel_weight_ty, 0, aircraft.b/2)[0]*abs(load_case-1)-quad(spanwise_wing_weight_ty, 0 , aircraft.b/2)[0]-aircraft.W_wl*aircraft.g0*aircraft.b/2)
    Az = (Bz+quad(spanwise_wing_loading, 0, aircraft.b/2)[0]*load_case - spanwise_flap_weight()*(aircraft.b/2*0.6) - spanwise_aileron_weight()*(aircraft.b/2*0.4) - quad(spanwise_fuel_weight, 0, aircraft.b/2)[0]*abs(load_case-1) - quad(spanwise_wing_weight, 0 ,aircraft.b/2)[0] - aircraft.W_wl*aircraft.g0)
    By = Bz/tan(aircraft.strut_angle)
    Ay = -By
    
    P_truss = np.sqrt(By**2+Bz**2)
    P_truss *= aircraft.ST_SF

    # Plot force diagram
    if plot:
        plot_loading_diagram(Az, Bz, By, Ay, load_case)

    for mat in possible_materials:
        sigma_yield = material[mat]['yield stress']*10**6   # [Pa]
        mat_density = material[mat]['density']              # [kg/m^3]
        E = material[mat]['E']*10**9                        # [Pa]

        a = 0.0125       # semi_minor of ellipse [m]
        b = 2*a          # semi_major of ellipse [m]

        if load_case == 1:
            A_min = P_truss / sigma_yield
            m = A_min*mat_density*aircraft.l_strut
            ms.loc[mat, 'Case 1'] = m
        else:
            I_min = P_truss * aircraft.l_strut**2 / (np.pi**2 * E)
            t = 4*I_min/(np.pi*a**3 * (1+(3*b/a)))  # thickness of strut sheet [m]
            A = np.pi*(a+b)*t
            I_check = np.pi*a**3*t*(1+(3*b/a))/4
            m = A*mat_density*aircraft.l_strut
            ms.loc[mat, 'Case 2'] = m

print(ms)



