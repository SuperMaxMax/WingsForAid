import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize
from math import tan, atan, cos

def plot_loading_diagram(Ay, Az, By, Bz, load_case):
    plt.figure()
    y = np.linspace(0, aircraft.b/2, 1000)
    # add all forces
    plt.plot(y, [-spanwise_wing_weight(y)]*1000, label='Wing weight')
    plt.plot(y, [-spanwise_flap_weight()]*1000, label='Flap weight')
    plt.plot(y, [-spanwise_aileron_weight()]*1000, label='Aileron weight')
    plt.plot(y, [-spanwise_fuel_weight(y)*abs(load_case-1)]*1000, label='Fuel weight')
    plt.plot(y, spanwise_wing_loading(y)*load_case, label='Wing loading')
    # set axis limits such that all vector arrows are visible
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
aircraft.W_wl = 0.5                # [kg] for 1 winglet
plot = True

def calc_strut_weight(wsl):
    # Test for two loading cases -> 1: full wing loading no fuel, 2: no wing loading full fuel
    weights = []
    for i in range(2):
        if i == 0:
            load_case = 1
        else:
            load_case = 0

        # Strut_location
        aircraft.wing_strut_location = wsl # [m]
        aircraft.strut_angle = atan(1.1/aircraft.wing_strut_location)          # [deg]

        # Calculate reaction forces
        Bz = -1/aircraft.wing_strut_location*(quad(spanwise_wing_loading_ty, 0, aircraft.b/2)[0]*load_case-spanwise_flap_weight()*((aircraft.b/2*0.6)**2/2)+spanwise_aileron_weight()*aircraft.b/2*0.4*(0.8*aircraft.b/2)-quad(spanwise_fuel_weight_ty, 0, aircraft.b/2)[0]*abs(load_case-1)-quad(spanwise_wing_weight_ty, 0 , aircraft.b/2)[0]-aircraft.W_wl*aircraft.g0*aircraft.b/2)
        Az = (Bz+quad(spanwise_wing_loading, 0, aircraft.b/2)[0]*load_case - spanwise_flap_weight()*(aircraft.b/2*0.6) - spanwise_aileron_weight()*(aircraft.b/2*0.4) - quad(spanwise_fuel_weight, 0, aircraft.b/2)[0]*abs(load_case-1) - quad(spanwise_wing_weight, 0 ,aircraft.b/2)[0] - aircraft.W_wl*aircraft.g0)
        By = Bz/tan(aircraft.strut_angle)
        Ay = -By

        # Calculate strut weight
        sigma_yield = 683*10**6 # [Pa]
        mat_density = 2850      # [kg/m^3]
        E = 73.1*10**9          # [Pa]
        m_strut_t = aircraft.wing_strut_location/cos(aircraft.strut_angle)*(np.sqrt(By**2+Bz**2)/sigma_yield*mat_density)
        m_strut_c = (aircraft.wing_strut_location/cos(aircraft.strut_angle))**2*np.sqrt(np.sqrt(By**2+Bz**2)/(np.pi**2*E))*mat_density

        if load_case == 1:
            weights.append(m_strut_t)
        else:
            weights.append(m_strut_c)

        # Draw force diagram along span
        if plot:
            plot_loading_diagram(Az, Bz, By, Ay, load_case)
        
    return max(weights)

# # Optimize strut length for minimum weight
# res = minimize(calc_strut_weight, 2, bounds=((0, aircraft.b/2),), method='SLSQP')

# print(f"Optimal strut length: {res.x[0]} m")

print(calc_strut_weight(2.5))



