import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def spanwise_wing_weight():
    # kg / m
    return (aircraft.W_w/aircraft.b)*aircraft.g0

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

def spanwise_fuel_weight():
    # kg / m
    return aircraft.W_F/aircraft.b*aircraft.g0

aircraft = UAV('aircraft')

# Strut_location
aircraft.wing_strut_location = 1.5 # [m]
aircraft.strut_angle = 45          # [deg]
aircraft.W_wl = 0.5                # [kg] for 1 winglet
plot = False

# Draw force diagram along span
if plot:
    plt.figure()
    y = np.linspace(0, aircraft.b/2, 1000)
    plt.plot(y, [spanwise_wing_weight(aircraft, y)]*1000, label='Wing weight')
    plt.plot(y, [spanwise_flap_weight(aircraft, y)]*1000, label='Flap weight')
    plt.plot(y, [spanwise_aileron_weight(aircraft, y)]*1000, label='Aileron weight')
    plt.plot(y, [spanwise_fuel_weight(aircraft, y)]*1000, label='Fuel weight')
    plt.plot(y, spanwise_wing_loading(aircraft, y), label='Wing loading')
    plt.xlabel('Spanwise location [m]')
    plt.ylabel('Force [N/m]')
    plt.legend()
    plt.show()

# Create loading diagrams
Bz = -1/aircraft.wing_strut_location*(quad(spanwise_wing_loading_ty, 0, aircraft.b/2)[0]-spanwise_flap_weight()*((aircraft.b/2*0.6)**2/2)+spanwise_aileron_weight()*aircraft.b/2*0.4*(0.8*aircraft.b/2)-spanwise_fuel_weight()*aircraft.b/2*aircraft.b/4-spanwise_wing_weight()*aircraft.b/2*aircraft.b/4-aircraft.W_wl*aircraft.g0*aircraft.b/2)
Az = -(Bz+quad(spanwise_wing_loading, 0, aircraft.b/2)[0] - spanwise_flap_weight()*(aircraft.b/2*0.6) - spanwise_aileron_weight()*(aircraft.b/2*0.4) - spanwise_fuel_weight()*(aircraft.b/2) - spanwise_wing_weight()*(aircraft.b/2) - aircraft.W_wl*aircraft.g0)
print(Bz)

# first: make fuel and structure as function of span _> then Az 


# # integrate wing loading along span
# wing_loading = quad(spanwise_wing_loading_ty, 0, aircraft.b/2)
# print(wing_loading)



