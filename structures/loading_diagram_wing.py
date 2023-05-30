import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad, cumulative_trapezoid, romberg, quadrature
from scipy.optimize import minimize
from math import tan, atan, cos
from material_dictionary import material, material_df, rank_material

def plot_loading_diagram(Ay, Az, By, Bz, load_case):
    plt.figure()
    y = np.linspace(0, ac.b/2, 1000)
    # add all forces with hatch pattern
    plt.plot(y, [-spanwise_wing_weight(y)]*1000, label='Wing weight')
    plt.plot(y, -spanwise_flap_weight(y), label='Flap weight')
    plt.plot(y, -spanwise_aileron_weight(y), label='Aileron weight')
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
        plt.arrow(ac.wing_strut_location-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
    else:
        plt.arrow(-np.sign(Ay)*0.5, 0, np.sign(Ay)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.arrow(ac.wing_strut_location-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
    if np.sign(Az) == 1:
        plt.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    else:
        plt.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    if np.sign(Bz) == 1:
        plt.arrow(ac.wing_strut_location, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    else:
        plt.arrow(ac.wing_strut_location, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
    plt.xlabel('Spanwise location [m]')
    plt.ylabel('Force [N/m]')
    plt.legend()
    plt.grid()
    plt.show()

def spanwise_wing_weight(y):
    # kg / m
    return (ac.W_w/ac.b)*ac.g0

def spanwise_wing_weight_ty(y):
    # kg / m
    return spanwise_wing_weight(y)*y

def spanwise_wing_loading(y):
    # put if for worst case scenarios
    gamma_0 = (ac.W_TO-ac.W_F)*ac.g0*(4/np.pi)*(1/ac.rho_cruise)*(1/ac.V_cruise)*(1/ac.b)
    return ac.rho_cruise*ac.V_cruise*gamma_0*np.sqrt(1-(2*y/ac.b)**2)

def spanwise_wing_loading_ty(y):
    # put if for worst case scenarios
    return spanwise_wing_loading(y)*y

def spanwise_flap_weight():
    # kg / m
    spanwise_weight = (ac.W_sc/2)*0.6/((ac.flap_e-ac.flap_s)*(ac.b/2))*ac.g0
    y_span = np.linspace(0, ac.b/2, 1000)
    test = spanwise_weight*(np.heaviside(y_span-ac.flap_s*ac.b/2, 1)-np.heaviside(y_span-ac.flap_e*ac.b/2,1))
    return test # spanwise_weight*(np.heaviside(y-ac.flap_s, 1)-np.heaviside(y-ac.flap_e,1))

def spanwise_flap_weight_ty(y):
    # kg / m
    return spanwise_flap_weight(y)*y

def spanwise_aileron_weight(y):
    # kg / m
    spanwise_weight = (ac.W_sc/2)*0.4/((ac.ail_e-ac.ail_s)*(ac.b/2))*ac.g0
    return spanwise_weight*(np.heaviside(y-ac.ail_s*ac.b/2, 1)-np.heaviside(y-ac.ail_e*ac.b/2,1))

def spanwise_aileron_weight_ty(y):
    # kg / m
    return spanwise_aileron_weight(y)*y

def spanwise_fuel_weight(y):
    # kg / m
    return ac.W_F/ac.b*ac.g0

def spanwise_fuel_weight_ty(y):
    # kg / m
    return spanwise_fuel_weight(y)*y

def total_spanwise(y):
    return -spanwise_wing_weight(y)+spanwise_wing_loading(y)-spanwise_flap_weight(y)-spanwise_aileron_weight(y)-spanwise_fuel_weight(y)

ac = UAV('aircraft')

# Strut_location
ac.wing_strut_location = 0.437277492*ac.b/2                       # [m] | based on statistical data, value is now based on half span running from middle fuselage
ac.strut_angle = atan(1.1/(ac.wing_strut_location-ac.w_out/2))    # [rad]
ac.l_strut = ac.wing_strut_location/np.cos(ac.strut_angle)        # [m]
ac.W_wl = 0.5                                                     # [kg] for 1 winglet
ac.flap_s = 0.2 # fraction of b/2
ac.flap_e = 0.6 # fraction of b/2
ac.ail_s = 0.6  # fraction of b/2
ac.ail_e = 0.8  # fraction of b/2
plot = True

# Material choice for strut weight
# add column for material index: sqrt(E)/rho
material_df['sqrt(E)/rho'] = np.sqrt(material_df['E'])/material_df['density']
# density, raw cost, eco cost, co2, yield stress, E, Kc, sqrt(E)/rho
prop_weights_strut = [0, 0.4, 0.1, 0.05, 0, 0, 0.05, 0.4]
# False -> the higher the value the better
possible_materials = rank_material(prop_weights_strut, [False])[:5]
ms = pd.DataFrame(index=possible_materials, columns=['Case 1', 'Case 2'])

for i in range(2):
    if i == 0:
        load_case = 1
    else:
        load_case = 0

    # Calculate reaction forces # CHECK SIGNS
    print(cumulative_trapezoid(spanwise_flap_weight(), np.linspace(0, ac.b/2, 1000), initial=0))
    quit()
    Bz = -1/ac.wing_strut_location*(quad(spanwise_wing_loading_ty, 0, ac.b/2)[0]*load_case-quad(spanwise_flap_weight_ty, 0, ac.b/2)[0]-quad(spanwise_aileron_weight_ty, 0, ac.b/2)[0]-quad(spanwise_fuel_weight_ty, 0, ac.b/2)[0]*abs(load_case-1)-quad(spanwise_wing_weight_ty, 0 , ac.b/2)[0]-ac.W_wl*ac.g0*ac.b/2)
    Az = (Bz+quad(spanwise_wing_loading, 0, ac.b/2)[0]*load_case - quad(spanwise_flap_weight, 0, ac.b/2)[0] - quad(spanwise_aileron_weight, 0, ac.b/2)[0] - quad(spanwise_fuel_weight, 0, ac.b/2)[0]*abs(load_case-1) - quad(spanwise_wing_weight, 0 ,ac.b/2)[0] - ac.W_wl*ac.g0)
    By = Bz/tan(ac.strut_angle)
    Ay = -By
    
    P_truss = np.sqrt(By**2+Bz**2)
    P_truss *= ac.ST_SF

    # Plot force diagram
    if plot:
        plot_loading_diagram(Az, Bz, By, Ay, load_case)

        # Make diagrams showing shear and bending moment
        # shear diagram
        y_span = np.linspace(0, ac.b/2, 200)
        y_strut = min(y_span, key=lambda x:abs(x-ac.wing_strut_location))
        shear = []
        bending = []
        for i in y_span:
            shear_y = quad(total_spanwise, i, ac.b/2)[0] + Bz*(1-np.heaviside(i-ac.wing_strut_location, 1)) # + smth for reaction force? # CHECK SIGNS
            shear.append(shear_y)
        # get the bending moment by integrating the shear along the length
        moment = cumulative_trapezoid(shear, y_span, initial=0)
        plt.figure()
        plt.subplot(211)
        plt.plot(y_span, shear)
        plt.xlabel('Spanwise location [m]')
        plt.ylabel('Shear force [N]')
        plt.subplot(212)
        plt.plot(y_span, moment)
        plt.xlabel('Spanwise location [m]')
        plt.ylabel('Bending moment [Nm]')
        plt.grid()
        plt.show()


    for mat in possible_materials:
        sigma_yield = material[mat]['yield stress']*10**6   # [Pa]
        mat_density = material[mat]['density']              # [kg/m^3]
        E = material[mat]['E']*10**9                        # [Pa]

        a = 0.0125       # semi_minor of ellipse [m]
        b = 2*a          # semi_major of ellipse [m]

        if load_case == 1:
            A_min = P_truss / sigma_yield
            m = A_min*mat_density*ac.l_strut
            ms.loc[mat, 'Case 1'] = m
        else:
            I_min = P_truss * ac.l_strut**2 / (np.pi**2 * E)
            t = 4*I_min/(np.pi*a**3 * (1+(3*b/a)))  # thickness of strut sheet [m]
            A = np.pi*(a+b)*t
            I_check = np.pi*a**3*t*(1+(3*b/a))/4
            m = A*mat_density*ac.l_strut
            ms.loc[mat, 'Case 2'] = m

print(ms)

# tomorrow
# correct flap and aileron
# check if other formulas are still correct




