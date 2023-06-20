import sys
sys.path.append("..")

#########################################################
# LOADING DIAGRAM WING
#########################################################
from parameters import UAV_final
import numpy as np
import matplotlib.gridspec as gridspec
from scipy.integrate import quad, cumulative_trapezoid, trapezoid
from scipy.optimize import minimize
from math import tan, atan, cos
import matplotlib.pyplot as plt

def plot_diagrams(Ay, Az, By, Bz, load_case, normal, shear, moment, torque):
    # Make diagram
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(4, 2)
    f_ax1 = fig.add_subplot(gs[:, 0])
    plt.rcParams.update({'font.size': 22})

    # Plot forces
    f_ax1.plot(y_span, -spanwise_wing_weight(y_span), label='Wing weight')
    f_ax1.plot(y_span, -spanwise_flap_weight(), label='Flap weight')
    f_ax1.plot(y_span, -spanwise_aileron_weight(), label='Aileron weight')
    f_ax1.plot(y_span, -spanwise_fuel_weight(y_span)*abs(load_case-1), label='Fuel weight')
    f_ax1.plot(y_span, spanwise_wing_loading(y_span), label='Wing loading')
    # set axis limits such that all vector arrows are visible, vector are scaled based on axis limits
    if load_case == 1:
        f_ax1.set_title('Spanwise force diagram for full wing loading, no fuel')
        global y_range1
        y_range1 = f_ax1.axis()[3]-f_ax1.axis()[2]
        hl1, hw1 = 0.06, 132
        hl2, hw2 = 132, 0.06
        w1, w2 = 33, 0.02
        global fac
        fac = 660
    else:
        f_ax1.set_title('Spanwise force diagram for no wing loading, full fuel')
        y_range2 = f_ax1.axis()[3]-f_ax1.axis()[2]
        hl1, hw1 = 0.06, 132*y_range2/y_range1
        hl2, hw2 = 132*y_range2/y_range1, 0.06
        w1, w2 = 33*y_range2/y_range1, 0.02
        fac *= y_range2/y_range1
    # add reaction forces
    if np.sign(Ay) == 1:
        f_ax1.arrow(-np.sign(Ay)*0.5, 0, np.sign(Ay)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.annotate(f'{round(abs(Ay),1)} N', xy=(-np.sign(Ay)*0.5, 0), xytext=(-np.sign(Ay)*0.5, hw1))
        f_ax1.arrow(ac.ST_y_strut-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.annotate(f'{round(abs(By),1)} N', xy=(ac.ST_y_strut-np.sign(By)*0.5, 0), xytext=(ac.ST_y_strut-np.sign(By)*0.5, hw1))
    else:
        f_ax1.arrow(-np.sign(Ay)*0.5, 0, np.sign(Ay)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.annotate(f'{round(abs(Ay),1)} N', xy=(-np.sign(Ay)*0.5, 0), xytext=(-np.sign(Ay)*0.5, hw1))
        f_ax1.arrow(ac.ST_y_strut-np.sign(By)*0.5, 0, np.sign(By)*0.5, 0, color='r', length_includes_head=True, width=w1, head_length=hl1, head_width=hw1)
        plt.annotate(f'{round(abs(By),1)} N', xy=(ac.ST_y_strut-np.sign(By)*0.5, 0), xytext=(ac.ST_y_strut-np.sign(By)*0.5, hw1))
    if np.sign(Az) == 1:
        f_ax1.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
        plt.annotate(f'{round(abs(Az),1)} N', xy=(0, -np.sign(Az)*fac), xytext=(hw2, -np.sign(Az)*fac))
    else:
        f_ax1.arrow(0, -np.sign(Az)*fac, 0, np.sign(Az)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
        plt.annotate(f'{round(abs(Az),1)} N', xy=(0, -np.sign(Az)*fac), xytext=(hw2, -np.sign(Az)*fac))
    if np.sign(Bz) == 1:
        f_ax1.arrow(ac.ST_y_strut, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
        plt.annotate(f'{round(abs(Bz),1)} N', xy=(ac.ST_y_strut, -np.sign(Bz)*fac), xytext=(ac.ST_y_strut+hw2, -np.sign(Bz)*fac))
    else:
        f_ax1.arrow(ac.ST_y_strut, -np.sign(Bz)*fac, 0, np.sign(Bz)*fac, color='g', length_includes_head=True, width=w2, head_length=hl2, head_width=hw2)
        plt.annotate(f'{round(abs(Bz),1)} N', xy=(ac.ST_y_strut, -np.sign(Bz)*fac), xytext=(ac.ST_y_strut+hw2, -np.sign(Bz)*fac))
    f_ax1.set_xlabel('Spanwise location [m]')
    f_ax1.set_ylabel('Force [N/m]')
    f_ax1.legend()
    f_ax1.grid()

    # Plot normal, shear, bending moment and torque
    f_ax2 = fig.add_subplot(gs[0, 1])
    f_ax2.plot(y_span, normal, 'g')
    f_ax2.fill_between(y_span, normal, 0, color='g', alpha=0.2, hatch='//')
    f_ax2.set_ylabel('Normal force [N]')
    f_ax2.grid()
    
    f_ax3 = fig.add_subplot(gs[1, 1])
    f_ax3.plot(y_span, shear, 'b')
    f_ax3.fill_between(y_span, shear, 0, color='b', alpha=0.2, hatch='//')
    f_ax3.set_ylabel('Shear force [N]')
    f_ax3.grid()

    f_ax4 = fig.add_subplot(gs[2, 1])
    f_ax4.plot(y_span, moment, 'r')
    f_ax4.fill_between(y_span, moment, 0, color='r', alpha=0.2, hatch='//')
    f_ax4.set_ylabel('Bending moment [Nm]')
    f_ax4.grid()

    f_ax4 = fig.add_subplot(gs[3, 1])
    f_ax4.plot(y_span, torque, color='purple')
    f_ax4.fill_between(y_span, torque, 0, color='purple', alpha=0.2, hatch='//')
    f_ax4.set_xlabel('Spanwise location [m]')
    f_ax4.set_ylabel('Torque [Nm]')
    f_ax4.grid()

    # Show
    plt.show()

def chord(y):
    return ac.rootchord-((ac.rootchord-ac.tipchord)/(ac.b/2))*y

def spanwise_wing_weight(y):
    # N/m
    avg_chord = ac.tipchord+(ac.rootchord-ac.tipchord)/2
    return ((ac.W_w/ac.b)*ac.g0)*(chord(y)/avg_chord)

def spanwise_wing_weight_ty(y):
    # N/m
    return spanwise_wing_weight(y)*y

def spanwise_wing_loading(y):
    # # put if for worst case scenarios
    # gamma_0 = (ac.W_TO-ac.W_F)*ac.g0*(4/np.pi)*(1/ac.rho_cruise)*(1/ac.V_cruise)*(1/ac.b)
    # return ac.rho_cruise*ac.V_cruise*gamma_0*np.sqrt(1-(2*y/ac.b)**2) # when inputting XFLR -> check for matching coordinate systems
    c = ac.lift_coefficients[::-1]
    return np.polyval(c, y)*n_ult

def spanwise_wing_loading_ty(y):
    # put if for worst case scenarios
    return spanwise_wing_loading(y)*y

def spanwise_flap_weight():
    # this function returns an array of spanwise flap weight
    spanwise_weight = (ac.W_sc/2)*0.6/((ac.yend_flap/(ac.b/2)-ac.ystart_flap/(ac.b/2))*(ac.b/2))*ac.g0 # N/m
    return spanwise_weight*(np.heaviside(y_span-ac.ystart_flap/(ac.b/2)*ac.b/2, 1)-np.heaviside(y_span-ac.yend_flap/(ac.b/2)*ac.b/2,1))

def spanwise_flap_weight_ty():
    # N/m 
    return np.multiply(spanwise_flap_weight(), y_span)

def spanwise_aileron_weight():
    # this function returns an array of spanwise flap weight
    spanwise_weight = (ac.W_sc/2)*0.4/((ac.yend_ail/(ac.b/2)-ac.ystart_ail/(ac.b/2))*(ac.b/2))*ac.g0 # N/m
    return spanwise_weight*(np.heaviside(y_span-ac.ystart_ail/(ac.b/2)*ac.b/2, 1)-np.heaviside(y_span-ac.yend_ail/(ac.b/2)*ac.b/2,1))

def spanwise_aileron_weight_ty():
    # N/m
    return np.multiply(spanwise_aileron_weight(), y_span)

def spanwise_fuel_weight(y):
    # N/m
    avg_chord = ac.tipchord+(ac.rootchord-ac.tipchord)/2
    return (ac.W_F/ac.b*ac.g0)*(chord(y)/avg_chord)

def spanwise_fuel_weight_ty(y):
    # N/m
    return spanwise_fuel_weight(y)*y

def spanwise_func_1(y):
    return -spanwise_wing_weight(y)+spanwise_wing_loading(y)

def spanwise_func_2(y):
    return -spanwise_wing_weight(y)-spanwise_fuel_weight(y)+spanwise_wing_loading(y)

ac = UAV_final()
plot = True

ac.ST_angle_strut = np.arctan(ac.ST_h_fus/(ac.ST_y_strut-ac.ST_w_fus/2))    # [rad]
ac.W_wl = 0                                                                 # [kg] for 1 winglet

# Make an array for the span
y_span = np.linspace(0, ac.b/2, 1000)

for k in range(3):
    if k == 0:
        load_case = 1
        n_ult = ac.ST_n_ult_pos
    elif k == 1:
        load_case = 0
        n_ult = ac.ST_n_ult_neg
    else:
        load_case = 1
        n_ult = 2

    # Calculate reaction forces # important -> signs of the forces are with z is positive upwards and y is positive in spanwise direction (this is different from global coordinate system)
    Bz = -1/ac.ST_y_strut*(quad(spanwise_wing_loading_ty, 0, ac.b/2)[0]-trapezoid(spanwise_flap_weight_ty(), y_span)-trapezoid(spanwise_aileron_weight_ty(), y_span)-quad(spanwise_fuel_weight_ty, 0, ac.b/2)[0]*abs(load_case-1)-quad(spanwise_wing_weight_ty, 0 , ac.b/2)[0]-ac.W_wl*ac.g0*ac.b/2)
    Az = -(Bz+quad(spanwise_wing_loading, 0, ac.b/2)[0] - trapezoid(spanwise_flap_weight(), y_span) - trapezoid(spanwise_aileron_weight(), y_span) - quad(spanwise_fuel_weight, 0, ac.b/2)[0]*abs(load_case-1) - quad(spanwise_wing_weight, 0 ,ac.b/2)[0] - ac.W_wl*ac.g0)
    By = Bz/tan(ac.ST_angle_strut)
    Ay = -By

    # Calculate normal, shear, moment and torque
    normal = []
    shear = []
    torque = []
    
    # flexural axis assumption # -> replace 0.25 with front spar location variable, also add in location based on wingbox design -> for iteration
    flex_ax = ac.x_strut

    # define moment arms
    cop = 0.25                      # chord from leading edge
    ma_wing_loading = (flex_ax-cop)*chord(y_span) # moment arm wing loading
    torque_wing_loading = -ma_wing_loading*spanwise_wing_loading(y_span) # torque wing loading

    pos_flap = ac.xc_aft_spar+0.4*(1-ac.xc_aft_spar)
    ma_flap = (flex_ax-pos_flap)*chord(y_span)
    torque_flap = -ma_flap*-spanwise_flap_weight()

    pos_ail = ac.xc_aft_spar+0.4*(1-ac.xc_aft_spar)
    ma_ail = (flex_ax-pos_ail)*chord(y_span)
    torque_ail = -ma_ail*-spanwise_aileron_weight()

    pos_strut = ac.xc_front_spar
    ma_strut = (flex_ax-pos_strut)*chord(ac.ST_y_strut)
    torque_strut = -ma_strut*Bz

    for j, i in enumerate(y_span):
        normal.append(By*(1-np.heaviside(i-ac.ST_y_strut, 1)))
        if load_case == 1:
            shear.append(-(quad(spanwise_func_1, i, ac.b/2)[0] - trapezoid(spanwise_flap_weight()[j:]+spanwise_aileron_weight()[j:], y_span[j:]) + Bz*(1-np.heaviside(i-ac.ST_y_strut, 1)) - ac.W_wl*ac.g0)) # + smth for reaction force? # CHECK SIGNS
            torque.append(trapezoid(torque_wing_loading[j:]+torque_flap[j:]+torque_ail[j:], y_span[j:]) + torque_strut*(1-np.heaviside(i-ac.ST_y_strut, 1)))
        else:
            shear.append(-(quad(spanwise_func_2, i, ac.b/2)[0] - trapezoid(spanwise_flap_weight()[j:]+spanwise_aileron_weight()[j:], y_span[j:]) + Bz*(1-np.heaviside(i-ac.ST_y_strut, 1)) - ac.W_wl*ac.g0)) # + smth for reaction force? # CHECK SIGNS
            torque.append(trapezoid(torque_wing_loading[j:]+torque_flap[j:]+torque_ail[j:], y_span[j:]) + torque_strut*(1-np.heaviside(i-ac.ST_y_strut, 1)))
    # get the bending moment by integrating the shear along the length
    moment = cumulative_trapezoid(shear, y_span, initial=0)

    # Plot diagrams
    if plot:
        plot_diagrams(Ay, Az, By, Bz, load_case, normal, shear, moment, torque)

    if k == 0:
        loading_tension = [normal, shear, moment, torque]
    elif k == 1:
        loading_compression = [normal, shear, moment, torque]
    else:
        loading_custom = [normal, shear, moment, torque]