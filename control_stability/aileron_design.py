import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
from scipy.integrate import quad

sys.path.append('.')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from horizontal_tailplane import hor_run
hor_run(aircraft)

def aileron_design(aircraft):
    t_lim = 2.6
    t = 2.7
    p =0.94

    while t > t_lim:
        y0 = 0.95 * aircraft.b/2    #distance symmetry line aircraft to start aileron
        y1 = p * aircraft.b/2       #distance symmetry line aircraft to end aileron

        tau = 0.41          # factor based on the fraction of chord that is the aileron

        C_l_delta_a = (2*aircraft.CLa_w_cruise* tau *aircraft.rootchord/(aircraft.Sw * aircraft.b)) * ((y0**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * y0**3)- (y1**2/2 + 2/3 * ((aircraft.taper - 1)/aircraft.b) * y1**3))
        
        delta_a_max = 20 * (np.pi/180)

        C_l = C_l_delta_a * delta_a_max

        L_a = 1/2 * atm.rho0 * (aircraft.V_s_min * 1.3)**2 * aircraft.Sw * aircraft.b * C_l
        
        y_d = 0.4 * aircraft.b/2        # [m] This is the average distance between the centre of gravity of the aircraft and the centre of drag
        Cdr = 0.9
        Pss = np.sqrt(2*L_a/((aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw) * Cdr * y_d**3))
    
        Ixx = 10000      # This is a total random guess and thus shall be changed. 

        phi_1 = Ixx/(atm.rho0 * y_d**3* Cdr * (aircraft.Sw + aircraft.AE_Sh_S * aircraft.Sw + aircraft.AE_Sv_S * aircraft.Sw)) * np.log(Pss**2)

        P_roll_rate = Pss**2/(2*phi_1)

        phi_des = 30 * (np.pi/180)

        t = np.sqrt(2*phi_des/P_roll_rate)
        print(t)
        p = p - 0.01


    print("\n--------Aileron Design------------\n")
    print(f"The aileron will span from {p*100}% to 95% of the span of the wing")
    print(f"The time to roll 30 degrees will be {t} s")

    aircraft.y_a_0 = p
    aircraft.y_a_1 = 0.95

def elevator_design(aircraft):
    pass


def rudder_design(aircraft):
    pass

aileron_design(aircraft)