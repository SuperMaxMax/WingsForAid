import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from control_stability.Class_II_cg_estimation import cg_calc
from control_stability.horizontal_tailplane import hor_run
from control_stability.vertical_tailplane import ver_run
from control_stability.control_surfaces import main_control_surface
from control_stability.wing_planform import WingPlanform

# from parameters import UAV, atmosphere
# aircraft = UAV('aircraft')
# atm      = atmosphere()

def main_stab_control(aircraft, plot, print_value):
    cg_calc(aircraft, plot)
    hor_run(aircraft, plot)
    ver_run(aircraft)
    exec(open("control_stability/landing_gear.py").read(), {'aircraft':aircraft, 'plot':plot})
    main_control_surface(aircraft)
    WingPlanform(aircraft, plot)

    if print_value:
        print("\n-----------------Summary-Stabily-&-Control----------------")

        print("-------------------C.G. Estimation---------------------------")
        print(f"CG_range is from {round(aircraft.X_cg_fwd, 3)} to {round(aircraft.X_cg_aft, 3)}")
        print(f"X_LEMAC is {round(aircraft.X_LEMAC, 3)} [m]\n")
        print(f"l_h:{aircraft.l_h}")

        print("--------------------Horizontal tailplane----------------------")
        print(f"Sh_S is {round(aircraft.Sh_S, 3)}\n")
        print(f"L_h is {round(aircraft.L_h, 3)}N")

        print("--------------------Vertical tailplane-------------------------")
        print(f"Sv_V is {round(aircraft.Sv_S, 3)}\n")
        print(f"\nC_m_alpha is {round(aircraft.C_m_alpha, 3)} [1/rad]\n")

        print("--------------------Landing Gear Design------------------------")

        print("--------------------Control Surfaces----------------------------")

        print("\n-------------------Aileron Design-----------------------------")
        print(f"The aileron will span from {round(aircraft.ystart_ail, 3)}m to {round(aircraft.yend_ail, 3)}[m] of the span of the wing with respect to the rootchord")
        print(f"The area of the aileron on one wing will be {round(aircraft.S_aileron, 3)} [m^2]\n")

        print("--------------------Elevator Design------------------------------")

        print("---------------------Rudder Design--------------------------------")
        print(f"The maximum required rudder deflection will be {round(aircraft.delta_r_value, 3)} [deg]")
# main_stab_control(aircraft, True, True)
