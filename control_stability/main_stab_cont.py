import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects

sys.path.append('..')

from parameters import UAV, atmosphere
aircraft = UAV('aircraft')
atm      = atmosphere()

from Class_II_cg_estimation import cg_calc
from horizontal_tailplane import hor_run
from vertical_tailplane import ver_run
from control_surfaces import main_control_surface

def main_stab_control(aircraft, plot, print_value):

    # cg_calc(aircraft, plot)
    hor_run(aircraft, plot)
    # ver_run(aircraft)
    # exec(open("landing_gear.py").read(), {'aircraft':aircraft, 'plot':plot})
    # main_control_surface(aircraft)

    if print_value:
        print("\n-----------------Summary-Stabily-&-Control----------------\n")

        print("-------------------C.G. Estimation---------------------------\n")

        print(f"CG_range is from {round(aircraft.X_cg_fwd, 3)} to {round(aircraft.X_cg_aft, 3)}")
        print(f"X_LEMAC is {round(aircraft.X_LEMAC, 3)} m\n")

        print("--------------------Horizontal tailplane----------------------\n")

        print(f"Sh_S is {round(aircraft.AE_Sh_S, 3)}\n")

        print("--------------------Vertical tailplane-------------------------\n")

        print(f"Sv_V is {round(aircraft.AE_Sv_S, 3)}\n")
        print(f"\nC_m_alpha is {round(aircraft.C_m_alpha, 3)} [1/rad]\n")

        print("--------------------Landing Gear Design------------------------\n")

        print("--------------------Control Surfaces----------------------------")

        print("\n-------------------Aileron Design-----------------------------")
        print(f"The aileron will span from {round(aircraft.y_a_0 - aircraft.w_out/2, 3)}m to {round(aircraft.y_a_1 - aircraft.w_out/2, 3)}m of the span of the wing with respect to the rootchord")
        print(f"The area of the aileron on one wing will be {aircraft.S_aileron} m^2\n")

        print("--------------------Elevator Design------------------------------\n")

        print("---------------------Rudder Design--------------------------------\n")
        print(f"The maximum required rudder deflection will be {round(aircraft.delta_r_value*180/np.pi, 3)}")
main_stab_control(aircraft, True, True)
