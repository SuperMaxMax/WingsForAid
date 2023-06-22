import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import aerodynamics.wing_planform as wp
import aerodynamics.horizontal_tail_design as htd
import aerodynamics.vertical_tail_design as vtd
import aerodynamics.drag_estimation as de
import numpy as np
import matplotlib.pyplot as plt
import copy
import scipy.integrate as integrate

def run_aero(aircraft):
    plt.figure()
    tapers = np.arange(0.1,1.1,0.1)
    for taper in tapers:
        aircraft_temp = copy.deepcopy(aircraft)
        x, y = wp.main_wing_planform(aircraft_temp, taper)
        plt.plot(x, y, label = f"taper = {taper}")

    u = 0
    v = 0
    a = aircraft.b/2
    b = -integrate.simps(y,x)*4/np.pi/a
    t = np.linspace(0, 0.5*np.pi, 30)
    plt.plot(u+a*np.cos(t), v+b*np.sin(t), label = "Elliptical Lift Distribution", color = 'k', linestyle = 'dashed', linewidth=2)
    plt.legend()
    plt.show()

    aircraft.taper = float(input("Enter taper ratio: "))

    wp.main_wing_planform(aircraft, aircraft.taper)
    htd.horizontal_tail_planform(aircraft)
    vtd.horizontal_tail_planform(aircraft)
    de.zero_liftdrag(aircraft)
