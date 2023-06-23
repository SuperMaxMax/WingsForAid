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
import scipy.interpolate as interpolate

def run_aero(aircraft):
    tapers = np.arange(0.4,0.8,0.01)
    diffs = []
    plt.figure()
    for i, taper in enumerate(tapers):
        aircraft_temp = copy.deepcopy(aircraft)
        x, y = wp.main_wing_planform(aircraft_temp, taper)
        plt.plot(x, y, label = f"taper = {taper}")

        span = np.linspace(0, aircraft_temp.b/2, 100)

        if i == 0:
            u = 0
            v = 0
            a = aircraft.b/2
            b = -integrate.simps(y,x)*4/np.pi/a
            t = np.linspace(0, 0.5*np.pi, 30)
            x_ell = u+a*np.cos(t)
            y_ell = v+b*np.sin(t)
            ell = interpolate.interp1d(x_ell, y_ell, kind='cubic',fill_value="extrapolate")
        
        tap = interpolate.interp1d(x, y, kind='cubic',fill_value="extrapolate")
        # calculate difference between ell and tap
        diffs.append(sum(abs(ell(span) - tap(span))))
    plt.plot(x_ell, y_ell, label = "Elliptical Lift Distribution", color = 'k', linestyle = 'dashed', linewidth=2)
    plt.legend()
    # plt.show()

    aircraft.taper = tapers[np.argmin(diffs)]
    print('chosen optimum taper = ', aircraft.taper)

    wp.main_wing_planform(aircraft, aircraft.taper)
    htd.horizontal_tail_planform(aircraft)
    vtd.horizontal_tail_planform(aircraft)
    de.zero_liftdrag(aircraft)