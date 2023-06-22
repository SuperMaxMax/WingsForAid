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

def run_aero(aircraft):
    print(aircraft.tipchord, 'taper before')
    plt.figure()
    tapers = np.arange(0.1,1.1,0.1)
    for taper in tapers:
        aircraft_temp = copy.deepcopy(aircraft)
        print(aircraft.tipchord, 'taper of temp aircraft')
        x, y = wp.main_wing_planform(aircraft_temp, taper)
        print(aircraft.tipchord, 'taper of temp aircraft after')
        plt.plot(x, y, label = f"taper = {taper}")
    plt.legend()
    plt.show()

    aircraft.taper = float(input("Enter taper ratio: "))

    wp.main_wing_planform(aircraft, aircraft.taper)
    htd.horizontal_tail_planform(aircraft)
    vtd.horizontal_tail_planform(aircraft)
    de.zero_liftdrag(aircraft)
