import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import aerodynamics.wing_planform as wp
import aerodynamics.horizontal_tail_design as htd
import aerodynamics.vertical_tail_design as vtd

def run_aero(aircraft):
    wp.main_wing_planform(aircraft)
    htd.horizontal_tail_planform(aircraft)
    vtd.horizontal_tail_planform(aircraft)

