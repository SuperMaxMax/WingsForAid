import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from parameters import UAV

import wing_planform as wp
import horizontal_tail_design as htd
import vertical_tail_design as vtd

aircraft = UAV('aircraft')

wp.main_wing_planform(aircraft)
htd.horizontal_tail_planform(aircraft)
vtd.horizontal_tail_planform(aircraft)