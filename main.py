from parameters import UAV
import aerodynamics.main_aero as ae
import 

aircraft = UAV("aircraft")

ae.wp.main_wing_planform(aircraft)
ae.htd.horizontal_tail_planform(aircraft)
ae.vtd.horizontal_tail_planform(aircraft)

cs.





