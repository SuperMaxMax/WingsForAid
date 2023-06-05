# 24. Select the vertical tail volume coefficient, V V (Table 6.4).
# 25. Assume the vertical tail moment arm (lv) equal to the horizontal tail moment arm (l).
# 26. Calculate vertical tail planform area, Sv (Equation (6.74)).
# 27. Select vertical tail airfoil section (Section 6.8.2.4).
# 28. Select vertical tail aspect ratio, ARv(Section 6.8.2.6).
# 29. Select vertical tail taper ratio, λv(Section 6.8.2.7).
# 30. Determine the vertical tail incidence angle (Section 6.8.2.5).
# 31. Determine the vertical tail sweep angle (Section 6.8.2.8).
# 32. Determine the vertical tail dihedral angle (Section 6.8.2.9).
# 33. Calculate vertical tail span (bv), root chord (Cvroot), and tip chord (Cvtip ), and
# MACv(Equations (6.76)–(6.79)).
# 34. Check the spin recovery.
# 35. Adjust the location of the vertical tail relative to the horizontal tail by changing lv to
# satisfy the spin recovery requirements (Section 6.8.2.2).
# 36. Analyze directional trim (Section 6.8.1).
# 37. Analyze directional stability (Section 6.8.1).
# 38. Modify to meet the design requirements.
# 39. Optimize the tail.

import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from parameters import UAV
object = UAV('aircraft')

def vertical_wing_design(aircraft):
    #airfoil data
    #NACA 0009
    #data from https://digital.library.unt.edu/ark:/67531/metadc65459/m2/1/high_res_d/19930090937.pdf except alpha stall
    list_Cl0 = 0
    list_Cd_min = 0.0064
    list_Cm_0 = 0
    list_alpha_0 = 0
    list_alpha_s = 13.2
    list_Cl_max = 1.39
    list_Cl_alpha = 0.098 * 180/np.pi
    list_tc = 0.09
    
    Sv_Sw = aircraft.AE_Sv_S
    b = aircraft.AE_b
    l_v = aircraft.AE_l_h
    Sw = aircraft.AE_Sw

    V_v = Sv_Sw * l_v / b
    S_v = Sv_Sw * Sw
    print(V_v)


vertical_wing_design(object)

