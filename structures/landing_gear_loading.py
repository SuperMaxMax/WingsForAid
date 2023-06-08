import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
from fuselage_truss_loading import dF_v_gust, dF_h_gust, yield_stress, E, density, dF_h_man

z_cg =