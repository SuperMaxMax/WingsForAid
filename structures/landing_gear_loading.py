import numpy as np
import matplotlib.pyplot as plt
from parameters import UAV
ac=UAV("aircraft")
from parameters import atmosphere
at=atmosphere()
from fuselage_truss_loading import dF_v_gust, dF_h_gust, yield_stress, E, density, dF_h_man

masses = np.array([12*ac.boxweight, ac.ST_W_eng, ac.ST_W_tb, ac.ST_W_fus, ac.W_w, ac.W_F])
z_cg = np.array([0.3,0.5061,0.67,0.67/2,0.67,0.67])

z_cg_tot = sum(z_cg*masses)/sum(masses)
print(z_cg_tot)
z_cg_ground = z_cg_tot + ac.ST_z_ground
print(z_cg_ground)
