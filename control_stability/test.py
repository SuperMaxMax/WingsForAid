# import sympy as sp
# import numpy as np
# import matplotlib.pyplot as plt 
# from parameters import UAV
# import sys
# sys.path.append('..')
# aircraft = UAV("aircraft")

# # Definieer de symbolen x en y
# x, y = sp.symbols('x y')

# # Bepaal de coördinaten van punt B
# B = (1, 0)

# # Definieer de vergelijking van de cirkel
# circle_eq = (x - 4)**2 + (y)**2 - 13

# # Bepaal de afgeleide van de cirkelvergelijking
# dy_dx = sp.diff(circle_eq, y) / sp.diff(circle_eq, x)

# # Bepaal de helling van de raaklijnen in het punt B
# m = dy_dx.subs([(x, B[0]), (y, B[1])])
# print(m)
# # Bepaal de vergelijkingen van de raaklijnen n1 en n2
# y1 = - B[1] - m * (x - B[0])
# n2_eq = (y - B[1]) + m * (x - B[0])

# # Vereenvoudig de vergelijkingen
# #1_eq = sp.simplify(n1_eq)
# n2_eq = sp.simplify(n2_eq)
# y_1_range = []
# y_2_range = []

# x_range = np.arange(0, 5, 0.1)
# for i in x_range:
#     y_1 = B[1] - m * (i - B[0])
#     y_2 = B[1] + m * (i - B[0])
#     y_1_range.append(y_1)
#     y_2_range.append(y_2)

# plt.plot(x_range, y_1_range)
# plt.show()

# control_surfaces.py test bt Bram
# T= 2800 
# z_position_T = -(aircraft.ST_z_prop + aircraft.prop_radius - aircraft.ST_z_cg_ground)
# eta_h = 0.96
# V_range = np.arange(aircraft.V_s_min, aircraft.V_cruise*1.4, 2)

# color_r = ['black', 'steelblue']
# e = 0

# # Random constanten
# aircraft.CLa_Ah_cruise = 3
# aircraft.CLa_w_cruise = 3
# aircraft.C_m_alpha = 0.01

# for X_cg in [aircraft.X_cg_fwd, aircraft.X_cg_aft]:
#     delta_eq_req_range = []
#     for V in V_range:
#         l_h = aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom - 3/4 * aircraft.AE_rootchord_h - (aircraft.X_LEMAC+ X_cg*aircraft.MAC_length)
#         C_L1 = 2 * aircraft.W_TO * aircraft.g0 /(aircraft.rho_TO* V**2*aircraft.Sw)
#         tail_volume = l_h/aircraft.MAC_length * aircraft.Sh_S

#         C_L_delta_e = aircraft.CLa_Ah_cruise * eta_h * aircraft.Sh_S * 1
#         C_m_delta_e = -aircraft.CLa_Ah_cruise *eta_h * tail_volume * 1

#         delta_e_req = - ((T*z_position_T/(1/2 * aircraft.rho_TO * V**2 * aircraft.Sw * aircraft.MAC_length) + aircraft.af_cm0) * aircraft.CLa_w_cruise + (C_L1 - aircraft.af_Cl0)*aircraft.C_m_alpha)/(aircraft.CLa_w_cruise*C_m_delta_e - aircraft.C_m_alpha*C_L_delta_e)
#         delta_eq_req_range.append(delta_e_req*180/np.pi)
#     plt.scatter(V_range, delta_eq_req_range, marker='x', color=color_r[e])#(color_r[e] for e in range(len(color_r))))
#     e += 1

# plt.show()

# X_cg = np.array([aircraft.X_cg_fwd, aircraft.X_cg_aft])
# C_L_delta_e = aircraft.CLa_Ah_cruise * eta_h * aircraft.Sh_S * 1
# C_L1 = 2 * aircraft.W_TO * aircraft.g0 /(aircraft.rho_TO * (V_range**2)*aircraft.Sw)

# #Calculate for both cg positions, [fwd, aft]
# l_h = aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom - 3/4 * aircraft.AE_rootchord_h - (aircraft.X_LEMAC+ X_cg*aircraft.MAC_length)
# tail_volume = l_h / aircraft.MAC_length * aircraft.Sh_S
# C_m_delta_e = -aircraft.CLa_Ah_cruise *eta_h * tail_volume * 1

# delta_e_req_fwd = (-((T*z_position_T/(1/2 * aircraft.rho_TO * (V_range**2) * aircraft.Sw * aircraft.MAC_length) + aircraft.af_cm0) * aircraft.CLa_w_cruise + (C_L1 - aircraft.af_Cl0)*aircraft.C_m_alpha)/(aircraft.CLa_w_cruise*C_m_delta_e[0] - aircraft.C_m_alpha*C_L_delta_e)) * (180/np.pi)
# delta_e_req_aft = (-((T*z_position_T/(1/2 * aircraft.rho_TO * (V_range**2) * aircraft.Sw * aircraft.MAC_length) + aircraft.af_cm0) * aircraft.CLa_w_cruise + (C_L1 - aircraft.af_Cl0)*aircraft.C_m_alpha)/(aircraft.CLa_w_cruise*C_m_delta_e[1] - aircraft.C_m_alpha*C_L_delta_e)) * (180/np.pi)

# plt.scatter(V_range, delta_e_req_fwd, marker='x', color='black')
# plt.scatter(V_range, delta_e_req_aft, marker='x', color='steelblue')

# plt.show()
