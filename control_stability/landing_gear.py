import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
# import sympy as sp 


sys.path.append('.')

#from parameters import UAV, atmosphere
#aircraft = UAV('aircraft')
#atm = atmosphere()

# Position x_cg forward and backward. 
aircraft.x_cg_position_aft = aircraft.X_LEMAC + aircraft.X_cg_aft * aircraft.MAC_length
aircraft.x_cg_position_fwd = aircraft.X_LEMAC + aircraft.X_cg_fwd * aircraft.MAC_length

# initial guess for position. Is needed to run, but does not influence result. 
# aircraft.position_landing_fwd = [2, 0]
# aircraft.position_landing_back = [4, 1]

# Initial values for plotting
left_click_point = None
right_click_point = None

line_1 = None
line_2 = None
line_3 = None
line_4 = None
line_5 = None
symmetry_point = None
Z_position_cg = 0

def longitudinal_position_landing_gear(aircraft, x_point, y_point):
    global line_1
    global line_2
    global symmetry_point

    if line_1 is not None:
        line_1.remove()
    if line_2 is not None:
        line_2.remove()
    if symmetry_point is not None:
        symmetry_point.remove()

    weight_front_max = 0.15
    weight_front_min = 0.08

    symmetry_point = ax.scatter([x_point], [-y_point], s = 100, color='green', marker='x')

    def calc_x_front_max_and_min(x_cg_position):
        length_b_m = x_point - x_cg_position
        d_max = length_b_m / weight_front_min
        d_min = length_b_m / weight_front_max
        x_front_min = x_point - d_max  # Smallest distance from the nose
        x_front_max = x_point - d_min  # Largest distance from the nose
        return x_front_max, x_front_min

    x_front_max_aft, x_front_min_aft = calc_x_front_max_and_min(aircraft.x_cg_position_aft)
    x_front_max_fwd, x_front_min_fwd = calc_x_front_max_and_min(aircraft.x_cg_position_fwd)

    x_front_min = max(x_front_min_aft, x_front_min_fwd)
    x_front_max = min(x_front_max_aft, x_front_max_fwd)

    # Plot limit for the front undercarriage
    line_1 = ax.plot([x_front_max, x_front_max], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=75, length=0.7)])[0]
    line_2 = ax.plot([x_front_min, x_front_min], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])[0]
    return line_1, line_2, symmetry_point
    

def lateral_position_landing_gear(aircraft):
    "lateral limits on the aircraft"
    global line_3
    global line_4
    global line_5
    global Z_position_cg

    if line_3 is not None:
        line_3.remove()
    if line_4 is not None:
        line_4[0].remove()
    if line_5 is not None:
        line_5[0].remove()

    # Turnover angle limit of z-position of c.g.

    d = aircraft.position_landing_back[0] - aircraft.position_landing_fwd[0]

    angle_alpha = np.arctan(aircraft.position_landing_back[1] / d)

    c = (aircraft.x_cg_position_aft - aircraft.position_landing_fwd[0]) * np.sin(angle_alpha)
    angle_psi = 55 * (np.pi / 180)  # rad
    Z_position_cg = np.tan(angle_psi) * c
    parameter_text.set_text(f'Max Z - position C.G.: {Z_position_cg:.2f} [m]')

    
    # Pitch Angle limit NOTE: 15 degrees can change for our aircraft. 
    angle_theta_max = np.arctan(aircraft.h_out/(aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom-aircraft.position_landing_back[0]))
    # print('HHHHHHH',angle_theta_max * 180 / np.pi)
    # angle_theta_max = 15 * (np.pi/180) 
    x_cg_aft_limit = aircraft.x_cg_position_aft + aircraft.ST_z_cg_ground * np.tan(angle_theta_max)
    line_3 = ax.plot([x_cg_aft_limit, x_cg_aft_limit], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=75, length=0.7)])[0]

    # II - Limit for the main leg position for give N to attain stability against turnover NOTE: right now z_cg is assumed to be 0.8*z_cg_max. This should be an arbitrary value. 
    d_wheel = 0.3
    n_y = 0.5
    k_sg = 1
    e_s = 1/3 * d_wheel
    track_width = 2 * aircraft.position_landing_back[1]

    r = n_y * aircraft.ST_z_cg_ground * (1 + 4 * k_sg * e_s * aircraft.ST_z_cg_ground/ track_width)

    b = np.sqrt((aircraft.position_landing_fwd[0]-aircraft.x_cg_position_fwd)**2+(aircraft.position_landing_fwd[1]-0)**2)
    angle_theta = np.arccos(r/b)
    dir_angle = np.arctan2((aircraft.position_landing_fwd[1]-0), (aircraft.position_landing_fwd[0]-aircraft.x_cg_position_fwd))
    dir_angle1 = dir_angle + angle_theta
    dir_angle2 = dir_angle - angle_theta

    T1x = aircraft.x_cg_position_fwd + r *np.cos(dir_angle1)
    T1y = 0 + r *np.sin(dir_angle1)

    m = (T1y - aircraft.position_landing_fwd[1])/(T1x - aircraft.position_landing_fwd[0])

    y_1_range = []
    y_2_range = []
    
    x_range = np.arange(aircraft.position_landing_back[0]-0.3, aircraft.position_landing_back[0]+0.3, 0.02)
    for i in x_range:
        y_1 = m * i - m * aircraft.position_landing_fwd[0]
        y_2 = -m * i + m * aircraft.position_landing_fwd[0]
        y_1_range.append(y_1)
        y_2_range.append(y_2)

    # plotting turn_over limit II in torenbeek. 
    line_4 = ax.plot(x_range, y_1_range, linewidth = '0.8', color='red', path_effects=[patheffects.withTickedStroke(spacing=5, angle=75, length=0.7)])
    line_5 = ax.plot(x_range, y_2_range, linewidth = '0.8', color='red', path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])
    # ax.text(x_range[-1]+0.1, y_1_range[-1]+0.1, "II")
    # ax.text(x_range[-1]+0.1, y_2_range[-1]+0.1, "II")

    plt.draw()
    return line_3, line_4, line_5

##########################################################################################################
def on_click(event):
    """Clicking of position of forward and backward undercarriage"""
    global left_click_point
    global right_click_point

    if event.inaxes is not None:
        if event.button == 1:  # Left-click
            if left_click_point is None:
                # Add a new left-click point
                left_click_point = ax.scatter([event.xdata], [event.ydata], s = 100, color='green', marker='x')
                
            else:
                # Update the position of the left-click point
                left_click_point.set_offsets([[event.xdata, event.ydata]])
            aircraft.position_landing_back = [event.xdata, event.ydata]
            longitudinal_position_landing_gear(aircraft, event.xdata, event.ydata)
            lateral_position_landing_gear(aircraft)

        elif event.button == 3:  # Right-click
            if right_click_point is None:
                # Add a new right-click point
                right_click_point = ax.scatter([event.xdata], [0], s = 50, color='green', marker='^')
                
            else:
                # Update the position of the right-click point
                right_click_point.set_offsets([[event.xdata, 0]])

            # Call the lateral_position_landing_gear function to update Z_position_cg
            aircraft.position_landing_fwd = [event.xdata, 0]

            lateral_position_landing_gear(aircraft)
        # Redraw the plot
        plt.draw()

def on_key(event):
    if event.key == 'enter':
        # Print the values
        print(f"\nNosewheel position:\nX-position: {round(aircraft.position_landing_fwd[0],3)} [m]")
        print(f"\nMain landing gear position:\nX-position: {round(aircraft.position_landing_back[0],3)} [m],\
        Y-position:{round(aircraft.position_landing_back[1],3)} [m]")
        print(f"\nMax CG position Z-axis: {round(Z_position_cg,3)}")       

##########################################################################################################
"""Plotting"""
# Enable interactive mode
plt.ion()

# Create a figure and an axis
fig, ax = plt.subplots(figsize=(8, 8))

# Create a scatter plot with a single point
# point = ax.scatter([0], [0])

# Connect the 'button_press_event' to the 'on_click' function
cid = fig.canvas.mpl_connect('button_press_event', on_click)
cid = fig.canvas.mpl_connect('key_press_event', on_key)


# -- Plotting fuselage --
# Fuselage wall
ax.plot([0.4, aircraft.l_f-aircraft.l_fus_tail_cone], [aircraft.w_out/2, aircraft.w_out/2], color='0.25', linewidth=0.8)
ax.plot([0.4, aircraft.l_f-aircraft.l_fus_tail_cone], [-aircraft.w_out/2, -aircraft.w_out/2], color='0.25', linewidth=0.8)

# Fuselage nose
ax.plot([0, 0.1], [0, aircraft.w_out/4], color='0.25', linewidth=0.8)
ax.plot([0, 0.1], [0, -aircraft.w_out/4], color='0.25', linewidth=0.8)

ax.plot([0.4, 0.1], [aircraft.w_out/2, aircraft.w_out/4], color='0.25', linewidth=0.8)
ax.plot([0.4, 0.1], [-aircraft.w_out/2, -aircraft.w_out/4], color='0.25', linewidth=0.8)

# Fuselage tailcone
ax.plot([aircraft.l_f-aircraft.l_fus_tail_cone, aircraft.l_f], [aircraft.w_out/2, 0], color='0.25', linewidth=0.8)
ax.plot([aircraft.l_f-aircraft.l_fus_tail_cone, aircraft.l_f], [-aircraft.w_out/2, 0], color='0.25', linewidth=0.8)

# Boom
ax.plot([aircraft.l_f - aircraft.l_fus_tail_cone, aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom], [0.1, 0.1], color='0.25', linewidth=0.8)
ax.plot([aircraft.l_f - aircraft.l_fus_tail_cone, aircraft.l_f - aircraft.l_fus_tail_cone + aircraft.l_f_boom], [-0.1, -0.1], color='0.25', linewidth=0.8)

ax.set_xlim(-0.3, aircraft.l_f+4)
ax.set_ylim(-1, aircraft.b/2+1)
ax.axhline(0, color='blue', linestyle='dashdot', linewidth = 1)

# Plotting z-cg. max location in top right corner
parameter_text = ax.text(0.45, 0.98, f'Max Z - position C.G.: {Z_position_cg:.2f} [m]',
                        transform=ax.transAxes, ha='right', va='top')

plt.text(aircraft.x_cg_position_fwd - 0.07, -0.17, r'$x_{cg, fwd}$', fontsize=8, va='center', ha = 'left')
plt.text(aircraft.x_cg_position_aft + 0.05, -0.17, r'$x_{cg, aft}$', fontsize=8, va='center', ha = 'left')

# Plot Wing Position and planform
ax.plot([aircraft.X_LEMAC, aircraft.X_LEMAC+aircraft.MAC_length], [aircraft.y_mac + aircraft.w_out/2, aircraft.y_mac + aircraft.w_out/2],color='0.25', linewidth=0.8)

ax.plot([aircraft.X_LEMAC - aircraft.x_lemac, aircraft.X_LEMAC - aircraft.x_lemac+aircraft.rootchord], [aircraft.w_out/2, aircraft.w_out/2], color='0.25', linewidth=0.8)

x_root_quart = aircraft.X_LEMAC-aircraft.x_lemac + 1/4*aircraft.rootchord
ax.plot([x_root_quart-1/4*aircraft.tipchord, x_root_quart+3/4*aircraft.tipchord], [aircraft.b/2 + aircraft.w_out/2, aircraft.b/2 + aircraft.w_out/2], color='0.25', linewidth=0.8)
ax.plot([aircraft.X_LEMAC-aircraft.x_lemac, x_root_quart-1/4*aircraft.tipchord], [aircraft.w_out/2, aircraft.b/2 + aircraft.w_out/2], color='0.25', linewidth=0.8)
ax.plot([aircraft.X_LEMAC-aircraft.x_lemac+aircraft.rootchord, x_root_quart+3/4*aircraft.tipchord], [aircraft.w_out/2, aircraft.b/2 + aircraft.w_out/2], color='0.25', linewidth=0.8)

# Plotting CG
ax.plot([aircraft.x_cg_position_fwd, aircraft.x_cg_position_aft], [0, 0], \
        color = 'black', linewidth = 2, marker = '|', markersize = 14, label = 'C.G. Range')
 
ax.set_aspect('equal')

# plt.grid()
plt.legend(loc="upper right")
plt.xlabel("X-position [m]")
plt.ylabel("Y-position [m]")


# Show the plot
if plot:
    plt.show(block=True)
