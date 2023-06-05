import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patheffects
# import sympy as sp 

sys.path.append('.')

from parameters import UAV, atmosphere

aircraft = UAV('aircraft')
atm = atmosphere()

# Position x_cg forward and backward. 
aircraft.x_cg_position_aft = aircraft.X_LEMAC + aircraft.X_cg_aft * aircraft.MAC_length
aircraft.x_cg_position_fwd = aircraft.X_LEMAC + aircraft.X_cg_fwd * aircraft.MAC_length

# initial guess for position. Is needed to run, but does not influence result. 
aircraft.position_landing_fwd = [2, 0]
aircraft.position_landing_back = [4, 1]

# Initial values for plotting
line_1 = None
line_2 = None
line_3 = None
line_4 = None
line_5 = None

symmetry_point = None
left_click_point = None
right_click_point = None
Z_position_cg = 0

def longitudinal_position_landing_gear(aircraft, x_point, ypoint):
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

    symmetry_point = ax.scatter([x_point], [-ypoint], color='orange', marker='^')

    def calc_x_front_max_and_min(x_cg_position):
        length_b_m = x_point - x_cg_position
        d_max = length_b_m / weight_front_min
        d_min = length_b_m / weight_front_max
        x_front_min = x_point - d_max
        x_front_max = x_point - d_min
        return x_front_max, x_front_min

    x_front_max_aft, x_front_min_aft = calc_x_front_max_and_min(aircraft.x_cg_position_aft)
    x_front_max_fwd, x_front_min_fwd = calc_x_front_max_and_min(aircraft.x_cg_position_fwd)

    x_front_max = x_front_max_fwd
    x_front_min = x_front_min_aft

    # Plot limit for the front undercarriage
    line_1 = ax.plot([x_front_max, x_front_max], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=75, length=0.7)])[0]
    line_2 = ax.plot([x_front_min, x_front_min], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=-75, length=0.7)])[0]
    
    

def lateral_position_landing_gear(aircraft):
    "lateral limits on the aircraft"
    global Z_position_cg
    global line_3
    global line_4
    global line_5

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
    parameter_text.set_text(f'Max Z_position_c.g.: {Z_position_cg:.2f} [m]')

    
    # Pitch Angle limit NOTE: 15 degrees can change for our aircraft. 

    angle_theta_max = 30 * (np.pi/180) 
    x_cg_aft_limit = aircraft.x_cg_position_aft + Z_position_cg * np.tan(angle_theta_max)
    line_3 = ax.plot([x_cg_aft_limit, x_cg_aft_limit], [-1, 1], color='red', linewidth="0.8", path_effects=[patheffects.withTickedStroke(spacing=5, angle=75, length=0.7)])[0]

    # II - Limit for the main leg position for give N to attain stability against turnover NOTE: right now z_cg is assumed to be 0.8*z_cg_max. This should be an arbitrary value. 
    d_wheel = 0.3
    n_y = 0.5
    k_sg = 1
    e_s = 1/3 * d_wheel
    track_width = 2 * aircraft.position_landing_back[1]

    r = n_y * 0.8*Z_position_cg * (1 + 4 * k_sg * e_s * 0.8*Z_position_cg/ track_width)

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

    plt.draw()


##########################################################################################################

def on_click(event):
    """Clicking of position of forward and backward undercarriage"""
    global left_click_point
    global right_click_point

    if event.inaxes is not None:
        if event.button == 1:  # Left-click
            if left_click_point is None:
                # Add a new left-click point
                left_click_point = ax.scatter([event.xdata], [event.ydata], color='orange', marker='^')
                
            else:
                # Update the position of the left-click point
                left_click_point.set_offsets([[event.xdata, event.ydata]])
            aircraft.position_landing_back = [event.xdata, event.ydata]
            longitudinal_position_landing_gear(aircraft, event.xdata, event.ydata)
            lateral_position_landing_gear(aircraft)

        elif event.button == 3:  # Right-click
            if right_click_point is None:
                # Add a new right-click point
                right_click_point = ax.scatter([event.xdata], [0], color='orange', marker='^')
                
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
        print(f"Position Landing Fwd: {aircraft.position_landing_fwd}")
        print(f"Position Landing Back: {aircraft.position_landing_back}")
        print(f"Z_cg_position: {Z_position_cg}")       

##########################################################################################################
"""Plotting"""
# Enable interactive mode


plt.ion()

# Create a figure and an axis
fig, ax = plt.subplots(figsize=(8, 8))

# Create a scatter plot with a single point
point = ax.scatter([0], [0])

# Connect the 'button_press_event' to the 'on_click' function
cid = fig.canvas.mpl_connect('button_press_event', on_click)
cid = fig.canvas.mpl_connect('key_press_event', on_key)

# Set the plot limits
x_cg_point = ax.scatter([aircraft.x_cg_position_aft], [0], color='black', label='2: c.g. aft')
x_cg_point = ax.scatter([aircraft.x_cg_position_fwd], [0], color='black', label='1: c.g. fwd')

# Plotting fuselage
ax.plot([-2+aircraft.x_cg_position_fwd, 1.5+aircraft.x_cg_position_fwd], [aircraft.w_out/2, aircraft.w_out/2], color='0.25', linewidth=0.8)
ax.plot([-2+aircraft.x_cg_position_fwd, 1.5+aircraft.x_cg_position_fwd], [-aircraft.w_out/2, -aircraft.w_out/2], color='0.25', linewidth=0.8)

ax.plot([0, -2.4+aircraft.x_cg_position_fwd], [0, aircraft.w_out/4], color='0.25', linewidth=0.8)
ax.plot([0, -2.4+aircraft.x_cg_position_fwd], [0, -aircraft.w_out/4], color='0.25', linewidth=0.8)

ax.plot([-2+aircraft.x_cg_position_fwd, -2.4+aircraft.x_cg_position_fwd], [aircraft.w_out/2, aircraft.w_out/4], color='0.25', linewidth=0.8)
ax.plot([-2+aircraft.x_cg_position_fwd, -2.4+aircraft.x_cg_position_fwd], [-aircraft.w_out/2, -aircraft.w_out/4], color='0.25', linewidth=0.8)

ax.set_xlim(-0.3, aircraft.l_f-1)
ax.set_ylim(-aircraft.b/2, aircraft.b/2+1)
ax.axhline(0, color='blue', linestyle='dotted')
# End plotting Fuselage

# Plotting z-cg. max location in top right corner
parameter_text = ax.text(0.95, 0.95, f'Max Z_position c.g.: {Z_position_cg:.2f} [m]',
                        transform=ax.transAxes, ha='right', va='top')

plt.text(aircraft.x_cg_position_aft+0.03, 0.2, 2, fontsize=8, va='center')
plt.text(aircraft.x_cg_position_fwd+0.03, 0.2, 1, fontsize=8, va='center')

# Plot Wing Position and planform
spanwise_pos = (aircraft.MAC_length - aircraft.rootchord) / -((aircraft.rootchord-aircraft.tipchord)/(aircraft.b/2))
ax.plot([aircraft.X_LEMAC, aircraft.X_LEMAC+aircraft.MAC_length], [spanwise_pos, spanwise_pos],color='0.25', linewidth=0.8)
ax.plot([aircraft.X_LEMAC - aircraft.x_lemac, aircraft.X_LEMAC - aircraft.x_lemac+aircraft.rootchord], [0, 0], color='0.25', linewidth=0.8)

x_root_quart = aircraft.X_LEMAC-aircraft.x_lemac + 1/4*aircraft.rootchord
ax.plot([x_root_quart-1/4*aircraft.tipchord, x_root_quart+3/4*aircraft.tipchord], [aircraft.b/2, aircraft.b/2], color='0.25', linewidth=0.8)
ax.plot([aircraft.X_LEMAC-aircraft.x_lemac, x_root_quart-1/4*aircraft.tipchord], [0, aircraft.b/2], color='0.25', linewidth=0.8)
ax.plot([aircraft.X_LEMAC-aircraft.x_lemac+aircraft.rootchord, x_root_quart+3/4*aircraft.tipchord], [0, aircraft.b/2], color='0.25', linewidth=0.8)

plt.grid()
plt.legend(loc="upper left")
plt.xlabel("x-position")
plt.ylabel("y_position")


# Show the plot
plt.show(block=True)
