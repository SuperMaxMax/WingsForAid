import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('.')

from parameters import UAV, atmosphere

aircraft = UAV('aircraft')
atm = atmosphere()

aircraft.x_cg_position_aft = aircraft.X_LEMAC + aircraft.X_cg_aft * aircraft.MAC_length
aircraft.x_cg_position_fwd = aircraft.X_LEMAC + aircraft.X_cg_fwd * aircraft.MAC_length

aircraft.position_landing_fwd = [2, 0]
aircraft.position_landing_back = [4, 1]

# Initial values for plotting
line_1 = None
line_2 = None
line_3 = None
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

    line_1 = ax.plot([x_front_max, x_front_max], [-0.2, 0.2], color='blue')[0]
    line_2 = ax.plot([x_front_min, x_front_min], [-0.2, 0.2], color='blue')[0]
    
    

def lateral_position_landing_gear(aircraft):
    global Z_position_cg
    global line_3

    if line_3 is not None:
        line_3.remove()

    d = aircraft.position_landing_back[0] - aircraft.position_landing_fwd[0]
    print(aircraft.position_landing_fwd[1])
    angle_alpha = np.arctan(aircraft.position_landing_back[1] / d)
    print(angle_alpha)
    c = (aircraft.x_cg_position_aft - aircraft.position_landing_fwd[0]) * np.sin(angle_alpha)
    angle_psi = 55 * (np.pi / 180)  # rad
    Z_position_cg = np.tan(angle_psi) * c
    parameter_text.set_text(f'Max Z_position_c.g.: {Z_position_cg:.2f} [m]')

    # pitch angle limit
    angle_theta_max = 15 * (np.pi/180)
    x_cg_aft_limit = aircraft.x_cg_position_aft + Z_position_cg * np.tan(angle_theta_max)
    line_3 = ax.plot([x_cg_aft_limit, x_cg_aft_limit], [-0.2, 0.2], color='blue')[0]
    plt.draw()


def height_landing_gear(aircraft):
    pass


def run(aircraft):
    longitudinal_position_landing_gear(aircraft)
    lateral_position_landing_gear(aircraft)
    height_landing_gear(aircraft)


def on_click(event):
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
            longitudinal_position_landing_gear(aircraft, event.xdata, event.ydata)
            lateral_position_landing_gear(aircraft)
            aircraft.position_landing_back = [event.xdata, event.ydata]
            # Call the lateral_position_landing_gear function to update Z_position_cg
            

        elif event.button == 3:  # Right-click
            if right_click_point is None:
                # Add a new right-click point
                right_click_point = ax.scatter([event.xdata], [0], color='orange', marker='^')
                
            else:
                # Update the position of the right-click point
                right_click_point.set_offsets([[event.xdata, 0]])

            # Call the lateral_position_landing_gear function to update Z_position_cg
            lateral_position_landing_gear(aircraft)
            aircraft.position_landing_fwd = [event.xdata, event.ydata]

        # Redraw the plot
        plt.draw()


# Enable interactive mode
plt.ion()

# Create a figure and an axis
fig, ax = plt.subplots()

# Create a scatter plot with a single point
point = ax.scatter([0], [0])

# Connect the 'button_press_event' to the 'on_click' function
cid = fig.canvas.mpl_connect('button_press_event', on_click)

# Set the plot limits
x_cg_point = ax.scatter([aircraft.x_cg_position_aft], [0], color='black', label='c.g. aft')
x_cg_point = ax.scatter([aircraft.x_cg_position_fwd], [0], color='black', label='c.g. fwd')
ax.plot([-1.5+aircraft.x_cg_position_fwd, 1.5+aircraft.x_cg_position_fwd], [aircraft.w_out/2, aircraft.w_out/2], color='grey', linewidth=0.5)
ax.plot([-1.5+aircraft.x_cg_position_fwd, 1.5+aircraft.x_cg_position_fwd], [-aircraft.w_out/2, -aircraft.w_out/2], color='grey', linewidth=0.5)
ax.set_xlim(0, aircraft.l_f)
ax.set_ylim(-aircraft.b/2 - 0.5, aircraft.b/2 + 0.5)
ax.axhline(0, color='blue', linestyle='dotted')

parameter_text = ax.text(0.95, 0.95, f'Max Z_position c.g.: {Z_position_cg:.2f} [m]',
                         transform=ax.transAxes, ha='right', va='top')

# Show the plot
plt.show(block=True)
