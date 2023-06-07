import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

def plot_wingbox():
    """
    Plot the wingbox of the aircraft
    """
    # Plot airfoil shape
    plt.plot(x, y_u, 'b', zorder=0)
    plt.plot(x, y_l, 'b', zorder=0)
    # Plot front spar, aft spar, top skin and bottom skin (in that order)
    plt.gca().add_patch(patches.Rectangle((xp[3]-t_f_spar/2, yp[3]), t_f_spar, l_f_spar, angle=0.0, linewidth=0, edgecolor='none', facecolor='r'))
    plt.gca().add_patch(patches.Rectangle((xp[2]-t_a_spar/2, yp[2]), t_a_spar, l_a_spar, angle=0.0, linewidth=0, edgecolor='none', facecolor='r'))
    plt.gca().add_patch(patches.Rectangle((xp[0], yp[0]-t_top/2), l_top, t_top, angle=angle_top*180/np.pi, linewidth=0, edgecolor='none', facecolor='r'))
    plt.gca().add_patch(patches.Rectangle((xp[3], yp[3]-t_bot/2), l_bot, t_bot, angle=angle_bot*180/np.pi, linewidth=0, edgecolor='none', facecolor='r'))
    # Plot corner points
    plt.scatter(xp, yp, zorder=2)
    # Plot centroids
    plt.scatter(xc, yc, c='k')
    plt.scatter(xc_stringer_top, yc_stringer_top, c='k', zorder=10)
    plt.scatter(xc_stringer_bot, yc_stringer_bot, c='k', zorder=10)
    # Plot stringers
    for i in range(len(xco_stringer_top)):
        if i == 0:
            plot_stringer(xco_stringer_top[i], yco_stringer_top[i], 180+angle_top*180/np.pi, mirror=True)
        else:
            plot_stringer(xco_stringer_top[i], yco_stringer_top[i], 180+angle_top*180/np.pi)
    for i in range(len(xco_stringer_bot)):
        if i == len(xc_stringer_bot)-1:
            plot_stringer(xco_stringer_bot[i], yco_stringer_bot[i], angle_bot*180/np.pi, mirror=True)
        else:
            plot_stringer(xco_stringer_bot[i], yco_stringer_bot[i], angle_bot*180/np.pi)
    plt.axis('equal')
    plt.grid()
    plt.show()

def plot_stringer_separate():
    """
    Plots the L-stringer separately in a figure
    """
    plt.figure()
    face_1 = patches.Rectangle((0, 0), a, t2, angle=0.0, linewidth=0, edgecolor='none', facecolor='g')
    plt.gca().add_patch(face_1)
    face_2 = patches.Rectangle((0, 0), t1, b, angle=0.0, linewidth=0, edgecolor='none', facecolor='g')
    plt.gca().add_patch(face_2)
    plt.scatter(xc_s, yc_s, c='k', zorder=10, label=f'centroid: ({xc_s:.4f}, {yc_s:.4f}) m')
    plt.title('Cross section of L-stringer')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.grid()
    plt.legend()
    plt.show()

def rotate_centroid(angle):
    """
    Rotate centroid of L-stringer around bottom left corner
    angle = angle of rotation
    """
    rot_mat = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    c = np.array([xc_s, yc_s])
    c = np.matmul(rot_mat, c)

    return c[0], c[1]

def plot_stringer(x, y, ang, mirror=False):
    """
    Plot stringer with rectangles in wingbox plot at given x and y for centroid
    x = x coordinate of centroid
    y = y coordinate of centroid
    ang = angle of rotation of L-stringer (around bottom left corner)
    mirror = True if stringer is needs to be mirrored
    """

    if mirror == False:
        face_1 = patches.Rectangle((x, y), a, t2, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_1)
        face_2 = patches.Rectangle((x, y), t1, b, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_2)
    else:
        face_1 = patches.Rectangle((x, y), -a, t2, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_1)
        face_2 = patches.Rectangle((x, y), -t1, b, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_2)

def chord(y):
    """
    Calculate chord length at given spanwise location
    y = spanwise location
    """

    return ac.rootchord-((ac.rootchord-ac.tipchord)/(ac.b/2))*y

def airfoil(f_spar, a_spar, span):
    """
    Returns the upper and lower surface of an airfoil with a common x-axis
    f_spar = location of front spar (as fraction of chord)
    a_spar = location of aft spar (as fraction of chord)
    span = spanwise location of airfoil
    """
    c = np.linspace(0, 1, 1000)
    m = float(airfoil_num[0]) / 100
    p = float(airfoil_num[1]) / 10
    c1 = c[c <= p]
    c2 = c[c > p]
    
    y_c = m/p**2*(2*p*c1 - c1**2)
    theta = np.arctan((2*m/p**2)*(p-c1))
    y_t = 5*float('0.'+airfoil_num[2:])*(0.2969*np.sqrt(c1)-0.126*c1-0.3516*c1**2+0.2843*c1**3-0.1015*c1**4)

    y_u = y_c + y_t*np.cos(theta)
    x_u = c1 - y_t*np.sin(theta)
    y_l = y_c - y_t*np.cos(theta)
    x_l = c1 + y_t*np.sin(theta)

    y_c = m/((1-p)**2)*((1-2*p)+2*p*c2-c2**2)
    theta = np.arctan((2*m/(1-p)**2)*(p-c2))
    y_t = 5*float('0.'+airfoil_num[2:])*(0.2969*np.sqrt(c2)-0.126*c2-0.3516*c2**2+0.2843*c2**3-0.1015*c2**4)

    y_u = np.concatenate((y_u, y_c + y_t*np.cos(theta)))
    x_u = np.concatenate((x_u, c2 - y_t*np.sin(theta)))
    y_l = np.concatenate((y_l, y_c - y_t*np.cos(theta)))
    x_l = np.concatenate((x_l, c2 + y_t*np.sin(theta)))

    x_gen = np.linspace(0, 1, 1000)
    y_u = np.interp(x_gen, x_u, y_u)
    y_l = np.interp(x_gen, x_l, y_l)

    # get the 4 cornerpoints of the wingbox
    x = np.array([f_spar, a_spar, a_spar, f_spar])
    y = np.array([np.interp(f_spar, x_gen, y_u), np.interp(a_spar, x_gen, y_u), np.interp(a_spar, x_gen, y_l), np.interp(f_spar, x_gen, y_l)])

    x_gen *= chord(span)
    y_u *= chord(span)
    y_l *= chord(span)
    x *= chord(span)
    y *= chord(span)
    
    return x_gen, y_u, y_l, x, y

def L_stringer(a, b, t1, t2, angle):
    """
    Calculates the area, centroid and second moment of area of an L-stringer
    a = length of bottom flange
    b = length of side flange
    t1 = thickness of side flange
    t2 = thickness of bottom flange
    angle = angle of rotation of L-stringer (around bottom left corner)
    """

    # Area
    A_s = a*t2+(b-t2)*t1

    # Centroid -> wrt to bottom left corner, this is also the rotation point
    # centroid no angle
    xc_s = (a*t2*a/2+(b-t2)*t1*t1/2)/A_s
    yc_s = (a*t2*t2/2+(b-t2)*t1*((b-t2)/2+t2))/A_s
    c = np.array([xc_s, yc_s])
    # centroid of individual parts
    c1 = np.array([a/2, t2/2])
    c2 = np.array([t1/2, (b-t2)/2+t2])
    # define rotation matrix
    rot_mat = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    # rotate individual parts
    c = np.matmul(rot_mat, c)
    c1 = np.matmul(rot_mat, c1)
    c2 = np.matmul(rot_mat, c2)
    # save rotated centroid
    xc_s_a = c[0]
    yc_s_a = c[1]
    # save distance from centroid to individual parts for parallel axis term
    dx_c_c1 = c[0]-c1[0]
    dy_c_c1 = c[1]-c1[1]
    dx_c_c2 = c[0]-c2[0]
    dy_c_c2 = c[1]-c2[1]

    # Second moment of area about centroid
    Ix_s_a = a*t2*(t2**2*np.cos(angle)+a**2*np.sin(angle))/12+a*t2*(dy_c_c1)**2+(b-t2)*t1*(t1**2*np.cos(angle)+(b-t2)**2*np.sin(angle))/12+(b-t2)*t1*(dy_c_c2)**2
    Iy_s_a = a*t2*(a**2*np.cos(angle)+t2**2*np.sin(angle))/12+a*t2*(dx_c_c1)**2+(b-t2)*t1*((b-t2)**2*np.cos(angle)+t1**2*np.sin(angle))/12+(b-t2)*t1*(dx_c_c2)**2

    return A_s, xc_s_a, yc_s_a, Ix_s_a, Iy_s_a

ac = UAV("aircraft")

### INPUT ###
airfoil_num = '4415'

f_spar = 0.25
t_f_spar = 0.005
a_spar = 0.75
t_a_spar = 0.005

t_top = 0.005
t_bot = 0.005

# chose L stringer for now
n_stringers = 6          # this is additional to already having 4 stringers at the corners
a = 0.015
b = 0.04
t1 = 0.005
t2 = 0.005

n_ribs = 10
t_ribs = 0.01
spacing_ribs = 0.1

### CALCULATIONS ###
# xp and yp are the 4 corner points -> [top left, top right, bot right, bot left]
x, y_u, y_l, xp, yp = airfoil(f_spar, a_spar, 0)

# Calculate properties of L-stringer (no angle)
A_s, xc_s, yc_s, Ix_s, Iy_s = L_stringer(a, b, t1, t2, 0)

# Calculate centroid, zero point coincides with origin plot
# Spar contributions
l_f_spar = yp[0]-yp[3]
l_a_spar = yp[1]-yp[2]

xc_f_spar = xp[0]
yc_f_spar = yp[3]+l_f_spar/2
xc_a_spar = xp[1]
yc_a_spar = yp[2]+l_a_spar/2

# Skin top and bot contributions
l_top = ((xp[1]-xp[0])**2+(yp[1]-yp[0])**2)**0.5
l_bot = ((xp[2]-xp[3])**2+(yp[2]-yp[3])**2)**0.5

angle_top = np.arctan((yp[1]-yp[0])/(xp[1]-xp[0]))
angle_bot = np.arctan((yp[2]-yp[3])/(xp[2]-xp[3]))

xc_top = xp[0]+(xp[1]-xp[0])/2
yc_top = yp[0]+(yp[1]-yp[0])/2
xc_bot = xp[3]+(xp[2]-xp[3])/2
yc_bot = yp[3]+(yp[2]-yp[3])/2

# Stringer contributions
n_stringers_top = int(n_stringers/2)
n_stringers_bot = n_stringers - n_stringers_top

dy_stringer_left_top = -(t_top/2)/np.cos(angle_top)+np.tan(angle_top)*t_f_spar/2

x_stringer_spacing_top = (xp[1]-xp[0]-t_f_spar/2-t_a_spar/2)/(n_stringers_top+1)
y_stringer_spacing_top = x_stringer_spacing_top*np.tan(angle_top)
xco_stringer_top = [i*x_stringer_spacing_top+xp[0]+t_f_spar/2 for i in range(n_stringers_top+2)]
yco_stringer_top = [i*y_stringer_spacing_top+yp[0]+dy_stringer_left_top for i in range(n_stringers_top+2)]

xc_stringer_top = [xp[0]+t_f_spar/2+rotate_centroid(np.pi+angle_top)[0]+2*xc_s*np.cos(angle_top)] + [i*x_stringer_spacing_top+xp[0]+t_f_spar/2+rotate_centroid(np.pi+angle_top)[0] for i in range(1, n_stringers_top+2)]
yc_stringer_top = [yp[0]+dy_stringer_left_top+rotate_centroid(np.pi+angle_top)[1]-2*yc_s*np.sin(angle_top)] + [i*y_stringer_spacing_top+yp[0]+dy_stringer_left_top+rotate_centroid(np.pi+angle_top)[1] for i in range(1, n_stringers_top+2)]

dy_stringer_left_bot = (t_bot/2)/np.cos(angle_bot)+np.tan(angle_bot)*t_f_spar/2

x_stringer_spacing_bot = (xp[2]-xp[3]-t_f_spar/2-t_a_spar/2)/(n_stringers_bot+1)
y_stringer_spacing_bot = x_stringer_spacing_bot*np.tan(angle_bot)
xco_stringer_bot = [i*x_stringer_spacing_bot+xp[3]+t_f_spar/2 for i in range(n_stringers_bot+2)]
yco_stringer_bot = [i*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot for i in range(n_stringers_bot+2)]

xc_stringer_bot = [i*x_stringer_spacing_bot+xp[3]+t_f_spar/2+rotate_centroid(angle_bot)[0] for i in range(n_stringers_bot+1)] + [(n_stringers_bot+1)*x_stringer_spacing_bot+xp[3]+t_f_spar/2+rotate_centroid(angle_bot)[0]-2*xc_s*np.cos(angle_bot)]
yc_stringer_bot = [i*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot+rotate_centroid(angle_bot)[1] for i in range(n_stringers_bot+1)] + [(n_stringers_bot+1)*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot+rotate_centroid(angle_bot)[1]-2*xc_s*np.sin(angle_bot)]

# Final centroid
xc = (xc_f_spar*l_f_spar*t_f_spar+xc_a_spar*l_a_spar*t_a_spar+xc_top*l_top*t_top+xc_bot*l_bot*t_bot+sum([xc_stringer_top[i]*A_s for i in range(n_stringers_top+2)])+sum([xc_stringer_bot[i]*A_s for i in range(n_stringers_bot+2)]))/(l_f_spar*t_f_spar+l_a_spar*t_a_spar+l_top*t_top+l_bot*t_bot+sum([A_s for i in range(n_stringers_top+2)])+sum([A_s for i in range(n_stringers_bot+2)]))
yc = (yc_f_spar*l_f_spar*t_f_spar+yc_a_spar*l_a_spar*t_a_spar+yc_top*l_top*t_top+yc_bot*l_bot*t_bot+sum([yc_stringer_top[i]*A_s for i in range(n_stringers_top+2)])+sum([yc_stringer_bot[i]*A_s for i in range(n_stringers_bot+2)]))/(l_f_spar*t_f_spar+l_a_spar*t_a_spar+l_top*t_top+l_bot*t_bot+sum([A_s for i in range(n_stringers_top+2)])+sum([A_s for i in range(n_stringers_bot+2)]))

plot_wingbox()
plot_stringer_separate()

# list of things to do
# - fix relfection issue
# - calculate second moment of area as function of span
# - incorporate ribs
# - incorporate variable amount of stringers per rib section?
# - calculate torsional stiffness as function of span
# - define bounds for all parameters and optimize for minimum weight, given a deflection and twist constraint


