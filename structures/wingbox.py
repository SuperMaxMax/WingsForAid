import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
from material_dictionary import rank_material, material_df, material
from loading_diagram_wing import loading_tension, loading_compression
import numpy as np
from scipy import interpolate
from scipy.integrate import cumulative_trapezoid, trapezoid
import matplotlib.pyplot as plt
from matplotlib import patches

def plot_wingbox(spanwise_pos):
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
    # plt.scatter(x_mid_top, y_mid_top, zorder=20, c='k')
    # plt.scatter(x_mid_bot, y_mid_bot, zorder=20, c='k')
    # Plot centroids
    plt.scatter(xc, yc, c='k')
    plt.scatter(xc_stringer_top, yc_stringer_top, c='k', zorder=10)
    plt.scatter(xc_stringer_bot, yc_stringer_bot, c='k', zorder=10)
    # Plot stringers
    for i in range(len(xco_stringer_top)):
        if i == 0:
            plot_L_stringer(xco_stringer_top[i], yco_stringer_top[i], 180+angle_top*180/np.pi, mirror=True)
        else:
            plot_L_stringer(xco_stringer_top[i], yco_stringer_top[i], 180+angle_top*180/np.pi)
    for i in range(len(xco_stringer_bot)):
        if i == len(xc_stringer_bot)-1:
            plot_L_stringer(xco_stringer_bot[i], yco_stringer_bot[i], angle_bot*180/np.pi, mirror=True)
        else:
            plot_L_stringer(xco_stringer_bot[i], yco_stringer_bot[i], angle_bot*180/np.pi)
    plt.axis('equal')
    plt.title('Wingbox at spanwise location: {:.2f} m'.format(spanwise_pos))
    plt.grid()
    plt.show()

def plot_stringer_separate():
    """
    Plots the L-stringer separately in a figure
    """
    plt.figure()
    face_1 = patches.Rectangle((0, 0), a, t_stringer, angle=0.0, linewidth=0, edgecolor='none', facecolor='g')
    plt.gca().add_patch(face_1)
    face_2 = patches.Rectangle((0, 0), t_stringer, b, angle=0.0, linewidth=0, edgecolor='none', facecolor='g')
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

def plot_L_stringer(x, y, ang, mirror=False):
    """
    Plot stringer with rectangles in wingbox plot at given x and y for centroid
    x = x coordinate of centroid
    y = y coordinate of centroid
    ang = angle of rotation of L-stringer (around bottom left corner)
    mirror = True if stringer is needs to be mirrored
    """

    if mirror == False:
        face_1 = patches.Rectangle((x, y), a, t_stringer, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_1)
        face_2 = patches.Rectangle((x, y), t_stringer, b, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_2)
    else:
        face_1 = patches.Rectangle((x, y), -a, t_stringer, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
        plt.gca().add_patch(face_1)
        face_2 = patches.Rectangle((x, y), -t_stringer, b, angle=ang, linewidth=0, edgecolor='none', facecolor='g')
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

def L_stringer(a, b, t_stringer, angle):
    """
    Calculates the area, centroid and second moment of area of an L-stringer
    a = length of bottom flange
    b = length of side flange
    t_stringer = thickness of side flange
    t_stringer = thickness of bottom flange
    angle = angle of rotation of L-stringer (around bottom left corner)
    """

    # Area
    A_s = a*t_stringer+(b-t_stringer)*t_stringer

    # Centroid -> wrt to bottom left corner, this is also the rotation point
    # centroid no angle
    xc_s = (a*t_stringer*a/2+(b-t_stringer)*t_stringer*t_stringer/2)/A_s
    yc_s = (a*t_stringer*t_stringer/2+(b-t_stringer)*t_stringer*((b-t_stringer)/2+t_stringer))/A_s
    c = np.array([xc_s, yc_s])
    # centroid of individual parts
    c1 = np.array([a/2, t_stringer/2])
    c2 = np.array([t_stringer/2, (b-t_stringer)/2+t_stringer])
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
    Ix_s_a = a*t_stringer*(t_stringer**2*np.cos(angle)+a**2*np.sin(angle))/12+a*t_stringer*(dy_c_c1)**2+(b-t_stringer)*t_stringer*(t_stringer**2*np.cos(angle)+(b-t_stringer)**2*np.sin(angle))/12+(b-t_stringer)*t_stringer*(dy_c_c2)**2
    Iy_s_a = a*t_stringer*(a**2*np.cos(angle)+t_stringer**2*np.sin(angle))/12+a*t_stringer*(dx_c_c1)**2+(b-t_stringer)*t_stringer*((b-t_stringer)**2*np.cos(angle)+t_stringer**2*np.sin(angle))/12+(b-t_stringer)*t_stringer*(dx_c_c2)**2

    return A_s, xc_s_a, yc_s_a, Ix_s_a, Iy_s_a

def ks_values(filename):
    ks_x = np.array([])
    ks_y = np.array([])
    file = open(f"ks_values/{filename}.txt", "r")
    for line in file:
        line = line.rstrip().split(",")
        ks_x = np.append(ks_x, float(line[0]))
        ks_y = np.append(ks_y, float(line[1]))

    ks = interpolate.interp1d(ks_x, ks_y, kind='cubic',fill_value="extrapolate")

    return ks, ks_y[-1]

def increase_values(values, increase):
    for i, value in enumerate(values):
        values[i] += increase
    
    return values

ac = UAV("aircraft")
plot = True

#### INPUT ####
airfoil_num = '4415'        # [-] NACA airfoil number (4 digit)

f_spar = 0.25               # [-] location of front spar (as fraction of chord) 
t_f_spar = 0.0001           # [m] thickness of front spar
a_spar = 0.75               # [-] location of aft spar (as fraction of chord)
t_a_spar = 0.0001           # [m] thickness of aft spar

t_top = 0.0001              # [m] thickness of top panel
t_bot = 0.0001              # [m] thickness of bottom panel

n_stringers = 0             # [-] this is additional to already having 4 stringers at the corners
a = 0.00015                 # [m] length of bottom flange
b = 0.0003                  # [m] length of side flange
t_stringer = 0.00001        # [m] thickness of stringer (both flanges)

n_ribs = 1                  # [-] number of ribs
t_ribs = 0.001              # [m] thickness of ribs

n_stringer_decr_ribs = 2    # [-] decrease in number of stringers per x ribs
stringer_decr_x_ribs = 1    # [-] decrease stringers by n_stringer_decr_ribs every x ribs

#### MATERIAL PROPERTIES ####
# Material choice for strut weight, add column for material index: sqrt(E)/rho -> the higher the value the better -> False
material_df['sqrt(E)/rho'] = (material_df['E'])**(1/3)/material_df['density']
# density, raw cost, eco cost, co2, yield stress, E, Kc, sqrt(E)/rho
prop_weights_strut = [0, 0.4, 0.1, 0.05, 0, 0, 0.05, 0.4, 0]
mat = rank_material(prop_weights_strut, [False, False])[0]
print(mat)

#### CALCULATIONS ####
step_size = 0.0001          # [m] step size for integration
### Spanwise position ###
spanwise_pos = np.linspace(0, ac.b/2, 1000)
# split spanwise position in multiple parts based on ribs
spanwise_pos_sec = np.array_split(spanwise_pos, n_ribs+1)

### Values to save ###
Ixs = np.zeros(len(spanwise_pos))
Iys = np.zeros(len(spanwise_pos))
Js = np.zeros(len(spanwise_pos))
xcs = np.zeros(len(spanwise_pos))
ycs = np.zeros(len(spanwise_pos))
As = np.zeros(len(spanwise_pos))
MOFs = np.ones(18)*1.1

### Ks values buckling ###
ks_A_clamped, ks_A_clamped_min = ks_values('ks_A_clamped')
ks_A_ss, ks_A_ss_min = ks_values('ks_A_ss')
ks_clamped, ks_clamped_min = ks_values('ks_clamped')

failure_tot = True
run_once_more = 0

while failure_tot == True:
    failure_ele = False
    n_stringers = 6
    w_ribs = 0
    for rib_number, j in enumerate(spanwise_pos_sec):
        for i in j:
            ### Airfoil ###
            # get spar sizes at tip -> for stringer size limitation
            x_tip, y_u_tip, y_l_tip, xp_tip, yp_tip = airfoil(f_spar, a_spar, ac.b/2)

            # xp and yp are the 4 corner points -> [top left, top right, bot right, bot left]
            x, y_u, y_l, xp, yp = airfoil(f_spar, a_spar, i)

            ### Stringer ###
            # calculate properties of L-stringer (no angle)
            A_s, xc_s, yc_s, Ix_s, Iy_s = L_stringer(a, b, t_stringer, 0)

            ### Centroid ###
            ## Spar contributions
            # spar lenghts
            l_f_spar = yp[0]-yp[3]
            l_a_spar = yp[1]-yp[2]

            # spar centroids
            xc_f_spar = xp[0]
            yc_f_spar = yp[3]+l_f_spar/2
            xc_a_spar = xp[1]
            yc_a_spar = yp[2]+l_a_spar/2

            ## Panel contributions
            # top and bottom panel length
            l_top = ((xp[1]-xp[0])**2+(yp[1]-yp[0])**2)**0.5
            l_bot = ((xp[2]-xp[3])**2+(yp[2]-yp[3])**2)**0.5

            # top and bottom panel angles
            angle_top = np.arctan((yp[1]-yp[0])/(xp[1]-xp[0]))
            angle_bot = np.arctan((yp[2]-yp[3])/(xp[2]-xp[3]))

            # top and bottom panel centroids
            xc_top = xp[0]+(xp[1]-xp[0])/2
            yc_top = yp[0]+(yp[1]-yp[0])/2
            xc_bot = xp[3]+(xp[2]-xp[3])/2
            yc_bot = yp[3]+(yp[2]-yp[3])/2

            ## Stringer contributions
            # number of stringers on top and bottom
            n_stringers_top = int(2*n_stringers/3)
            n_stringers_bot = n_stringers - n_stringers_top

            # y-distance between cornerpoint (xp[0], yp[0]) and intersection of inner edge of the front spar and top panel, this point coincides with the cornerpoint of the first stringer
            dy_stringer_left_top = -(t_top/2)/np.cos(angle_top)+np.tan(angle_top)*t_f_spar/2

            # get cornerpoints (xco, yco) of all stringers on top, stringer will be placed with bottom left corner at this point, and then rotated such that it aligns with the top panel
            x_stringer_spacing_top = (xp[1]-xp[0]-t_f_spar/2-t_a_spar/2)/(n_stringers_top+1)
            y_stringer_spacing_top = x_stringer_spacing_top*np.tan(angle_top)
            xco_stringer_top = [i*x_stringer_spacing_top+xp[0]+t_f_spar/2 for i in range(n_stringers_top+2)]
            yco_stringer_top = [i*y_stringer_spacing_top+yp[0]+dy_stringer_left_top for i in range(n_stringers_top+2)]

            # convert the cornerpoints to centroids, this is done by rotating the centroid of a straight L-stringer around the bottom left corner, the first stringer is a special case because of the corner
            xc_stringer_top = [xp[0]+t_f_spar/2+rotate_centroid(np.pi+angle_top)[0]+2*xc_s*np.cos(angle_top)] + [i*x_stringer_spacing_top+xp[0]+t_f_spar/2+rotate_centroid(np.pi+angle_top)[0] for i in range(1, n_stringers_top+2)]
            yc_stringer_top = [yp[0]+dy_stringer_left_top+rotate_centroid(np.pi+angle_top)[1]-2*yc_s*np.sin(angle_top)] + [i*y_stringer_spacing_top+yp[0]+dy_stringer_left_top+rotate_centroid(np.pi+angle_top)[1] for i in range(1, n_stringers_top+2)]

            # y-distance between cornerpoint (xp[3], yp[3]) and intersection of inner edge of the front spar and bottom panel, this point coincides with the cornerpoint of the first stringer
            dy_stringer_left_bot = (t_bot/2)/np.cos(angle_bot)+np.tan(angle_bot)*t_f_spar/2

            # get cornerpoints (xco, yco) of all stringers on bottom, stringer will be placed with bottom left corner at this point, and then rotated such that it aligns with the bottom panel
            x_stringer_spacing_bot = (xp[2]-xp[3]-t_f_spar/2-t_a_spar/2)/(n_stringers_bot+1)
            y_stringer_spacing_bot = x_stringer_spacing_bot*np.tan(angle_bot)
            xco_stringer_bot = [i*x_stringer_spacing_bot+xp[3]+t_f_spar/2 for i in range(n_stringers_bot+2)]
            yco_stringer_bot = [i*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot for i in range(n_stringers_bot+2)]

            # convert the cornerpoints to centroids, this is done by rotating the centroid of a straight L-stringer around the bottom left corner, the last stringer is a special case because of the corner
            xc_stringer_bot = [i*x_stringer_spacing_bot+xp[3]+t_f_spar/2+rotate_centroid(angle_bot)[0] for i in range(n_stringers_bot+1)] + [(n_stringers_bot+1)*x_stringer_spacing_bot+xp[3]+t_f_spar/2+rotate_centroid(angle_bot)[0]-2*xc_s*np.cos(angle_bot)]
            yc_stringer_bot = [i*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot+rotate_centroid(angle_bot)[1] for i in range(n_stringers_bot+1)] + [(n_stringers_bot+1)*y_stringer_spacing_bot+yp[3]+dy_stringer_left_bot+rotate_centroid(angle_bot)[1]-2*xc_s*np.sin(angle_bot)]

            ## Final centroid | (xc*A)/(A_tot) = (x1*A1+x2*A2+...+xn*An)/(A1+A2+...+An)
            A_tot = (l_f_spar*t_f_spar+l_a_spar*t_a_spar+l_top*t_top+l_bot*t_bot+sum([A_s for i in range(n_stringers_top+2)])+sum([A_s for i in range(n_stringers_bot+2)]))
            xc = (xc_f_spar*l_f_spar*t_f_spar+xc_a_spar*l_a_spar*t_a_spar+xc_top*l_top*t_top+xc_bot*l_bot*t_bot+sum([xc_stringer_top[i]*A_s for i in range(n_stringers_top+2)])+sum([xc_stringer_bot[i]*A_s for i in range(n_stringers_bot+2)]))/A_tot
            yc = (yc_f_spar*l_f_spar*t_f_spar+yc_a_spar*l_a_spar*t_a_spar+yc_top*l_top*t_top+yc_bot*l_bot*t_bot+sum([yc_stringer_top[i]*A_s for i in range(n_stringers_top+2)])+sum([yc_stringer_bot[i]*A_s for i in range(n_stringers_bot+2)]))/A_tot

            xcs[np.where(spanwise_pos == i)] = xc
            ycs[np.where(spanwise_pos == i)] = yc
            As[np.where(spanwise_pos == i)] = A_tot

            ### Second moment of area ###
            ## Spar contributions
            # spar second moments of area
            Ix_f_spar = t_f_spar*l_f_spar**3/12
            Iy_f_spar = l_f_spar*t_f_spar**3/12
            Ix_a_spar = t_a_spar*l_a_spar**3/12
            Iy_a_spar = l_a_spar*t_a_spar**3/12

            # spar parallel axis terms
            Ix_f_spar += l_f_spar*t_f_spar*(yc_f_spar-yc)**2
            Iy_f_spar += t_f_spar*l_f_spar*(xc_f_spar-xc)**2
            Ix_a_spar += l_a_spar*t_a_spar*(yc_a_spar-yc)**2
            Iy_a_spar += t_a_spar*l_a_spar*(xc_a_spar-xc)**2

            ## Panel contributions
            # top and bottom panel second moments of area
            Ix_top = t_top*l_top*(l_top**2*np.cos(angle_top)+t_top**2*np.cos(angle_top))/12
            Iy_top = t_top*l_top*(l_top**2*np.sin(angle_top)+t_top**2*np.sin(angle_top))/12
            Ix_bot = t_bot*l_bot*(l_bot**2*np.cos(angle_bot)+t_bot**2*np.cos(angle_bot))/12
            Iy_bot = t_bot*l_bot*(l_bot**2*np.sin(angle_bot)+t_bot**2*np.sin(angle_bot))/12

            # top and bottom panel parallel axis terms
            Ix_top += t_top*l_top*(yc_top-yc)**2
            Iy_top += t_top*l_top*(xc_top-xc)**2
            Ix_bot += t_bot*l_bot*(yc_bot-yc)**2
            Iy_bot += t_bot*l_bot*(xc_bot-xc)**2

            ## Stringer contributions
            # stringer second moments of area
            Ix_s_t, Iy_s_t  = L_stringer(a, b, t_stringer, angle_top)[3:5]
            Ix_s_b, Iy_s_b  = L_stringer(a, b, t_stringer, angle_bot)[3:5]
            Ix_s = (n_stringers_top+2)*Ix_s_t+(n_stringers_bot+2)*Ix_s_b
            Iy_s = (n_stringers_top+2)*Iy_s_t+(n_stringers_bot+2)*Iy_s_b

            # stringer parallel axis terms
            Ix_s += sum([A_s*(yc_stringer_top[i]-yc)**2 for i in range(n_stringers_top+2)])+sum([A_s*(yc_stringer_bot[i]-yc)**2 for i in range(n_stringers_bot+2)])
            Iy_s += sum([A_s*(xc_stringer_top[i]-xc)**2 for i in range(n_stringers_top+2)])+sum([A_s*(xc_stringer_bot[i]-xc)**2 for i in range(n_stringers_bot+2)])

            ## Final second moment of area
            Ix = Ix_f_spar+Ix_a_spar+Ix_top+Ix_bot+Ix_s
            Iy = Iy_f_spar+Iy_a_spar+Iy_top+Iy_bot+Iy_s

            Ixs[np.where(spanwise_pos == i)] = Ix
            Iys[np.where(spanwise_pos == i)] = Iy

            ### Torsional constant ###
            # horizontal length wingbox
            l_wb = (xp[1]-xp[0])
            
            # enclosed area | stringers not taken into account -> conservative
            A_enc = l_a_spar*l_wb+abs(yp[2]-yp[3])*l_wb/2+abs(yp[0]-yp[1])*l_wb/2
            J = (4*A_enc**2)/(l_f_spar/t_f_spar+l_top/t_top+l_a_spar/t_a_spar+l_bot/t_bot)

            Js[np.where(spanwise_pos == i)] = J

            ### Failure analysis ###
            ## Compressive and tensile strength failure ##
            # stress due to normal force
            sigma_n = loading_tension[0][np.where(spanwise_pos == i)[0][0]]/(l_top*t_top+l_bot*t_bot+l_f_spar*t_f_spar+l_a_spar*t_a_spar+(n_stringers+4)*A_s)

            # stress due to bending (always taking worst case scenario)
            # front spar | calculate at both ends of the spar
            sigma_b_f_spar_top = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_f_spar+l_f_spar/2-yc)/Ix
            sigma_b_f_spar_bot = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_f_spar-l_f_spar/2-yc)/Ix
            # calculate minimum and maximum stress from bending and add normal stress
            sigma_f_spar_c = min(sigma_b_f_spar_top, sigma_b_f_spar_bot)+sigma_n
            sigma_f_spar_t = max(sigma_b_f_spar_top, sigma_b_f_spar_bot)+sigma_n

            # aft spar | calculate at both ends of the spar
            sigma_b_a_spar_top = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_a_spar+l_a_spar/2-yc)/Ix
            sigma_b_a_spar_bot = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_a_spar-l_a_spar/2-yc)/Ix
            # calculate minimum and maximum stress from bending and add normal stress
            sigma_a_spar_c = min(sigma_b_a_spar_top, sigma_b_a_spar_bot)+sigma_n
            sigma_a_spar_t = max(sigma_b_a_spar_top, sigma_b_a_spar_bot)+sigma_n

            # top panel | calculate at both ends of the panel
            sigma_b_top_panel_1 = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_top+(yp[0]-yp[1])/2-yc)/Ix
            sigma_b_top_panel_2 = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_top-(yp[0]-yp[1])/2-yc)/Ix
            # calculate minimum and maximum stress from bending and add normal stress
            sigma_top_panel_c = min(sigma_b_top_panel_1, sigma_b_top_panel_2)+sigma_n
            sigma_top_panel_t = max(sigma_b_top_panel_1, sigma_b_top_panel_2)+sigma_n

            # bottom panel | calculate at both ends of the panel
            sigma_b_bot_panel_1 = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_bot+(yp[2]-yp[3])/2-yc)/Ix
            sigma_b_bot_panel_2 = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(yc_bot-(yp[2]-yp[3])/2-yc)/Ix
            # calculate minimum and maximum stress from bending and add normal stress
            sigma_bot_panel_c = min(sigma_b_bot_panel_1, sigma_b_bot_panel_2)+sigma_n
            sigma_bot_panel_t = max(sigma_b_bot_panel_1, sigma_b_bot_panel_2)+sigma_n

            if sigma_f_spar_c < 0:
                MOF_sigma_f_spar_c = material[mat]['yield stress']*10**6/abs(sigma_f_spar_c)
                MOFs[0] = MOF_sigma_f_spar_c
                if MOF_sigma_f_spar_c < 1:
                    print(f"Front spar compressive strength failure, MOF: {MOF_sigma_f_spar_c:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
            
            if sigma_f_spar_t > 0:
                MOF_sigma_f_spar_t = material[mat]['yield stress']*10**6/abs(sigma_f_spar_t)
                MOFs[1] = MOF_sigma_f_spar_t
                MOF_fracture_f_spar = (material[mat]['Kc']*10**6/np.sqrt(np.pi*0.005))/abs(sigma_f_spar_t)
                MOFs[2] = MOF_fracture_f_spar
                if MOF_sigma_f_spar_t < 1:
                    print(f"Front spar tensile strength failure, MOF: {MOF_sigma_f_spar_t:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break

                if MOF_fracture_f_spar < 1:
                    print(f"Front spar fracture failure, MOF: {MOF_fracture_f_spar:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
            
            if sigma_a_spar_c < 0:
                MOF_sigma_a_spar = material[mat]['yield stress']*10**6/abs(sigma_a_spar_c)
                MOFs[3] = MOF_sigma_a_spar
                if MOF_sigma_a_spar < 1:
                    print(f"Aft spar compressive strength failure, MOF: {MOF_sigma_a_spar:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break

            if sigma_a_spar_t > 0:
                MOF_sigma_a_spar_t = material[mat]['yield stress']*10**6/abs(sigma_a_spar_t)
                MOFs[4] = MOF_sigma_a_spar_t
                MOF_fracture_a_spar = (material[mat]['Kc']*10**6/np.sqrt(np.pi*0.005))/abs(sigma_a_spar_t)
                MOFs[5] = MOF_fracture_a_spar
                if MOF_sigma_a_spar_t < 1:
                    print(f"Aft spar tensile strength failure, MOF: {MOF_sigma_a_spar_t:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
                if MOF_fracture_a_spar < 1:
                    print(f"Aft spar fracture failure, MOF: {MOF_fracture_a_spar:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
            
            if sigma_top_panel_c < 0:
                MOF_sigma_top_panel = material[mat]['yield stress']*10**6/abs(sigma_top_panel_c)
                MOFs[6] = MOF_sigma_top_panel
                if MOF_sigma_top_panel < 1:
                    print(f"Top panel compressive strength failure, MOF: {MOF_sigma_top_panel:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
            
            if sigma_bot_panel_t > 0:
                MOF_sigma_bot_panel_t = material[mat]['yield stress']*10**6/abs(sigma_bot_panel_t)
                MOFs[7] = MOF_sigma_bot_panel_t
                MOF_fracture_bot_panel = (material[mat]['Kc']*10**6/np.sqrt(np.pi*0.005))/abs(sigma_bot_panel_t)
                MOFs[8] = MOF_fracture_bot_panel
                if MOF_sigma_bot_panel_t < 1:
                    print(f"Bottom panel tensile strength failure, MOF: {MOF_sigma_bot_panel_t:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
                if MOF_fracture_bot_panel < 1:
                    print(f"Bottom panel fracture failure, MOF: {MOF_fracture_bot_panel:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break

            if sigma_bot_panel_c < 0:
                MOF_sigma_bot_panel = material[mat]['yield stress']*10**6/abs(sigma_bot_panel_c)
                MOFs[9] = MOF_sigma_bot_panel
                if MOF_sigma_bot_panel < 1:
                    print(f"Bottom panel compressive strength failure, MOF: {MOF_sigma_bot_panel:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break

            if sigma_top_panel_t > 0:
                MOF_sigma_top_panel_t = material[mat]['yield stress']*10**6/abs(sigma_top_panel_t)
                MOFs[10] = MOF_sigma_top_panel_t
                MOF_fracture_top_panel = (material[mat]['Kc']*10**6/np.sqrt(np.pi*0.005))/abs(sigma_top_panel_t)
                MOFs[11] = MOF_fracture_top_panel
                if MOF_sigma_top_panel_t < 1:
                    print(f"Top panel tensile strength failure, MOF: {MOF_sigma_top_panel_t:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break
                if MOF_fracture_top_panel < 1:
                    print(f"Top panel fracture failure, MOF: {MOF_fracture_top_panel:.2f}")
                    failure_ele = True
                    t_f_spar, t_a_spar, t_top, t_bot, a = increase_values([t_f_spar, t_a_spar, t_top, t_bot, a], step_size)
                    b = 2*a
                    break

            if i == spanwise_pos_sec[rib_number][0]:
                w_ribs += A_enc*t_ribs*material[mat]['density']
                if plot == True and run_once_more == 1:
                    # plot provided is how the wingbox will look up until the next cross section
                    plot_wingbox(i)

                # this is done for every section between ribs, values used for failure analysis are the values at the end of the section (most critical)
                ## Define poisson's ratio
                v = 0.33
                ## Spar buckling (shear)
                # a/b values and their respective ks values
                ab_f_spar = (spanwise_pos_sec[rib_number][-1]-spanwise_pos_sec[rib_number][0])/l_f_spar
                if ab_f_spar < 5:
                    ks_f_spar = ks_clamped(ab_f_spar)
                else:
                    ks_f_spar = ks_clamped_min
                ab_a_spar = (spanwise_pos_sec[rib_number][-1]-spanwise_pos_sec[rib_number][0])/l_a_spar
                if ab_a_spar < 5:
                    ks_a_spar = ks_clamped(ab_a_spar)
                else:
                    ks_a_spar = ks_clamped_min

                # critical shear buckling stress
                tau_cr_f_spar = ((np.pi**2*material[mat]['E']*10**9*ks_f_spar)/(12*(1-v**2)))*(t_f_spar/l_f_spar)**2
                tau_cr_a_spar = ((np.pi**2*material[mat]['E']*10**9*ks_a_spar)/(12*(1-v**2)))*(t_a_spar/l_a_spar)**2

                # shear flow
                tau_avg = loading_tension[1][np.where(spanwise_pos == i)[0][0]]/(l_f_spar*t_f_spar+l_a_spar*t_a_spar)
                tau_max = 1.5*tau_avg

                # shear flow due to torsion
                q = loading_tension[3][np.where(spanwise_pos == i)[0][0]]/(2*A_enc)

                # combine | sign convention -> positive shear force -> positive shear stress (pointing downwards), 
                tau_f_spar = tau_max - q/t_f_spar
                tau_a_spar = tau_max + q/t_a_spar

                # check for failure
                MOFs[12] = abs(tau_cr_f_spar/tau_f_spar)
                if abs(tau_f_spar) > tau_cr_f_spar:
                    print(f"Front spar shear buckling, MOF: {abs(tau_cr_f_spar/tau_f_spar):.2f}")
                    failure_ele = True
                    t_f_spar += step_size
                    break
                
                MOFs[13] = abs(tau_cr_a_spar/tau_a_spar)
                if abs(tau_a_spar) > tau_cr_a_spar:
                    print(f"Aft spar shear buckling, MOF: {abs(tau_cr_a_spar/tau_a_spar):.2f}")
                    failure_ele = True
                    t_a_spar += step_size
                    break

                ## Panel buckling (compression) -> look at section between centroids of stringers
                # a/b values and their respective ks values
                ab_top = (spanwise_pos_sec[rib_number][-1]-spanwise_pos_sec[rib_number][0])/(x_stringer_spacing_top**2+y_stringer_spacing_top**2)**0.5
                if ab_top < 5:
                    ks_top = ks_A_clamped(ab_top)
                else:
                    ks_top = ks_A_clamped_min
                ab_bot = (spanwise_pos_sec[rib_number][-1]-spanwise_pos_sec[rib_number][0])/(x_stringer_spacing_bot**2+y_stringer_spacing_bot**2)**0.5
                if ab_bot < 5:
                    ks_bot = ks_A_clamped(ab_bot)
                else:
                    ks_bot = ks_A_clamped_min

                # critical compression buckling stress
                sigma_cr_top = ((np.pi**2*material[mat]['E']*10**9*ks_top)/(12*(1-v**2)))*(t_top/l_top)**2
                sigma_cr_bot = ((np.pi**2*material[mat]['E']*10**9*ks_bot)/(12*(1-v**2)))*(t_bot/l_bot)**2
                # points to evaluate stress at top and bottom
                x_mid_top = np.array([(xc_stringer_top[i-1]+xc_stringer_top[i])/2+(yc_s+t_top/2)*np.sin(angle_top) for i in range(1, n_stringers_top+2)])
                y_mid_top = np.array([(yc_stringer_top[i-1]+yc_stringer_top[i])/2+(yc_s+t_top/2)*np.cos(angle_top) for i in range(1, n_stringers_top+2)])
                x_mid_bot = np.array([(xc_stringer_bot[i-1]+xc_stringer_bot[i])/2+(yc_s+t_bot/2)*np.sin(angle_bot) for i in range(1, n_stringers_bot+2)])
                y_mid_bot = np.array([(yc_stringer_bot[i-1]+yc_stringer_bot[i])/2-(yc_s+t_bot/2)*np.cos(angle_bot) for i in range(1, n_stringers_bot+2)])

                # stress due to bending
                sigma_top_b_p = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(y_mid_top-yc)/Ix
                sigma_bot_b_p = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(y_mid_bot-yc)/Ix

                sigma_top_p = sigma_n + sigma_top_b_p
                sigma_bot_p = sigma_n + sigma_bot_b_p

                # print if wingbox meets buckling requirements, in future wingbox should automatically update if it does not meet requirements
                sigma_top_p_c = sigma_top_p[sigma_top_p < 0]
                sigma_bot_p_c = sigma_bot_p[sigma_bot_p < 0]

                # check for failure
                if len(sigma_top_p_c) != 0:
                    MOFs[14] = abs(sigma_cr_top/abs(sigma_top_p.min()))
                    if not np.all(abs(sigma_top_p_c) < sigma_cr_top):
                        print(f"Panel compression buckling top, MOF: {abs(sigma_cr_top/abs(sigma_top_p.min())):.2f}")
                        failure_ele = True
                        t_top += step_size
                        break
                        
                if len(sigma_bot_p_c) != 0:
                    MOFs[15] = abs(sigma_cr_bot/abs(sigma_bot_p.min()))
                    if not np.all(abs(sigma_bot_p_c) < sigma_cr_bot):
                        print(f"Panel compression buckling bottom, MOF: {abs(sigma_cr_bot/abs(sigma_bot_p.min())):.2f}")
                        failure_ele = True
                        t_bot += step_size
                        break

                ## Stringer column buckling (compression) ## , calculate this once for entire wingbox (stringers run entire length)
                # K factor for end conditions
                K = 1/4
                d_rib = spanwise_pos_sec[rib_number][-1]-spanwise_pos_sec[rib_number][0]
                
                # critical compression buckling stress, assumption: all stringers run full length (will not be true, but it is the worst case)
                sigma_cr_top_s = K*np.pi**2*material[mat]['E']*10**9*((n_stringers_top+2)*Ix_s_t)/((d_rib)**2*(n_stringers_top+2)*A_s)
                sigma_cr_bot_s = K*np.pi**2*material[mat]['E']*10**9*((n_stringers_bot+2)*Ix_s_b)/((d_rib)**2*(n_stringers_bot+2)*A_s)

                # stress due to bending
                sigma_top_b_s = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(np.array(yc_stringer_top)-yc)/Ix
                sigma_bot_b_s = -loading_tension[2][np.where(spanwise_pos == i)[0][0]]*(np.array(yc_stringer_bot)-yc)/Ix

                sigma_top_s = sigma_n + sigma_top_b_s
                sigma_bot_s = sigma_n + sigma_bot_b_s

                # print if stringers meets column buckling requirements, in future wingbox should automatically update if it does not meet requirements
                sigma_top_s_c = sigma_top_s[sigma_top_s < 0]
                sigma_bot_s_c = sigma_bot_s[sigma_bot_s < 0]

                # check for failure
                if len(sigma_top_s_c) != 0:
                    MOFs[16] = abs(sigma_cr_top_s/min(sigma_top_s_c))
                    if MOFs[16] < 1:
                        print(f"Stringer column buckling top, MOF: {abs(sigma_cr_top_s/min(sigma_top_s_c)):.2f}")
                        failure_ele = True
                        inc = (b+step_size)/b
                        b *= inc
                        a *= inc
                        t_stringer *= inc
                        break
                
                if len(sigma_bot_s_c) != 0:
                    MOFs[17] = abs(sigma_cr_bot_s/min(sigma_bot_s_c))
                    if MOFs[17] < 1:
                        print(f"Stringer column buckling bottom, MOF: {abs(sigma_cr_bot_s/max(abs(sigma_bot_s[sigma_bot_s < 0]))):.2f}")
                        failure_ele = True
                        inc = (b+step_size)/b
                        b *= inc
                        a *= inc
                        t_stringer *= inc
                        break

            if failure_ele == True:
                break

        # remove stringers if necessary
        if stringer_decr_x_ribs != 0:
            if (rib_number+1) % stringer_decr_x_ribs == 0:
                n_stringers -= n_stringer_decr_ribs

    if np.all(MOFs > 1) and i == spanwise_pos[-1]:
        if run_once_more == 1:
            failure_tot = False
        elif n_ribs != 0:
            failure_tot = True
            run_once_more += 1
        else:
            failure_tot = False

print(f"t_f_spar: {t_f_spar:.5f} m")
print(f"t_a_spar: {t_a_spar:.5f} m")
print(f"t_top: {t_top:.5f} m")
print(f"t_bot: {t_bot:.5f} m")
print(f"a: {a:.5f} m")
print(f"b: {b:.5f} m")
print(f"t_stringer: {t_stringer:.5f} m")
plot_wingbox(i)

# Calculate tip deflection
d2v_dy2 = -loading_tension[2]/(material[mat]['E']*10**9*Ixs)
dv_dy = cumulative_trapezoid(spanwise_pos, d2v_dy2, initial=0)
v = cumulative_trapezoid(spanwise_pos, dv_dy, initial=0)
max_deflection = max(v)

# Calculate tip twist
dtheta_dy = loading_tension[3]/((material[mat]['E']*10**9/(2*(1+0.33)))*Js)
theta = cumulative_trapezoid(spanwise_pos, dtheta_dy, initial=0)
max_twist = max(theta)

# Calculate weight
weight = trapezoid(As, spanwise_pos)*material[mat]['density'] + w_ribs
# weight = ((As[0]+As[-1])/2)*(ac.b/2)*material[mat]['density']
print(f"Weight: {weight:.2f} kg")

# Plotting
if plot == True:
    plt.figure()
    plt.subplot(221)
    plt.plot(spanwise_pos, Ixs, label='Ix')
    plt.plot(spanwise_pos, Iys, label='Iy')
    plt.plot(spanwise_pos, Js, label='Torsional constant')
    plt.ylabel('[m^4]')
    plt.legend()
    plt.grid()

    plt.subplot(222)
    plt.plot(spanwise_pos, xcs/chord(spanwise_pos), label='xc')
    plt.plot(spanwise_pos, ycs/chord(spanwise_pos), label='yc')
    plt.ylabel('Centroid [m]')
    plt.legend()
    plt.grid()

    plt.subplot(223)
    plt.plot(spanwise_pos, v)
    plt.xlabel('Spanwise location [m]')
    plt.ylabel('Tip deflection [m]')
    plt.grid()
    plt.title(f"Tip deflection: {max_deflection:.6f} m")

    plt.subplot(224)
    plt.plot(spanwise_pos, theta)
    plt.xlabel('Spanwise location [m]')
    plt.ylabel('Tip twist [rad]')
    plt.grid()
    plt.title(f"Tip twist: {max_twist:.6f} rad")
    plt.show()

# checks for failures to be implemented:
# - spar buckling compression
# - panel buckling compression
# - column buckling stringers compression
# - shear buckling spars
# - shear buckling panels
# - compressive strength failure of each component
# - tensile strength failure of each component
# - crack propagation


# list of things to do