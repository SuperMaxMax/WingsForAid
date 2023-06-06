import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt

def plot_wingbox():
    """Plots the wingbox of the aircraft"""
    # Plot airfoil shape
    plt.plot(x, y_u, 'b')
    plt.plot(x, y_l, 'b')
    # Plot spars and top and bottom of wingbox
    plt.scatter(xp, yp)
    plt.plot([xp[0], xp[1]], [yp[0], yp[1]], 'r')
    plt.plot([xp[1], xp[2]], [yp[1], yp[2]], 'r')
    plt.plot([xp[2], xp[3]], [yp[2], yp[3]], 'r')
    plt.plot([xp[3], xp[0]], [yp[3], yp[0]], 'r')
    # Plot centroid
    plt.scatter(xc, yc, c='g')
    plt.axis('equal')
    plt.grid()
    plt.show()

def chord(y):
    return ac.rootchord-((ac.rootchord-ac.tipchord)/(ac.b/2))*y

def airfoil(f_spar, a_spar, span):
    """Returns the upper and lower surface of an airfoil with a common x-axis"""
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

ac = UAV("aircraft")
airfoil_num = '4415'

# Input: define wingbox parameters
f_spar = 0.25
t_f_spar = 0.01
a_spar = 0.75
t_a_spar = 0.01

t_top = 0.01
t_bottom = 0.01

n_stringers = 10
A_stringer = 0.0001
I_stringer = 0.0001

n_ribs = 10
t_ribs = 0.01
spacing_ribs = 0.1

# xp and yp -> [top left, top right, bottom right, bottom left]
x, y_u, y_l, xp, yp = airfoil(f_spar, a_spar, 0)


# Calculate centroid, zero point coincides with origin plot
# Spar contributions
l_f_spar = yp[0]-yp[3]
l_a_spar = yp[1]-yp[2]
yc_f_spar = l_f_spar*t_f_spar*(yp[3]+l_f_spar/2)
xc_f_spar = l_f_spar*t_f_spar*(xp[0])
yc_a_spar = l_a_spar*t_a_spar*(yp[2]+l_a_spar/2)
xc_a_spar = l_a_spar*t_a_spar*(xp[1])
# Skin top and bottom contributions
l_top = ((xp[1]-xp[0])**2+(yp[1]-yp[0])**2)**0.5
l_bottom = ((xp[2]-xp[3])**2+(yp[2]-yp[3])**2)**0.5
xc_top = l_top*t_top*(xp[0]+(xp[1]-xp[0])/2)
yc_top = l_top*t_top*(yp[0]+(yp[1]-yp[0])/2)
xc_bottom = l_bottom*t_bottom*(xp[3]+(xp[2]-xp[3])/2)
yc_bottom = l_bottom*t_bottom*(yp[3]+(yp[2]-yp[3])/2)
# Stringer contributions
# to come
# Rib contributions
# to come
# Final centroid
xc = (xc_f_spar+xc_a_spar+xc_top+xc_bottom)/((l_f_spar*t_f_spar)+(l_a_spar*t_a_spar)+(l_top*t_top)+(l_bottom*t_bottom))
yc = (yc_f_spar+yc_a_spar+yc_top+yc_bottom)/((l_f_spar*t_f_spar)+(l_a_spar*t_a_spar)+(l_top*t_top)+(l_bottom*t_bottom))

plot_wingbox()


