import matplotlib.pyplot as plt
import numpy as np
from parameters import UAV


def VC_lim_low(obj):
    return 33 * (obj.WS * 0.020885) ** 0.5 * 0.51444 #[m/s]

def VD_lim_low(obj):
    return 1.5 * VC_lim_low(obj)    #[m/s]

def VA_lim_low(obj): #Cant be greater than VC
    V, n_pos, n_neg = stall_req(obj)
    V_S = np.interp(1, n_pos, V)

    n = CS23_max(obj)

    return V_S * n**0.5 #[m/s]

def VS1(obj):
    return ((2 * obj.WS) / (1.225 * obj.CL_max_clean )) ** 0.5


def stall_req(obj):
    step = 0.5
    V_D = VD_lim_low(obj)
    V = np.arange(0, V_D*1.1 + step, step)
    n_pos = 0.5 * 1.225 * V**2 * obj.CL_max_land / obj.WS               #DENSITY?  #WHICH CL_max?
    n_neg = n_pos * -1

    return V, n_pos, n_neg


def VB_lim_low(obj): #Cant be greater than VC

    V = VC_lim_low(obj)              # DEFINE SPEEDS
    u_hat_fs = 50              # [f/s]  
    u_hat_ms = u_hat_fs * 0.3048                                  
    rho_cruise = 0.9
    MGC = 13 / 10
    mu = 2 * obj.WS / (rho_cruise * obj.g0 * MGC * obj.CLa)  # Airplane mass ratio []

    K = 0.88 * mu / (5.3 + mu)                              

    u = K * u_hat_ms

    delta_n = obj.rho0 * V * obj.CLa * u / (2 * obj.WS)
    
    n_c = 1 + delta_n

    V, n_pos, n_neg = stall_req(obj)
    V_S1 = VS1(obj)
    
    VB1 = V_S1 * n_c**0.5 #[m/s]
    
    V, n_pos, n_neg = stall_req(obj)
    
    
          # DEFINE SPEEDS

    u_hat_fs = 66              # [f/s]  
    u_hat_ms = u_hat_fs * 0.3048                                  
    rho_cruise = 0.9
    MGC = 13 / 10
    mu = 2 * obj.WS / (rho_cruise * obj.g0 * MGC * obj.CLa)  # Airplane mass ratio []

    K = 0.88 * mu / (5.3 + mu)                              

    u = K * u_hat_ms

    delta_n = obj.rho0 * VB1 * obj.CLa * u / (2 * obj.WS)
    
    n_b = 1 + delta_n

    y2 = n_b
    y1 = 1
    x2 = VB1
    n_line = (y2-y1) / x2 * V + 1

    VB2 = 0
    for i in range(len(V)):
        if n_line[i] - 0.1 < n_pos[i] < n_line[i] + 0.1:
            VB2 = V[i]

    return min(VB1, VB2)


def CS23_max(obj):
    if obj.type == "normal":
        y = min(2.1 + (24000 / (obj.W_TO * 2.205 + 10000)), 3.8)
    elif obj.type == "utility":
        y = 4.4
    else:
        print("Aircraft not properly defined in parameters. Don't take results seriously!")
        y = 0

    return y

def CS23_min(obj):
    y = -0.4 * CS23_max(obj)

    return y




def plot_Vn(obj):
    V, n_pos, n_neg = stall_req(obj)
    n_max = CS23_max(obj)
    n_min = CS23_min(obj)
   
    # X-Values of points. Points are from typical V-n diagram
    A_x = np.interp(n_max, n_pos, V)
    D_x = VD_lim_low(obj)
    E_x = D_x
    F_x = VC_lim_low(obj)
    H_x = np.interp(-n_min, n_pos, V)

    # V_S
    V_S = np.interp(1, n_pos, V)
    
    # Setting plot limits:
    left, right = 0, obj.V_D * 1.2
    plt.xlim(left, right)
    bottom, top = -1 * np.ceil(n_min * -1.2 * 2) / 2, np.ceil(n_max*1.1 * 2) / 2
    y_range = top - bottom
    plt.ylim(bottom, top)

    yticks = np.arange(bottom, top, 0.5)
    plt.yticks(yticks)
    
    # X-axis:
    plt.axhline(y = 0, color = "black", linewidth = '0.7')
    plt.ylabel("Load factor (n)")
    plt.xlabel("Airspeed (V) [m/s]")

    # Plotting points:
    plt.plot(A_x, n_max, 'ko')      #A
    plt.text(A_x, n_max + 0.35 , s = 'A', ha='center', va='top')

    plt.plot(D_x, n_max, 'ko')      #D
    plt.text(D_x, n_max + 0.35 , s = 'D', ha='center', va='top')

    plt.plot(E_x, n_min, 'ko')          #E
    plt.text(E_x + 2,  n_min + 0.35 , s = 'E', ha='center', va='top')

    plt.plot(F_x, n_min, 'ko')      #F
    plt.text(F_x, n_min - 0.15 , s = 'F', ha='center', va='top')

    plt.plot(H_x, n_min, 'ko')      #H
    plt.text(H_x, n_min - 0.15 , s = 'H', ha='center', va='top')

    # Plotting straight lines
    plt.plot([A_x, D_x], [n_max, n_max], color = 'black')   # Between A and D
    plt.plot([D_x, E_x], [n_max, n_min], color = 'black')       # Between D and E
    plt.plot([E_x, F_x], [n_min, n_min], color = 'black')       # Between E and F
    plt.plot([F_x, H_x], [n_min, n_min], color = 'black')   # Between F and H


    # Plot stall speed requirement until point A for positive and until point H for negative:
    V_pos = V[V<=A_x]
    n_pos = n_pos[0:V_pos.size]
    V_neg = V[V<=H_x]
    n_neg = n_neg[0:V_neg.size]

    plt.plot(V_pos, n_pos, color = 'black')
    plt.plot(V_neg, n_neg, color = 'black')

    # Plot V_S
    y0_frac_VS = abs(bottom) / y_range
    y1_frac_VS = (abs(bottom) + 1) / y_range

    plt.axvline(x = V_S, ymin = y0_frac_VS, ymax = y1_frac_VS, linestyle = "--", linewidth = "1", color = "black")
    plt.text(V_S, -0.15 , s = 'Vs', ha='center', va='top')

    # Plot V_A
    y0_frac_VA = (abs(bottom) + n_max) / y_range
    y1_frac_VA = (abs(bottom)) / y_range
    plt.axvline(x = A_x, ymin = y0_frac_VA, ymax = y1_frac_VA, linestyle = "--", linewidth = "1", color = "black")
    plt.text(A_x, -0.15 , s = 'Va', ha='center', va='top')

    # Plot V_C
    y0_frac_VA = (abs(bottom) + n_min) / y_range
    y1_frac_VA = (abs(bottom)) / y_range
    plt.axvline(x = F_x, ymin = y0_frac_VA, ymax = y1_frac_VA, linestyle = "--", linewidth = "1", color = "black")
    plt.text(F_x, 0.3 , s = 'Vc', ha='center', va='top')
    plt.text(E_x - 3, 0.3 , s = 'Vd', ha='center', va='top')



    # Plot y=1
    plt.axhline(y = 1, linestyle = "--", linewidth = "1", color = "black")

    return n_max

def gust_points(obj):
    # ORDER: VB pos, VB neg, VC pos, VC, VD pos, VD neg
    V = np.array([VB_lim_low(obj), VB_lim_low(obj), VC_lim_low(obj), VC_lim_low(obj), VD_lim_low(obj), VD_lim_low(obj)]) # DEFINE SPEEDS
    u_hat_fs = np.array([66, 66, 50, 50, 25, 25])   # [f/s]  
    u_hat_ms = u_hat_fs * 0.3048                                  
    rho_cruise = 0.9
    MGC = 13 / 10
    mu = 2 * obj.WS / (rho_cruise * obj.g0 * MGC * obj.CLa)  # Airplane mass ratio []

    K = 0.88 * mu / (5.3 + mu)                              

    u = K * u_hat_ms

    delta_n_mag = obj.rho0 * V * obj.CLa * u / (2 * obj.WS)
    delta_n = np.array([1, -1, 1, -1, 1, -1]) * delta_n_mag
    
    n_peak = 1 + delta_n

    return V, n_peak

def max_n(obj):

    V, n_peak = gust_points(obj)
    n_peak = np.append(n_peak, CS23_max(obj))

    return max(n_peak)

def plot_gust(obj):

    V, n_peak = gust_points(obj)
    labels = ["B'","G'","C'","F'","D'","E'"]

    for i in range(len(V)):
        plt.plot(V[i], n_peak[i], 'o', color = 'pink')
        if labels[i]!="F'" and labels[i]!="D'" and labels[i]!="E'" and labels[i]!="B'":
            plt.text(V[i], n_peak[i] + 0.35 , s = labels[i], ha='center', va='top')
        elif labels[i]=="B'":
            plt.text(V[i]-3, n_peak[i] + 0.35 , s = labels[i], ha='center', va='top')
        else:
            plt.text(V[i]+3, n_peak[i] + 0.35 , s = labels[i], ha='center', va='top')
        

    # Plotting straight lines
    plt.plot([0, V[0]], [1, n_peak[0]], color = 'pink')                        #Between 1 and B'
    plt.plot([0, V[1]], [1, n_peak[1]], color = 'pink')                        #Between 1 and G'
    plt.plot([V[0], V[2]], [n_peak[0], n_peak[2]], color = 'pink')             #Between B' and C'
    plt.plot([V[1], V[3]], [n_peak[1], n_peak[3]], color = 'pink')             #Between G' and F'
    plt.plot([V[2], V[4]], [n_peak[2], n_peak[4]], color = 'pink')             #Between C' and D'
    plt.plot([V[3], V[5]], [n_peak[3], n_peak[5]], color = 'pink')             #Between F' and E'
    plt.plot([V[4], V[5]], [n_peak[4], n_peak[5]], color = 'pink')             #Between D' and E'


    n_max = CS23_max(obj)
    n_min = CS23_min(obj)

    # Plot V_C
    bottom, top = -1 * np.ceil(n_min * -1.2 * 2) / 2, np.ceil(n_max*1.1 * 2) / 2
    y_range = top - bottom
    y0_frac_VA = (abs(bottom)) / y_range
    y1_frac_VA = (abs(bottom) + n_peak[0]) / y_range
    
    plt.axvline(x = V[0], ymin = y0_frac_VA, ymax = y1_frac_VA, linestyle = "--", linewidth = "1", color = "black")
    plt.text(V[0] - 3, 0.3 , s = 'Vb', ha='center', va='top')

    plt.axhline(y = max_n(obj), linestyle = "--", linewidth = "1.5", color = "red")
    plt.text(10, max_n(obj) + 0.3 , s = "max n", ha='center', va='top', color = 'red', fontsize = 12)   
    

    

def plot_all(obj):
    V, n_pos, n_neg = stall_req(obj)
   
    V_S = np.interp(1, n_pos, V)

    obj.V_s_min             = V_S
    obj.V_cruise            = VC_lim_low(obj)
    obj.V_D                 = VD_lim_low(obj)
    obj.V_B                 = VB_lim_low(obj)
    obj.V_A                 = VA_lim_low(obj)

    plot_Vn(obj)
    plot_gust(obj)
    plt.title(f"V-n diagram for {obj.name}")
#     plt.show()

# concept = UAV('DET_CON_1', 'tractor', boom=False, braced_wing=False)
# plot_all(concept)
