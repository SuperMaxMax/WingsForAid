import matplotlib.pyplot as plt
import numpy as np
from parameters import UAV

concept = UAV('naam')

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


def stall_req(obj):
    step = 0.5
    V = np.arange(0, obj.V_D*1.1 + step, step)
    n_pos = 0.5 * 1.225 * V**2 * obj.CL_max_TO / obj.WS               #DENSITY?  #WHICH CL_max?
    n_neg = n_pos * -1

    return V, n_pos, n_neg


def plot(obj):
    V, n_pos, n_neg = stall_req(obj)
    n_max = CS23_max(obj)
    n_min = CS23_min(obj)
   
    #X-Values of points. Points are from typical V-n diagram
    A_x = np.interp(n_max, n_pos, V)
    D_x = obj.V_D
    E_x = D_x
    F_x = obj.V_cruise
    H_x = np.interp(-n_min, n_pos, V)

    #V_S
    V_S = np.interp(1, n_pos, V)
    

    #Setting plot limits:
    left, right = 0, obj.V_D * 1.2
    plt.xlim(left, right)
    bottom, top =  -1 * np.ceil(n_min * -1.2 * 2) / 2 , np.ceil(n_max*1.1 * 2) / 2
    y_range = top - bottom
    plt.ylim(bottom, top)

    yticks = np.arange(bottom, top, 0.5)
    plt.yticks(yticks)
    
    #X-axis:
    plt.axhline(y = 0, color = "black", linewidth = '0.7')


    #Plotting points:
    plt.plot(A_x, n_max, 'ro', color = 'black')    #A
    plt.text(A_x, n_max + 0.35 , s = 'A', ha='center', va='top')

    plt.plot(D_x, n_max, 'ro', color = 'black')    #D
    plt.text(D_x, n_max + 0.35 , s = 'D', ha='center', va='top')

    plt.plot(E_x, 0, 'ro', color = 'black')        #E
    plt.text(E_x + 2,  0 + 0.35 , s = 'E', ha='center', va='top')

    plt.plot(F_x, n_min, 'ro', color = 'black')    #F
    plt.text(F_x, n_min - 0.15 , s = 'F', ha='center', va='top')

    plt.plot(H_x, n_min, 'ro', color = 'black')    #H
    plt.text(H_x, n_min - 0.15 , s = 'H', ha='center', va='top')




    #Plotting straight lines
    plt.plot([A_x, D_x], [n_max, n_max], color = 'black')            #Between A and D
    plt.plot([D_x, E_x], [n_max, 0], color = 'black')    #Between D and E
    plt.plot([E_x, F_x], [0, n_min], color = 'black')    #Between E and F
    plt.plot([F_x, H_x], [n_min, n_min], color = 'black')    #Between F and H


    #Plot stall speed requirement until point A for positive and until point H for negative:
    V_pos = V[V<=A_x]
    n_pos = n_pos[0:V_pos.size]
    V_neg = V[V<=H_x]
    n_neg = n_neg[0:V_neg.size]

    plt.plot(V_pos, n_pos, color = 'black')
    plt.plot(V_neg, n_neg, color = 'black')

    #Plot V_S
    y0_frac_VS = abs(bottom) / y_range
    y1_frac_VS = (abs(bottom) + 1) / y_range
    plt.axvline(x = V_S, ymin = y0_frac_VS, ymax = y1_frac_VS, linestyle = "--", linewidth = "1", color = "black")
    

    #Plot V_A
    y0_frac_VA = (abs(bottom) + n_max) / y_range
    y1_frac_VA = (abs(bottom)) / y_range
    plt.axvline(x = A_x, ymin = y0_frac_VA, ymax = y1_frac_VA, linestyle = "--", linewidth = "1", color = "black")
    
    #Plot V_C
    y0_frac_VA = (abs(bottom) + n_min) / y_range
    y1_frac_VA = (abs(bottom)) / y_range
    plt.axvline(x = F_x, ymin = y0_frac_VA, ymax = y1_frac_VA, linestyle = "--", linewidth = "1", color = "black")
    

    #Plot V_C



    #Plot y=1
    plt.axhline(y = 1, linestyle = "--", linewidth = "1", color = "black")



    plt.show()


plot(concept)
 
