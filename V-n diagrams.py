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
    print(n_min)
    print(n_neg)
    print(V)
   
    #X-Values of points
    A = np.interp(n_max, n_pos, V)
    D = obj.V_D
    E = D
    F = obj.V_cruise
    H = np.interp(n_min, n_neg, V)
    print(H)

    #V_S
    V_S = np.interp(1, n_pos, V)
    

    #Setting plot limits:
    left, right = 0, obj.V_D * 1.2
   # plt.xlim(left, right)
    bottom, top = n_min * 1.2, n_max * 1.2
    #plt.ylim(bottom, top)

    #Plotting points:
    plt.plot(A, n_max, 'ro')    #A
    plt.plot(D, n_max, 'ro')    #D
    plt.plot(D, 0, 'ro')        #E
    plt.plot(F, n_min, 'ro')    #F
    plt.plot(H, n_min, 'ro')    #H



    #Plotting straight lines
    plt.plot([A, D], [n_max, n_max])            #Between A and D
    plt.plot([D, E], [n_max, 0])    #Between D and E
    plt.plot([E, F], [0, n_min])    #Between E and F



    plt.axvline(x = V_S, ymin = -1 * n_min / bottom)
    plt.axhline(y = 1)



    plt.plot(V, n_pos)
    plt.plot(V, n_neg)
    plt.show()



plot(concept)
 
