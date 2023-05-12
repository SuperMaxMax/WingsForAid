import  numpy as np
from    parameters import *
import  matplotlib.pyplot as plt

def metertofeet(meter):
    feet = meter*3.2808399
    return feet

# W/S and W/P diagrams
# drag polar
def dragpolar(CL, CD0, e, A):
    CD = CD0 + (CL**2/(np.pi*A*e))
    return CD

def stallWS(V, rho, CL_max):
    WoS = 1/2 * rho * V**2 * CL_max
    return WoS

def TOP_calc(WS, sigma, CL_TO, BHP, W_TO):
    TOP = WS/(sigma*CL_TO*(BHP/W_TO))
    return TOP

def altitude_effects(h, Lambda, R, g, T0, rho0, BHP0):
    rho     = rho0*(1+ ((Lambda*h)/T0))**(-((g/(R*Lambda))+1))
    sigma   = rho/rho0
    BHP     = BHP0*(sigma)**(3/4)
    return rho, sigma, BHP

def create_line(x1, y1, x2, y2, num_points):
    x = np.linspace(x1, x2, num_points)
    y = np.linspace(y1, y2, num_points)
    line = np.vstack((x, y))
    return line

def WP_WSdiagrams(h, plot=True, CLrange=True):
    design_points = np.empty(0)
    for i in range(7):
        CL_max_clean    = 1.3+(0.1*i)   #1.3 - 1.9
        CL_max_TO       = 1.3+(0.1*i)
        CL_max_land     = 1.6+(0.1*i)
        CL_TO           = CL_max_TO/(1.1**2)
        CL_LDG          = CL_max_land/(1.1**2)
        # find atmospheric properties and power at altitude according to ISA
        rho, sigma, BHP = altitude_effects(h, Lambda, R, g0, T0, rho0, P_max)
        
        # define a W/S and W/P array to make graphs later on
        WS  = np.arange(100.0, 2201.0, 1.0)
        WP  = np.arange(0.0, 1.501, (1.5/len(WS)))
        
        # sizing for stall, maximum W/S values based on a CL max and a Vstall of 50 kts for dropping
        WS_stall = stallWS(V_s_min, rho, CL_max_clean)
        WS_stall = np.full(len(WP), WS_stall)
        WS_stall = np.vstack((WS_stall, WP))
        
        # Take off sizing using Take-off parameter (TOP) 
        TOP_req = 250 
        WP_TO = (TOP_req/WS)*CL_TO*sigma                                     #read from graph of ADSEE 1 slides, slide 29 lecture 3
        WP_TO = np.vstack((WP_TO, WS))

        # Landing W/S
        WS_landing = (CL_LDG*rho*(LDG_dist/0.5915))/(2*f)
        WS_landing = np.full(len(WP), WS_landing)
        WS_landing = np.vstack((WS_landing, WP))
        
        # W/P cruise
        rho_cruise, sigma_cruise, BHP_cruise = altitude_effects(h_cruise, Lambda, R, g0, T0, rho0, P_max)
        WP_cruise   = (power_setting/cruise_frac)*eta_p*sigma**(3/4)*(((CD0*1/2*rho_cruise*V_cruise**3)/(WS*cruise_frac))+(WS*(1/(np.pi*A*e*rho_cruise*V_cruise))))**(-1)
        WP_cruise   = np.vstack((WP_cruise, WS))

        #Climb performance
        climb_rate  = (8.3/100)*V_climb                      #8.3% climb gradient found in CS23
        WP_Climb    = eta_p/(climb_rate+((np.sqrt(WS)*np.sqrt(2/rho0))/(1.345*(A*e)**(3/4)/CD0**(1/4))))
        WP_Climb    = np.vstack((WP_Climb, WS))
        
        #Evaluating if stall WS or landing WS is limiting
        stall_limit     = 0
        landing_limit   = 0
        plot_stall_WS   = False
        plot_both_WS    = False
        plot_LDG_WS     = False

        if WS_stall[0][0] < WS_landing[0][0]:
            plot_stall_WS = True
            WS_limit    = WS_stall
        if WS_stall[0][0] == WS_landing[0][0]:
            plot_both_WS = True
            WS_limit    = WS_stall
        if WS_stall[0][0] > WS_landing[0][0]:
            plot_LDG_WS = True
            WS_limit    = WS_landing
        
        #find the design point
        design_point = None
        WS_designpoint = WS_limit[0][0:2101]
        tolerance = 0.5
        index = np.where(np.abs(WP_TO[1] - WS_designpoint) <= tolerance)
        WS_designpoint = WS_designpoint[0]
        WP_values = np.array([WP_TO[0][index], WP_Climb[0][index], WP_cruise[0][index]])
        WP_designpoint = np.min(WP_values)
        design_point = (WS_designpoint, WP_designpoint)
        design_points = np.append(design_points, design_point)

        #plotting
        if plot:
            lab = "Take-off, CL TO = " + str(np.round(CL_TO, decimals=2))
            plt.plot(WP_TO[1], WP_TO[0], label=lab)
            if plot_stall_WS:
                plt.plot(WS_stall[0], WS_stall[1])
            if plot_both_WS:
                plt.plot(WS_landing[0], WS_landing[1])
                plt.plot(WS_stall[0], WS_stall[1])
            if plot_LDG_WS:
                plt.plot(WS_landing[0], WS_landing[1])
            plt.plot(WP_cruise[1], WP_cruise[0], label="Cruise requirement")
            plt.plot(WP_Climb[1], WP_Climb[0], label="Climb requirement")
            if design_point is not None:
                plt.plot(design_point[0], design_point[1], 'ro', label='Design Point')
    
    if plot:
        plt.xlabel("W/S [N/m^2]")
        plt.ylabel("W/P [N/W]")
        plt.ylim((0, 1.5))
        #plt.legend(loc='upper right')
        plt.show()
    
    return design_points

#W/P and W/S are now known, define wing surface area
def geometry_determination(MTOW, plot=True):
    WPWS = WP_WSdiagrams(0, plot=False)
    W_TO = MTOW*g0
    WS_values = np.array(WPWS[0:12:2])
    WP_values = np.array(WPWS[1:13:2])
    Sw = W_TO/WS_values
    P_values = W_TO/WP_values*0.00134102209                 #convert to horsepower
    b = np.sqrt(A*Sw)
    cos_lambda_c04 = 1
    lambda_co4 = np.arccos(cos_lambda_c04)                  #rad, As Mcruise < 0.7, use 0 sweep angle
    taper = 0.2*(2-lambda_co4)
    rootchord = (2*Sw)/((1+taper)*b)
    tipchord = taper*rootchord
    wings = np.empty(0)
    for i in range(len(Sw)):
        points = np.array([[0, rootchord[i]/4, tipchord[i]/4, 0, -tipchord[i]/4, -3*tipchord[i]/4, -3*rootchord[i]/4, -rootchord[i]/4],
                            [0, 0, b[i]/2, b[i]/2, b[i]/2, b[i]/2, 0, 0]])
        LE = create_line(points[0][1], points[1][1], points[0][2], points[1][2], 1000) #leading edge
        ct = create_line(points[0][2], points[1][2], points[0][4], points[1][4], 1000)  #tip chord
        TE = create_line(points[0][4], points[1][4], points[0][6], points[1][6], 1000) #trailing edge
        cr = create_line(points[0][6], points[1][6], points[0][1], points[1][1], 1000)  #root chord
        qc = create_line(points[0][0], points[1][0], points[0][3], points[1][3], 1000) #quarter chord line
        hc = create_line(points[0][-1], points[1][-1], points[0][4], points[1][4], 1000) #half chord line
        #points used to create MAC
        point_tip = (points[0][2]+rootchord[i], points[1][2])
        point_root= (points[0][6]-tipchord[i], points[1][6])
        constr_line = create_line(point_root[0], point_root[1], point_tip[0], point_tip[1], 1000)
        mask = np.abs(constr_line - hc)
        mask = mask[0] + mask[1]
        min = np.min(mask)
        index = np.where(mask == min)
        y_mac = hc[1][index]
        tolerance = 0.0000005
        x_lemac = LE[0][np.where(np.abs(LE[1]-y_mac)<=tolerance)][0]
        x_temac = TE[0][np.where(np.abs(TE[1]-y_mac)<=tolerance)][0]
        MAC = np.array([np.linspace(x_temac, x_lemac, 1000), np.full(1000, y_mac)])
        if plot:
            plt.plot(LE[0], LE[1], color='black')
            plt.plot(ct[0], ct[1], color='black')
            plt.plot(TE[0], TE[1], color='black')
            plt.plot(cr[0], cr[1], color='black')
            plt.plot(qc[0], qc[1], color='black')
            plt.plot(hc[0], hc[1], color='black')
            #plt.plot(constr_line[0], constr_line[1], color='red')
            plt.plot(MAC[0], MAC[1], color='black')
            plt.gca().set_aspect('equal', adjustable = 'box')
            plt.xlabel("y [m]")
            plt.show()
        
    
        
        
    return Sw, P_values, b, lambda_co4, taper, rootchord, tipchord, 

a=geometry_determination(W_TO)





    









