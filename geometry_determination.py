import  numpy as np
from    parameters import *
import  matplotlib.pyplot as plt
import  pandas as pd

# def metertofeet(meter):
#     feet = meter*3.2808399
#     return feet

# def dragpolar(obj):
#     CD = obj.CD0 + (obj.CL**2/(np.pi*obj.A*obj.e))
#     return CD

def stallWS(obj, CL_max, h):
    rho = altitude_effects(obj, h)[0]
    WoS = 1/2 * rho * obj.V_s_min**2 * CL_max
    return WoS

def altitude_effects(obj, h):
    rho     = obj.rho0*(1 + ((obj.Lambda*h)/obj.T0))**(-((obj.g0/(obj.R_gas*obj.Lambda))+1))
    sigma   = rho/obj.rho0
    BHP     = obj.power*(sigma)**(3/4)
    return rho, sigma, BHP

def create_line(x1, y1, x2, y2, num_points):
    x = np.linspace(x1, x2, num_points)
    y = np.linspace(y1, y2, num_points)
    line = np.vstack((x, y))
    return line

def geometry_determination(obj, plot=False):
    # create empty array for capturing design points from W/P - W/S diagrams
    design_points = np.empty(0)
    for i in range(len(obj.CL_max_clean)):                  
        CL_max_clean    = obj.CL_max_clean[i]   
        CL_max_TO       = obj.CL_max_TO[i]
        CL_max_land     = obj.CL_max_land[i]
        CL_TO           = obj.CL_TO[i]
        CL_LDG          = obj.CL_LDG[i]
        
        # find atmospheric properties and power at altitude according to ISA
        rho_TO, sigma_TO, BHP_TO = altitude_effects(obj, obj.h_TO)
        obj.rho_TO = rho_TO
        obj.sigma_TO = sigma_TO
        obj.BHP_TO = BHP_TO
        
        # define a W/S and W/P array to make graphs later on
        WS  = np.arange(100.0, 2201.0, 1.0)
        WP  = np.arange(0.0, 1.501, (1.5/len(WS)))
        
        # sizing for stall, maximum W/S values based on a CL max and a Vstall of 50 kts for dropping
        WS_stall = stallWS(obj, CL_max_clean, obj.h_TO)
        WS_stall = np.full(len(WP), WS_stall)
        WS_stall = np.vstack((WS_stall, WP))
        
        # take off sizing using Take-off parameter (TOP) 
        obj.TOP_req = 250                               # read from graph of ADSEE 1 slides, slide 29 lecture 3
        WP_TO = (obj.TOP_req/WS)*CL_TO*obj.sigma_TO         
        WP_TO = np.vstack((WP_TO, WS))

        # landing W/S
        WS_landing = (CL_LDG*obj.rho_TO*(obj.LDG_dist/0.5915))/(2*obj.Mff)
        WS_landing = np.full(len(WP), WS_landing)
        WS_landing = np.vstack((WS_landing, WP))
        
        # W/P cruise
        rho_cruise, sigma_cruise, BHP_cruise = altitude_effects(obj, obj.h_cruise)
        obj.rho_cruise = rho_cruise
        obj.sigma_cruise = sigma_cruise
        obj.BHP_cruise = BHP_cruise
        WP_cruise   = (obj.power_setting/obj.cruise_frac)*obj.eta_p*obj.sigma_cruise**(3/4)*(((obj.CD0*1/2*obj.rho_cruise*obj.V_cruise**3)/(WS*obj.cruise_frac))+(WS*(1/(np.pi * obj.A * obj.e* obj.rho_cruise * obj.V_cruise))))**(-1)
        WP_cruise   = np.vstack((WP_cruise, WS))

        # climb performance
        obj.climb_rate  = (8.3/100) * obj.V_climb    #8.3% climb gradient found in CS23
        WP_Climb    = obj.eta_p/(obj.climb_rate+((np.sqrt(WS)*np.sqrt(2/obj.rho0))/(1.345*(obj.A*obj.e)**(3/4)/obj.CD0**(1/4))))
        WP_Climb    = np.vstack((WP_Climb, WS))
        
        # evaluating if stall WS or landing WS is limiting
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
        
        # find the design point
        design_point = None
        WS_designpoint = WS_limit[0][0:2101]
        tolerance = 0.5
        index = np.where(np.abs(WP_TO[1] - WS_designpoint) <= tolerance)
        WS_designpoint = WS_designpoint[0]
        WP_values = np.array([WP_TO[0][index], WP_Climb[0][index], WP_cruise[0][index]])
        WP_designpoint = np.min(WP_values)
        design_point = (WS_designpoint, WP_designpoint)
        design_points = np.append(design_points, design_point)
        design_point_high = None
        WS_intersection=WP_Climb[0][np.argmin(np.abs(WP_Climb[0]-WP_cruise[0]))]
        print(WS_intersection)
        


        # plotting
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
            plt.plot(design_point_high[0], design_point_high[1], 'ro', label='High W/S Design Point')
    
    if plot:
        plt.xlabel("W/S [N/m^2]")
        plt.xlim((300, 1000))
        plt.ylabel("W/P [N/W]")
        plt.ylim((0, 0.7))
        # plt.legend(loc='upper right')
        plt.show()
    
    # --- Fuselage parameters
    cumulative_box_length   = obj.n_boxes*0.4                   # box 40x40x60, cumulative length in meter
    length_between_boxes    = ((obj.n_boxes/2)-1)*0.2           # 20 cm in between boxes
    engine_length           = 0.6651                            # engine length in cm, EASA type certificate data sheet ROTAX 912 series
    engine_fairing          = 0.2                               # 20 cm room around the engine 
    obj.d_eff               = np.sqrt(1.10*1.40)                # from cross sectional drawing with width 1.40 m and height 1.10 m
    d_engine_boxes          = 0.4                               # 40 cm, leaves room for possible fire wall
    if obj.boom:                                                # assume one effective diameter after last box
        l_tc = 0.8                                              # from drawing [m]
    else:                                                       
        l_tc = 1.75*obj.d_eff                                   # Based on drawing of the aircraft.

    # Fuselage dimensions
    obj.h_out = 1.10                                            # meter, from cross sectional drawing
    obj.h_in  = 0.90
    obj.w_out = 1.40
    obj.w_in  = 1.20
    obj.l_f   = cumulative_box_length + length_between_boxes + engine_length + engine_fairing + l_tc + d_engine_boxes
    if obj.boom:
        obj.S_G = 19.77
    else:
        obj.S_G = 35.66
    # --- Wing parameters
    Weight_TO = obj.W_TO*obj.g0                                 # find the take off weight in newtons
    WS_values = design_points[0:12:2]                           # take the S/W values
    WP_values = design_points[1:13:2]                           # take the P/W values
    obj.WS    = WS_values[0]
    obj.WP    = WP_values[0]
    
    # find surface area
    obj.Sw = Weight_TO/WS_values
    
    
    # find power value                 
    obj.P_values = Weight_TO/WP_values*0.00134102209            #convert to horsepower
    
    # wingspan
    obj.b = np.sqrt(obj.A*obj.Sw)
    

    obj.MGC = obj.Sw/obj.b                                      # Mean geometric chord [m]
    
    # quarter chord sweep angle (0 as the cruise speed is around 100-110 knots which equates to M<0.2)
    obj.cos_lambda_c04 = 1
    obj.lambda_co4 = np.arccos(obj.cos_lambda_c04)              #rad, As Mcruise < 0.7, use 0 sweep angle
    
    # taper ratio
    obj.taper = 0.2*(2-obj.lambda_co4)
    
    # root and tipchord
    obj.rootchord = (2*obj.Sw)/((1+obj.taper)*obj.b)
    obj.tipchord = obj.taper*obj.rootchord

    for i in range(len(obj.Sw)):
        # define important points on the wing planform (corners, quarther and half chord points)
        points = np.array([[0, obj.rootchord[i]/4, obj.tipchord[i]/4, 0, -obj.tipchord[i]/4, -3*obj.tipchord[i]/4, -3*obj.rootchord[i]/4, -obj.rootchord[i]/4],
                           [0, 0, obj.b[i]/2, obj.b[i]/2, obj.b[i]/2, obj.b[i]/2, 0, 0]])
        LE = create_line(points[0][1], points[1][1], points[0][2], points[1][2], 1000)      # leading edge
        ct = create_line(points[0][2], points[1][2], points[0][4], points[1][4], 1000)      # tip chord
        TE = create_line(points[0][4], points[1][4], points[0][6], points[1][6], 1000)      # trailing edge
        cr = create_line(points[0][6], points[1][6], points[0][1], points[1][1], 1000)      # root chord
        qc = create_line(points[0][0], points[1][0], points[0][3], points[1][3], 1000)      # quarter chord line
        hc = create_line(points[0][-1], points[1][-1], points[0][4], points[1][4], 1000)    # half chord line
        
        # points used to create MAC
        point_tip = (points[0][2]+obj.rootchord[i], points[1][2])
        point_root= (points[0][6]-obj.tipchord[i], points[1][6])
        constr_line = create_line(point_root[0], point_root[1], point_tip[0], point_tip[1], 1000)
        mask = np.abs(constr_line - hc)
        mask = mask[0] + mask[1]
        min = np.min(mask)
        index = np.where(mask == min)
        obj.y_mac = hc[1][index][0]
        tolerance = 0.01
        obj.x_lemac = LE[0][np.where(np.abs(LE[1]-obj.y_mac)<=tolerance)][0]
        x_temac = TE[0][np.where(np.abs(TE[1]-obj.y_mac)<=tolerance)][0]
        MAC = np.array([np.linspace(x_temac, obj.x_lemac, 1000), np.full(1000, obj.y_mac)])
        obj.MAC_length = obj.x_lemac-x_temac
        
        # half chord sweep angle
        tan_lambda_co2 = (points[0][-1]-points[0][4])/(points[1][4]-points[1][-1])
        obj.lambda_co2 = np.arctan(tan_lambda_co2)

        # plot the wings
        if plot:
            plt.plot(LE[0], LE[1], color='black')
            plt.plot(ct[0], ct[1], color='black')
            plt.plot(TE[0], TE[1], color='black')
            plt.plot(cr[0], cr[1], color='black')
            plt.plot(qc[0], qc[1], color='black')
            plt.plot(hc[0], hc[1], color='black')
            # plt.plot(constr_line[0], constr_line[1], color='red')
            plt.plot(MAC[0], MAC[1], color='black')
            plt.gca().set_aspect('equal', adjustable = 'box')
            plt.xlabel("y [m]")
            plt.show()
        
            
    obj.t_c = np.full(np.shape(obj.Sw), 0.12)                                    # thickness over chord. 0.18 ADSEE 1, Lecture 6, Slide 22; no supercritical airfoils considered
    obj.dihedral = np.full(np.shape(obj.Sw), 1)                                  # degree, high wing value, same slides as above.

    obj.Qw_wing = ((obj.kq*obj.t_c)/(np.sqrt(1+obj.taper)))*obj.Sw*np.sqrt((obj.Sw/obj.A))

    if obj.braced_wing:
        obj.CD0  /= obj.Drag_increase
        kq_strut = 1
        l_strut  = 1.1/np.cos((45*np.pi/180))
        strut_taper = 1
        c_strut  = 0.1
        t_strut  = 0.04
        tc_strut = t_strut/c_strut
        S_strut  = c_strut * l_strut
        A_strut  = l_strut**2/S_strut
        interference_penalty = 2
        Qw_strut = 2*((kq_strut*tc_strut)/(np.sqrt(1+strut_taper)))*S_strut*np.sqrt(S_strut/A_strut)*interference_penalty
        obj.Drag_increase = 1 + Qw_strut/obj.Qw_wing
        obj.CD0  *= obj.Drag_increase                                                                         #kg, estimate using length and density of AL2024 t3
    else:
        obj.Drag_increase = 1
    
    