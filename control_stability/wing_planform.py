import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad    
import math
import sys

sys.path.append('..')

from parameters import UAV
ac = UAV("aircraft")

# Adjust font size
# plt.rcParams.update({'font.size': 22})

# --------------------------------------------- Parameters ---------------------------------------------
def WingPlanform(ac, plot):
    line_width = 3
    # Parameters
    AR = ac.A                                    # [-]
    S = ac.Sw                                    # [m^2]
    c4_sweep = ac.sweep_co4                     # [rad]
    span = ac.b                                  # [m]
    taper_ratio = ac.taper                       # [-]
    c_r = ac.rootchord                           # [m]
    c_t = ac.tipchord                            # [m]

    # For aileron
    y1_a = ac.ystart_ail             # Start of aileron - [m]
    y2_a = ac.yend_ail               # End of aileron - [m]
    c_a = 1 - ac.xc_aft_spar                        # aileron chord - [c]

    # For flaps
    y1_f = 0                         # Start of flaps - [m]
    y2_f = ac.yend_flap              # End of flaps - [m]
    c_f = 1 - ac.xc_aft_spar         # flap chord - [c]

    # Spars
    f_spar = 0.15                    # front spar location - [c]
    a_spar = ac.xc_aft_spar          # aft spar location - [c]

    # -------------------------------------------- Calculations --------------------------------------------

    # Calculation MAC: chord and location
    # Integrate from root to tip 
    ymin = 0
    ymax = span/2

    # Define function to integrate
    def function(y): 
        return (c_r - ((c_r-c_t)/(span/2))*y)**2

    # Integrate 
    result, error = quad(function, ymin, ymax)

    # Apply formulas to find MAC and spanwise position
    mac = (2/S) * result
    spanwise_pos = (mac - c_r) / -((c_r-c_t)/(span/2))

    # Root chord position
    q_c_r = 0
    le_c_r = q_c_r + 0.25 * c_r
    te_c_r = q_c_r - 0.75 * c_r

    # Tip chord position
    q_c_t = q_c_r - (span/2) * math.tan(c4_sweep)
    le_c_t = q_c_t + 0.25 * c_t
    te_c_t = q_c_t - 0.75 * c_t

    # MAC position
    q_c_mac = q_c_r - (spanwise_pos) * math.tan(c4_sweep)
    le_mac = q_c_mac + 0.25 * mac
    te_mac = q_c_mac - 0.75 * mac

    # Calculating XLEMAC
    XLEMAC = spanwise_pos * math.tan(c4_sweep) + (1/4)*c_r - (1/4)*mac

    # Calculate leading and trailing edge sweep angle
    slope_le = (le_c_r - le_c_t) / (span/2)
    le_sweep = math.atan(slope_le)

    slope_te = (te_c_r - te_c_t) / (span/2)
    te_sweep = math.atan(slope_te)

    # Calculate c/2 sweep angle
    slope_c2 = ((le_c_r-0.5*c_r) - (le_c_t-0.5*c_t)) / (span/2)
    c2_sweep = math.atan(slope_c2)

    # Calculate LE spar sweep angle
    slope_le_sp = ((le_c_r-0.2*c_r) - (le_c_t-0.2*c_t)) / (span/2)
    le_sp_sweep = math.atan(slope_le_sp)

    # Calculate TE spar sweep angle
    slope_te_sp = ((le_c_r-0.75*c_r) - (le_c_t-0.75*c_t)) / (span/2)
    te_sp_sweep = math.atan(slope_te_sp)

    # TE aileron positions
    slope_ai = ((le_c_r-(1-c_a)*c_r) - (le_c_t-(1-c_a)*c_t)) / (span/2)
    ai1 = -(0.75-c_a)*c_r-y1_a*slope_ai
    ai2 = -0.75*c_r-y1_a*slope_te
    ai3 = -(0.75-c_a)*c_r-y2_a*slope_ai
    ai4 = -0.75*c_r-y2_a*slope_te

    # TE flaps positions
    slope_fl = ((le_c_r-(1-c_f)*c_r) - (le_c_t-(1-c_f)*c_t)) / (span/2)
    fl1 = -(0.75-c_f)*c_r-y1_f*slope_fl
    fl2 = -0.75*c_r-y1_f*slope_te
    fl3 = -(0.75-c_f)*c_r-y2_f*slope_fl
    fl4 = -0.75*c_r-y2_f*slope_te

    # ------------------------------------------ Printing results ------------------------------------------

    # print("-------------------- Summary ----------------------")
    # print("Span:                        ", round(span,2), "[m]")
    # print("Taper ratio:                 ", round(taper_ratio,3), "[-]")
    # print("Root chord:                  ", round(c_r,2), "[m]")
    # print("Tip chord:                   ", round(c_t,2), "[m]")
    # print("MAC:                         ", round(mac,2), "[m]")
    # print("Spanwise position MAC:       ", round(spanwise_pos,2), "[m]")
    # print("XLEMAC:                      ", round(XLEMAC,2), "[m]")
    # print("Leading edge sweep angle:    ", round(le_sweep*(180/math.pi),2), "[deg]")
    # print("Front spar sweep angle:      ", round(le_sp_sweep*(180/math.pi),2), "[deg]")
    # print("Half chord sweep angle:      ", round(c2_sweep*(180/math.pi),2), "[deg]")
    # print("Aft spar sweep angle:        ", round(te_sp_sweep*(180/math.pi),2), "[deg]")
    # print("Trailing edge sweep angle:   ", round(te_sweep*(180/math.pi),2), "[deg]")
    # print("---------------------------------------------------")

    # ---------------------------------------------- Plotting ----------------------------------------------

    # Make plot
    fig1, ax1 = plt.subplots(figsize=(10,6))
    # plt.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    
    ax1.set_xlim(-span/2 - 2, span/2 + 2)
    ax1.set_ylim(te_c_r - 0.5, le_c_r + 0.5)

    # Fuselage
    ax1.plot((-ac.w_out/2, -ac.w_out/2), (-3, 3), 'black', label = "Fuselage", linewidth = line_width)
    ax1.plot((ac.w_out/2, ac.w_out/2), (-3, 3), 'black', linewidth = line_width)

    # Quarter chord line left side
    ax1.plot((0 - ac.w_out/2, -span/2 - ac.w_out/2),(q_c_r,q_c_t), label="Quarter Chord line", linewidth=line_width)

    # Spars
    ax1.plot((0 + ac.w_out/2, span/2 + ac.w_out/2),(le_c_r-f_spar*c_r,le_c_t-f_spar*c_t), 'm', label="Spars", linewidth=line_width)
    ax1.plot((0 - ac.w_out/2, -span/2 - ac.w_out/2),(le_c_r-f_spar*c_r,le_c_t-f_spar*c_t), 'm', linewidth=line_width)
    ax1.plot((0 + ac.w_out/2, span/2 + ac.w_out/2),(le_c_r-a_spar*c_r,le_c_t-a_spar*c_t), 'm', linewidth=line_width)
    ax1.plot((0 - ac.w_out/2, -span/2 - ac.w_out/2),(le_c_r-a_spar*c_r,le_c_t-a_spar*c_t), 'm', linewidth=line_width)

    # Aileron
    ax1.plot((y1_a + ac.w_out/2, y2_a + ac.w_out/2),(ai1,ai3),'y', label="Ailerons", linewidth=line_width)
    ax1.plot((y1_a + ac.w_out/2, y1_a + ac.w_out/2),(ai1,ai2),'y', linewidth=line_width)
    ax1.plot((y2_a + ac.w_out/2, y2_a + ac.w_out/2),(ai3,ai4),'y', linewidth=line_width)
    ax1.plot((-y1_a - ac.w_out/2, -y2_a - ac.w_out/2),(ai1,ai3),'y', linewidth=line_width)
    ax1.plot((-y1_a - ac.w_out/2, -y1_a - ac.w_out/2),(ai1,ai2),'y', linewidth=line_width)
    ax1.plot((-y2_a - ac.w_out/2, -y2_a - ac.w_out/2),(ai3,ai4),'y', linewidth=line_width) 

    # Flaps
    ax1.plot((y1_f + ac.w_out/2, y2_f + ac.w_out/2),(fl1,fl3),'c', label="Flaps", linewidth=line_width)
    ax1.plot((y1_f + ac.w_out/2, y1_f + ac.w_out/2),(fl1,fl2),'c', linewidth=line_width)
    ax1.plot((y2_f + ac.w_out/2, y2_f + ac.w_out/2),(fl3,fl4),'c', linewidth=line_width)
    ax1.plot((-y1_f - ac.w_out/2, -y2_f - ac.w_out/2),(fl1,fl3),'c', linewidth=line_width)
    ax1.plot((-y1_f - ac.w_out/2, -y1_f - ac.w_out/2),(fl1,fl2),'c', linewidth=line_width)
    ax1.plot((-y2_f - ac.w_out/2, -y2_f - ac.w_out/2),(fl3,fl4),'c', linewidth=line_width) 

    # Root chord
    ax1.plot((0 + ac.w_out/2, 0 + ac.w_out/2),(le_c_r,te_c_r), 'limegreen', label="Root chord", linewidth=line_width)
    ax1.plot((0 - ac.w_out/2, 0 - ac.w_out/2),(le_c_r,te_c_r), 'limegreen', linewidth=line_width)

    # Tip chords
    ax1.plot((span/2 + ac.w_out/2, span/2 + ac.w_out/2),(le_c_t,te_c_t), 'r', linewidth=line_width)
    ax1.plot((-span/2 - ac.w_out/2, -span/2 - ac.w_out/2),(le_c_t,te_c_t), 'r', label="Tip chord", linewidth=line_width)

    # MAC on left side
    mac_plot = ax1.plot((-spanwise_pos - ac.w_out/2, -spanwise_pos - ac.w_out/2),(te_mac,le_mac), label="MAC", linewidth=line_width)

    # Connecting
    ax1.plot((-span/2 - ac.w_out/2, 0 - ac.w_out/2),(le_c_t,le_c_r), 'b', linewidth=line_width)
    ax1.plot((0 + ac.w_out/2, span/2 + ac.w_out/2),(le_c_r,le_c_t), 'b', linewidth=line_width)
    ax1.plot((-span/2 - ac.w_out/2, 0 - ac.w_out/2),(te_c_t,te_c_r), 'b', linewidth=line_width)
    ax1.plot((0 + ac.w_out/2, span/2 + ac.w_out/2),(te_c_r,te_c_t), 'b', linewidth=line_width)

    # Finalizing
    fig1.gca().set_aspect('equal', adjustable='box')
    # ax1.set_title("Wing planform")
    ax1.set_xlabel("Spanwise position [m]")
    ax1.set_ylabel("Longitudinal position [m]")
    ax1.legend(loc = "upper center", fontsize = "medium", bbox_to_anchor = (0.5, 1.4), ncol = 4)
    ax1.grid()

    plt.savefig("WingPlanform.png", bbox_inches = "tight")

    if plot:
        plt.show(block=True)