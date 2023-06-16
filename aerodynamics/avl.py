import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize


def export(aircraft):
    df = pd.read_csv("AVL_mass.csv")
    masses = np.savetxt("AVL_mass.txt", df, delimiter=" ", fmt='%s')

    f = open('AVL_mass.txt', 'r')
    file_contents = f.read()

    avl_mass = """#-------------------------------------------------
    #  Wings For Aid Aircraft
    #
    #  Dimensional unit and parameter data.
    #  Mass & Inertia breakdown.
    #-------------------------------------------------

    #  Names and scalings for units to be used for trim and eigenmode calculations.
    #  The Lunit and Munit values scale the mass, xyz, and inertia table data below.
    #  Lunit value will also scale all lengths and areas in the AVL input file.
    Lunit = 1.0   m
    Munit = 1.0   kg
    Tunit = 1.0   s

    #------------------------- 
    #  Gravity and density to be used as default values in trim setup (saves runtime typing).
    #  Must be in the unit names given above (m,kg,s).
    g   = 9.81
    rho = 0.904637

    #-------------------------
    #  Mass & Inertia breakdown.
    #  x y z  is location of item's own CG.
    #  Ixx... are item's inertias about item's own CG.
    #
    #  x,y,z system here must be exactly the same one used in the .avl input file
    #     (same orientation, same origin location, same length units)
    #
    #  mass      x      y        z       Ixx   Iyy   Izz    Ixy  Ixz  Iyz
    # 
    """,file_contents,"""

    # Created by Jan Vonken, 11-06-2023"""


    avl_file = """Wings For Aid

    0.1842                   !   Mach
    0     0     0.0          !   iYsym  iZsym  Zsym
    11.7113 1.248  9.527     !   Sref   Cref   Bref   reference area, chord, span
    0.00  0.0   0.0          !   Xref   Yref   Zref   moment reference location (arb.)
    0.02595                  !   CDp

    #
    #==============================================================
    #
    SURFACE
    wing
    #Horseshoe Vortex Distribution
    """,10, ' ', 1.0, ' ' ,   20, ' ' ,   1.0, ' ' , """   ! Nchord   Cspace   Nspan  Sspace
    #
    # reflect image wing about y=0 plane
    YDUPLICATE
        """,0.0,""" 
    #
    # twist angle bias for whole surface
    ANGLE
        """,aircraft.i_w*180/np.pi, """    

    # x,y,z bias for whole surface
    TRANSLATE
        """,aircraft.X_LEMAC, ' ' ,      0.00000, ' ' ,      aircraft.ST_z_ground+aircraft.ST_h_fus+0.075*aircraft.rootchord, """

    # Here the section starts
    #---------------Root Section-------------------------------------------
    #    Xle         Yle         Zle         chord       angle   
    SECTION
        """,0.0, ' ' ,           0.0, ' ' ,          0.0, ' ' ,          aircraft.rootchord, ' ' ,        0.000*180/np.pi, """ 

    NACA
    4415

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    flap     """,1.0, ' ' ,     aircraft.xc_aft_spar, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,     1.0,"""

    #---------------End Flap Section-------------------------------------------
    #    Xle         Yle         Zle         chord       angle   
    SECTION
        """,np.tan(aircraft.sweep_LE)*aircraft.yend_flap, ' ' ,       aircraft.yend_flap, ' ' ,       0.0, ' ' ,         2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.yend_flap)), ' ' ,     -aircraft.yend_flap/(aircraft.b/2) * aircraft.wing_twist*180/np.pi,"""   

    NACA
    4415

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    flap     """,1.0, ' ' ,     aircraft.xc_aft_spar, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,     1.0,"""


    #---------------Begin Aileron Section-------------------------------------------
    #    Xle         Yle         Zle         chord       angle   
    SECTION
        """,np.tan(aircraft.sweep_LE)*aircraft.ystart_ail, ' ' ,       aircraft.ystart_ail, ' ' ,        0.0, ' ' ,          2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.ystart_ail)), ' ' ,       -aircraft.ystart_ail/(aircraft.b/2) * aircraft.wing_twist*180/np.pi,"""   

    NACA
    4415

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    aileron  """,1.0, ' ' ,     aircraft.xc_aft_spar, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,    -1.0,"""

    #---------------End Aileron Section-------------------------------------------
    #    Xle         Yle         Zle         chord       angle   
    SECTION
        """,np.tan(aircraft.sweep_LE)*aircraft.yend_ail, ' ' ,       aircraft.yend_ail, ' ' ,        0.0, ' ' ,          2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.yend_ail)), ' ' ,       -aircraft.yend_ail/(aircraft.b/2) * aircraft.wing_twist*180/np.pi,"""   

    NACA
    4415

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    aileron  """,1.0, ' ' ,     aircraft.xc_aft_spar, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,   -1.0,"""

    #-----------Tip Section------------------
    SECTION
        """,np.tan(aircraft.sweep_LE)*aircraft.b/2, ' ' ,     aircraft.b/2, ' ' ,     0.0, ' ' ,         aircraft.tipchord, ' ' ,         -aircraft.wing_twist*180/np.pi,"""   

    NACA
    4415

    #
    #==============================================================
    #
    SURFACE
    H-Stab
    """,8, ' ' ,   1.0, ' ' ,   5, ' ' ,   1.0,"""  !  Nchord   Cspace
    #
    # reflect image wing about y=0 plane
    YDUPLICATE
        """,0.0,""" 

    # twist angle bias for whole surface
    ANGLE
        """,aircraft.i_w_h,"""  

    # x,y,z bias for whole surface
    TRANSLATE
    """,aircraft.X_LEMAC + 0.25 * aircraft.rootchord + aircraft.l_h - 0.25*aircraft.AE_rootchord_h, ' ' ,      0.0, ' ' ,      aircraft.ST_z_ground+aircraft.ST_h_fus+0.06*aircraft.AE_rootchord_h ,"""

    #-----------------------Root Section-----------------------------------
    #    Xle         Yle         Zle         chord       angle 
    SECTION
    """,0.0, ' ' ,           0.0, ' ' ,         0.0, ' ' ,          aircraft.AE_rootchord_h, ' ' ,       0.000*180/np.pi,"""   

    NACA
    0012

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    elevator  """,1.0, ' ' ,     0.001, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,     1.0,"""

    #-------------------Tip Section-------------------------
    #    Xle         Yle         Zle         chord       angle 
    SECTION
        """,0.0, ' ' ,          aircraft.b_h/2, ' ' ,        0.0, ' ' ,        aircraft.AE_tipchord_h, ' ' ,      0.000*180/np.pi,"""  

    NACA
    0012

    CONTROL
    #Cname   Cgain  Xhinge  HingeVec       SgnDup
    elevator  """,1.0, ' ' ,     0.001, ' ' ,     0.0, ' ' ,  1.0, ' ' ,  0.0, ' ' ,     1.0,"""

    #
    #==============================================================
    #
    SURFACE
    V-Stab
    """,8, ' ' ,   1.0, ' ' ,   10, ' ' ,   aircraft.xc_aft_spar,"""  ! Nchord   Cspace

    # twist angle bias for whole surface
    ANGLE
        """,aircraft.AE_i_w_v*180/np.pi,"""  

    # x,y,z bias for whole surface
    TRANSLATE
    """,aircraft.X_LEMAC + 0.25 * aircraft.rootchord + aircraft.l_h - 0.25*aircraft.AE_rootchord_h, ' ' ,      0.00000, ' ' ,      aircraft.ST_z_ground+aircraft.ST_h_fus,"""
    #-------------------Lower Section-------------------------------------------
    #    Xle         Yle         Zle         chord       angle   
    SECTION
    """,0.0, ' ' ,            0.0, ' ' ,          0.0, ' ' ,         aircraft.AE_rootchord_v, ' ' ,        0.000*180/np.pi,"""   

    NACA 
    0012

    #CONTROL
    #rudder    """,1.0, ' ' ,     0.40, ' ' ,    0.0, ' ' ,  0.0, ' ' ,  1.0,"""
    #-----------------Upper Section----------------------
    SECTION
    """,np.tan(aircraft.sweep_LE_v) * aircraft.AE_b_v, ' ' ,           0.0, ' ' ,         aircraft.AE_b_v, ' ' ,         aircraft.AE_tipchord_v, ' ' ,       0.000*180/np.pi,"""  

    NACA
    0012

    #CONTROL
    #rudder    """,1.0, ' ' ,     0.40, ' ' ,     0.0, ' ' ,  0.0, ' ' ,  1.0,"""
    #==============================================================

    #Created by Jan Vonken 10-06-2023"""

    f.close()

    string = str()
    for i in range(len(avl_file)):
        string = string + str(avl_file[i])

    #print(string)

    path = r"/Users/janvonken/Desktop/TU/3rd Year/DSE/simulating/Avl/OwnCraft/WingsForAid2.avl"""
    f = open(path, "w")
    f.write(string)
    f.close()

    string = str()
    for i in range(len(avl_mass)):
        string = string + str(avl_mass[i])

    #print(string)

    path = r"/Users/janvonken/Desktop/TU/3rd Year/DSE/simulating/Avl/OwnCraft/WingsForAid2.mass"""
    f = open(path, "w")
    f.write(string)
    f.close()
