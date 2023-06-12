import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from parameters import UAV
aircraft = UAV('aircraft')

import wing_planform as wp
import horizontal_tail_design as htd
import vertical_tail_design as vtd



#wp.main_wing_planform(aircraft)
##htd.horizontal_tail_planform(aircraft)
#vtd.horizontal_tail_planform(aircraft)



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
""",10,  1.0,  20,  1.0,"""   ! Nchord   Cspace   Nspan  Sspace
#
# reflect image wing about y=0 plane
YDUPLICATE
     """,0.0,""" 
#
# twist angle bias for whole surface
ANGLE""",
     aircraft.i_w, """    

# x,y,z bias for whole surface
TRANSLATE
    """,aircraft.X_LEMAC,     0.00000,     0.00000, """

# Here the section starts
#---------------Root Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,0.0,          0.0,         0.0,         aircraft.rootchord,       0.000, """ 

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
flap     """,1.0,    0.75,    0.0, 1.0, 0.0,    1.0,"""

#---------------End Flap Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,np.tan(aircraft.AE_sweep_LE)*aircraft.yend_flap,      aircraft.yend_flap,      0.0,        2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.yend_flap)),    -aircraft.yend_flap/(aircraft.b/2) * aircraft.wing_twist,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
flap     """,1.0,    0.75,    0.0, 1.0, 0.0,    1.0,"""


#---------------Begin Aileron Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,np.tan(aircraft.AE_sweep_LE)*aircraft.ystart_ail,      aircraft.ystart_ail,       0.0,         2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.ystart_ail)),      -aircraft.ystart_ail/(aircraft.b/2) * aircraft.wing_twist,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
aileron  """,1.0,    0.75,    0.0, 1.0, 0.0,   -1.0,"""

#---------------End Aileron Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,np.tan(aircraft.AE_sweep_LE)*aircraft.yend_ail,      aircraft.yend_ail,       0.0,         2*aircraft.Sw/(1+aircraft.taper)/aircraft.b*(1-(1-aircraft.taper)/aircraft.b*abs(2*aircraft.yend_ail)),      -aircraft.yend_ail/(aircraft.b/2) * aircraft.wing_twist,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
aileron  """,1.0,    0.75,    0.0, 1.0, 0.0,  -1.0,"""

#-----------Tip Section------------------
SECTION
    """,np.tan(aircraft.AE_sweep_LE)*aircraft.b/2,    aircraft.b/2,    0.0,        aircraft.tipchord,        -aircraft.wing_twist,"""   

NACA
4415

#
#==============================================================
#
SURFACE
H-Stab
""",8,  1.0,  5,  1.0,"""  !  Nchord   Cspace
#
# reflect image wing about y=0 plane
YDUPLICATE
     """,0.0,""" 

# twist angle bias for whole surface
ANGLE
     """,aircraft.AE_i_w_h,"""  

# x,y,z bias for whole surface
TRANSLATE
   """,aircraft.X_LEMAC + 0.25 * aircraft.rootchord + aircraft.l_h - 0.25*aircraft.AE_rootchord_h,     0.0,     0.0,"""

#-----------------------Root Section-----------------------------------
#    Xle         Yle         Zle         chord       angle 
SECTION
   """,0.0,          0.0,        0.0,         aircraft.AE_rootchord_h,      0.000,"""   

NACA
0012

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
elevator  """,1.0,    0.001,    0.0, 1.0, 0.0,    1.0,"""

#-------------------Tip Section-------------------------
#    Xle         Yle         Zle         chord       angle 
SECTION
     """,0.0,         aircraft.AE_b_h/2,       0.0,       aircraft.AE_tipchord_h,     0.000,"""  

NACA
0012

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
elevator  """,1.0,    0.001,    0.0, 1.0, 0.0,    1.0,"""

#
#==============================================================
#
SURFACE
V-Stab
""",8,  1.0,  10,  0.75,"""  ! Nchord   Cspace

# twist angle bias for whole surface
ANGLE
     """,aircraft.AE_i_w_v,"""  

# x,y,z bias for whole surface
TRANSLATE
   """,aircraft.X_LEMAC + 0.25 * aircraft.rootchord + aircraft.l_h - 0.25*aircraft.AE_rootchord_h,     0.00000,     0.00000,"""
#-------------------Lower Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
   """,0.0,           0.0,         0.0,        aircraft.AE_rootchord_v,       0.000,"""   

NACA 
0012

#CONTROL
#rudder    """,1.0,    0.40,   0.0, 0.0, 1.0,"""
#-----------------Upper Section----------------------
SECTION
   """,np.tan(aircraft.AE_lambda_LE_v) * aircraft.b_v,          0.0,        aircraft.b_v,        aircraft.AE_tipchord_v,      0.000,"""  

NACA
0012

#CONTROL
#rudder    """,1.0,    0.40,    0.0, 0.0, 1.0,"""
#==============================================================

#Created by Jan Vonken 10-06-2023"""


f = open(r"WFA.avl", "w")

f.write(avl_file)

