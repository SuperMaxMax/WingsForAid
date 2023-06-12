print("""Wings For Aid

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
     0.5089, """    

# x,y,z bias for whole surface
TRANSLATE
    """,2.276,     0.00000,     0.00000, """

# Here the section starts
#---------------Root Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,0.0,          0.0,         0.0,         1.490,       0.000, """ 

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
flap     """,1.0,    0.75,    0.0, 1.0, 0.0,    1.0,"""

#---------------End Flap Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,0.04278,      1.5629,      0.0,         1.3189,      0.000,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
flap     """,1.0,    0.75,    0.0, 1.0, 0.0,    1.0,"""


#---------------Begin Aileron Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,0.09054,      3.308,       0.0,         1.1278,      0.000,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
aileron  """,1.0,    0.75,    0.0, 1.0, 0.0,   -1.0,"""

#---------------End Aileron Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
    """,0.10879,      3.975,       0.0,         1.0407,      0.000,"""   

NACA
4415

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
aileron  """,1.0,    0.75,    0.0, 1.0, 0.0,  -1.0,"""

#-----------Tip Section------------------
SECTION
    """,0.130375,    4.7635,    0.0,        0.9685,        0.000,"""   

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
     """,0.0,"""  

# x,y,z bias for whole surface
TRANSLATE
   """,6.56,     0.0,     0.0,"""

#-----------------------Root Section-----------------------------------
#    Xle         Yle         Zle         chord       angle 
SECTION
   """,0.0,          0.0,        0.0,         0.931167,      0.000,"""   

NACA
0012

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
elevator  """,1.0,    0.0,    0.0, 1.0, 0.0,    1.0,"""

#-------------------Tip Section-------------------------
#    Xle         Yle         Zle         chord       angle 
SECTION
     """,0.0,         1.0942,       0.0,       0.931167,     0.000,"""  

NACA
0012

CONTROL
#Cname   Cgain  Xhinge  HingeVec       SgnDup
elevator  """,1.0,    0.0,    0.0, 1.0, 0.0,    1.0,"""

#
#==============================================================
#
SURFACE
V-Stab
""",8,  1.0,  10,  0.75,"""  ! Nchord   Cspace

# twist angle bias for whole surface
ANGLE
     """,0.0,"""  

# x,y,z bias for whole surface
TRANSLATE
   """,6.438067,     0.00000,     0.00000,"""
#-------------------Lower Section-------------------------------------------
#    Xle         Yle         Zle         chord       angle   
SECTION
   """,0.0,           0.0,         0.0,        0.9530,       0.000,"""   

NACA 
0012

#CONTROL
#rudder    """,1.0,    0.40,   0.0, 0.0, 1.0,"""
#-----------------Upper Section----------------------
SECTION
   """,0.8402,          0.0,        1.2,        0.6672,      0.000,"""  

NACA
0012

#CONTROL
#rudder    """,1.0,    0.40,    0.0, 0.0, 1.0,"""
#==============================================================

#Created by Jan Vonken 10-06-2023"""
)