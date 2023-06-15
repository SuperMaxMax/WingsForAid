import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np

from parameters import UAV
aircraft = UAV('aircraft')

### Defining Material - GFRP ###
yield_strength = 251E6 * 0.33                                       # yield strength from material book [Pa]
elastic_modulus = 9E9                                               # E from material book [Pa]
strain = yield_strength/elastic_modulus                             # Strain 
T_GFRP = 0.5 * yield_strength * strain                              # toughness of GRP [J/m3]
rho_GFRP = 1500                                                     # density of GRP [kg/m3] from https://www.idemitsu.com/en/business/ipc/engineering/grade/glass/grade.html
                
### Defining Material - Steel 4130 ###              
ultimate_strength = 560E6                                           # ultimate strength from material book [Pa]
yield_strength = 460E6                                              # yield strength from material book [Pa]
elastic_modulus = 205E9                                             # E from material book [Pa]
strain = yield_strength/elastic_modulus                             # Strain [-]
strain_at_break = 0.215                                             # Strain at failure strength [-] 
T_steel = 0.5 * yield_strength * strain + (strain_at_break-strain) *  yield_strength + 0.5 * (ultimate_strength-yield_strength) * (strain_at_break-strain)                             # toughness of GRP [J/m3]
rho_steel = 7850                                                    # density of GRP [kg/m3] from https://www.idemitsu.com/en/business/ipc/engineering/grade/glass/grade.html

### Defining Hail Stone ### 
#info from https://www.nssl.noaa.gov/education/svrwx101/hail/
d_hail = 0.0254                                     # diameter hail stone [m]
#d_hail = 0.04445
r_hail = d_hail/2                                   # radius hail stone [m]
v_hail = 11.176                                     # velocity hail stone [m/s]
#v_hail = 17.8816
rho_hail = 200                                      # denisty hail stone [kg/m3]
A_hail = np.pi * (r_hail**2)                        # impact area [m2]

### Defining Chicken ### 
d_bird = 0.2                                        # diameter chicken [m]
r_bird = d_bird/2                                   # radius chicken [m]
v_bird = aircraft.V_cruise + 25/3.6                 # velocity cruise [m/s]
M_bird = 1.0                                        # mass chicken [kg]
A_bird = np.pi * (r_bird**2)                        # impact area [m2]


### Kinematic Energy Hail Stone Calculation ###
V_hail = (4/3) * np.pi * (r_hail**3)                # volume hail stone [m3]
M_hail = V_hail * rho_hail                          # mass hail stone [kg]
E_kin_hail = 0.5 * M_hail * (v_hail**2)             # kinematic energy hail [J]

### Kinematic Energy bird Calculation ###
E_kin_bird = 0.5 * M_bird * (v_bird**2)             # kinematic energy hail [J]

### Thickness sizing ###
t_fairing = (E_kin_hail)/(A_hail*T_GFRP)            # required thickness fairing

### Thickness sizing ###
t_nose = (E_kin_bird)/(A_bird*T_steel)              # required thickness nose

### Main cabin ###
circumference = 2.4                                                             # circumfernece of fuselage from CAD
surface_main_fuselage_cone = aircraft.l_fus_main_cone * circumference           # wetted area main cone of fuselage [m2]
m_main_fuselage_cone = surface_main_fuselage_cone * t_fairing * rho_GFRP        # mass fairing of main cone of fuselage [kg]
surface_main_fuselage_cone_bottom = aircraft.l_fus_main_cone * 1.0 - 12 * 0.4 * 0.4
m_main_fuselage_cone_bottom = surface_main_fuselage_cone_bottom * t_fairing * rho_GFRP * 2

### Tail cone ###
surface_tail_fuselage_cone = 1.55385                                            # wetted area tail cone of fuselage [m2] 
m_tail_fuselage_cone = surface_tail_fuselage_cone * t_fairing * rho_GFRP        # mass fairing of tail cone [kg]

### Nose cone ###
surface_nose_fuselage_cone = 2.72                                               # wetted area tail cone of fuselage [m2] 
m_nose_fuselage_cone = surface_nose_fuselage_cone * t_nose * rho_steel          # mass fairing of tail cone [kg]

### Mass fairing ###
mass_fairing = m_main_fuselage_cone + m_tail_fuselage_cone + m_nose_fuselage_cone + m_main_fuselage_cone_bottom
