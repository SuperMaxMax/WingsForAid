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

### Defining Hail Stone ### 
#info from https://www.nssl.noaa.gov/education/svrwx101/hail/
d_hail = 0.0254                                     # diameter hail stone [m]
#d_hail = 0.04445
r_hail = d_hail/2                                   # radius hail stone [m]
v_hail = 11.176                                     # velocity hail stone [m/s]
#v_hail = 17.8816
rho_hail = 200                                      # denisty hail stone [kg/m3]
A_hail = np.pi * (r_hail**2)                        # impact area [m2]

### Kinematic Energy Hail Stone Calculation ###
V_hail = (4/3) * np.pi * (r_hail**3)                # volume hail stone [m3]
M_hail = V_hail * rho_hail                          # mass hail stone [kg]
E_kin_hail = 0.5 * M_hail * (v_hail**2)             # kinematic energy hail [J]

### Thickness sizing ###
t_sheet = (E_kin_hail)/(A_hail*T_GFRP)            # required thickness fairing

t_nose_tail = 0.002

a = 0.005
t_stringer = 0.001

b = 0.96666
h = 0.67
w = 1

F = 10 * 9.80665

def stress():
    stress_panel = 100000000000   
    n = 0
    while stress_panel > yield_strength:
        y = ((t_sheet*0.5 + a)*8*a*t_stringer*n)/(8*a*t_stringer*n+t_sheet*b)
        #Ixx = ((1/12) * t_stringer * ((2*a)**3) + 4 * t_stringer * (a**3)) * n + (1/12) * b * (t_sheet**3) + b * t_sheet * (a**2)
        Ixx = ((1/6) * t_stringer * (2*a)**3 + 4 * a * t_stringer * (y-a-0.5*t_sheet)**2 + 2*a*t_stringer*(2*a+0.5*t_stringer+0.5*t_sheet-y)**2 + 2*a*t_stringer*(0.5*t_sheet+0.5*t_stringer-y)**2) * n + (1/12)*b*t_sheet**3 + b*t_sheet*(0.5*t_sheet-y)**2
        r = max(0.5*t_sheet+y, t_sheet*0.5+2*t_stringer+2*a-y)
        stress_panel = F * h * r / Ixx
        mass = (a*8*t_stringer*n*h + b*h*t_sheet + a*8*t_stringer*2*b)*rho_GFRP
        n += 1 
    return mass

main_sides_mass = 6*stress()

bottom_mass = (t_sheet * b*3 * w + 4*a*8*t_stringer*b*3 + 24*a*8*t_stringer*0.45 - t_sheet*0.45*0.45*12) * rho_GFRP

top_mass = (3*b*w*t_sheet) * rho_GFRP

tail_mass = 1.55385 * t_nose_tail * rho_GFRP

nose_mass = 2.72 * t_nose_tail * rho_GFRP

mass_total = main_sides_mass + bottom_mass + top_mass + tail_mass + nose_mass
print(mass_total)
print(main_sides_mass,bottom_mass,top_mass,nose_mass,tail_mass)