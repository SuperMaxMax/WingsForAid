import sys
sys.path.append("..")

# Start your import below this
# import packages
import numpy as np
import scipy as sp
# Import other files
from parameters import UAV

aircraft = UAV('aircraft')

# print all attributes of aircraft
# print(aircraft.__dict__)

def print_name(obj):
    # print aircraft name
    print(obj.name)

# def overwrite_value(obj):
#     # overwrite aircraft name
#     obj.name = 'new name'

print_name(aircraft)
# overwrite_value(aircraft)
print_name(aircraft)

# ---------------------------------------------------------------------------------
"""Take-off Simulation"""


