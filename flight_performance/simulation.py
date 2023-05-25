import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

# Start your import below this
# import packages
import numpy as np
import scipy as sp
# Import other files
from parameters import UAV

aircraft = UAV("aircraft")

print(aircraft.name)