import sys
sys.path.append("..")

# Start your import below this
from parameters import UAV
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad



I_min =
a = np.linspace(0,0.1,100) #semi_minor of ellipse
b = 2*a #semi_major of ellipse
t = 4*I_min/(np.pi*a**3 * (1+(3*b/a))) #thickness of

