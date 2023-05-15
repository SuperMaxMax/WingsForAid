
from parameters import UAV
import Class_I_weight_estimation as C1
import geometry_determination as geo

concept1 = UAV()

C1.iteration(concept1)

geo.wing_sizing(concept1, False)

print(concept1.L_D)
