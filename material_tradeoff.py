import numpy as np 


materials = np.array((),
                     (),
                     (),
                     ()
                     )

weights_fuel_tank = [0, 0.5 0.7, ]

score = 0
for i in len(materials):
    score_new = materials[i]*weights_fuel_tank
    if score_new > score: 
        score = score_new
        material_choice = materials[i]

print("yey you have selected", material_choice)