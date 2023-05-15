
from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
# import geometry_determination as geo

concept_1 = UAV('concept_1')
concept_2 = UAV('concept_2')
concept_3 = UAV('concept_3')
concept_4 = UAV('concept_4')
concept_5 = UAV('concept_5')

# iteration
it = True

for concept in [concept_1, concept_2, concept_3, concept_4, concept_5]:
    while it:
        c1.iteration(concept)
        W_OE_c1 = concept.W_OE

        # geometry determination stuff

        c2.weight_empty(concept)
        W_OE_c2 = concept.W_OE
    
        change = (W_OE_c2 - W_OE_c1)/W_OE_c1
        if abs(change) < 0.001:
            it = False
        else:
            continue
    
    print(f"Concept {concept.name} has a weight of {concept.W_OE} kg")