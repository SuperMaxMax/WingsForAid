
from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import geometry_determination as geo

concept_1 = UAV('concept_1')
concept_2 = UAV('concept_2')
concept_3 = UAV('concept_3')
concept_4 = UAV('concept_4')
concept_5 = UAV('concept_5')

# iteration
it = True

for concept in [concept_1]:#[concept_1, concept_2, concept_3, concept_4, concept_5]:
    while it:
        # class 1
        c1.iteration(concept)
        W_TO_c1 = concept.W_TO
        print(W_TO_c1)
        
        # geometry determination
        geo.geometry_determination(concept)

        # class 2
        c2.weight_empty(concept)
        W_TO_c2 = concept.W_TO
        print(W_TO_c1)

        # iterate between class 1 and class 2
        # change = (W_TO_c2 - W_TO_c1)/W_TO_c1
        # if abs(change) < 0.01:
        #     it = False
        # else:
        #     continue

        # cg calculation
        c2.cg_calc(concept)
        
        it = False
    
    print(f"Concept {concept.name} has a weight of {concept.W_TO} kg")