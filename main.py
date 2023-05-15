
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
        print("After class 1 iteration:")
        print(f"Mff: {concept.Mff}, L/D: {concept.L_D}, W_TO: {concept.W_TO}, W_OE: {concept.W_OE}, W_F: {concept.W_F}")
        W_TO_c1 = concept.W_TO
        
        # geometry determination
        geo.geometry_determination(concept)
        print("After geometry determination:")
        print(f"b: {concept.b} Sw: {concept.Sw} S_G: {concept.S_G}")

        # class 2
        c2.weight_empty(concept)
        print("After class 2 iteration:")
        print(f"W_OE: {concept.W_OE}, W_TO: {concept.W_TO}")
        W_TO_c2 = concept.W_TO

        print("")
        print("----------------------------------")
        print("")

        # iterate between class 1 and class 2
        change = (W_TO_c2 - W_TO_c1)/W_TO_c1
        if abs(change) < 0.01:
            it = False
        else:
            continue

    # cg calculation
    c2.cg_calc(concept)
    
    print(f"Concept {concept.name} has a weight of {concept.W_TO} kg")
