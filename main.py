from parameters import UAV
import Class_I_weight_estimation as c1
import Class_II_weight_estimation as c2
import Class_II_cg_estimation as c2cg
import geometry_determination as geo
import pandas as pd

concept_1 = UAV('concept_1', 'tractor', boom=True)
concept_2 = UAV('concept_2', 'tractor', boom=False)
concept_3 = UAV('concept_3', 'pusher', boom=False)
concept_4 = UAV('concept_4', 'pusher', boom=False)
concept_5 = UAV('concept_5', 'fuselage', boom=False)

# create dataframe with members and values, to save all concepts in
df = pd.DataFrame()

for concept in [concept_1, concept_2, concept_3, concept_4, concept_5]:
    # iteration
    it = True
    W_TO_c2_old = 6000
    while it:
        # class 1
        c1.run(concept)
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
        change = (W_TO_c2 - W_TO_c2_old)/W_TO_c2_old

        if abs(change) < 0.00001:
            it = False
        else:
            W_TO_c2_old = W_TO_c2

    # cg calculation
    c2cg.cg_calc(concept)

    # save all attributes of object to csv file
    members = [attr for attr in dir(concept) if not callable(getattr(concept, attr)) and not attr.startswith("__")]
    values = [getattr(concept, member) for member in members]

    # round values
    values = [round(value, 3) if isinstance(value, float) else value for value in values]

    # add to dataframe
    df[concept.name] = values

# set index of dataframe
df.index = members

# save dataframe to csv file
df.to_csv('output.csv', sep=';')