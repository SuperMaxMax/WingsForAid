import pandas as pd
import copy
import numpy as np

# import excel file with material properties
material_excel = pd.read_excel('material_properties.xlsx', sheet_name='Sheet2', header=0, index_col=0)

def rank_material(material_excel, weights):
    material_ranking = copy.deepcopy(material_excel)
    print(material_ranking.columns)
    

    ascending_tf = [True, True, True, True, False, False, False]
    for i, property in enumerate(material_ranking.columns):
        property_ranking = material_ranking[property].rank(ascending=ascending_tf[i])
        material_ranking[property] = property_ranking

    material_ranking['Ranking'] = material_ranking.apply(lambda row: sum(row * weights), axis=1)
    material_ranking.sort_values(by=['Ranking'], inplace=True)
    material_ranking.drop(material_ranking.columns[6:], axis=1, inplace=True)

    print(material_ranking)


test = [1, 0, 0, 0, 0, 0, 0]

rank_material(material_excel, test)
# # example of how to print a property
# print(material['6061 T6']['Density'])

# # all materials in dictionary
# print(material.keys())

# # all properties of a material
# print(material['6061 T6'].keys())



