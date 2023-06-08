import pandas as pd
import copy
import numpy as np

def rank_material(weights, asc):
    material_ranking = copy.deepcopy(material_df)

    ascending_tf = [True, True, True, True, False, False, False] + asc
    for i, property in enumerate(material_ranking.columns):
        property_ranking = material_ranking[property].rank(ascending=ascending_tf[i])
        material_ranking[property] = property_ranking

    material_ranking['Ranking'] = material_ranking.apply(lambda row: sum(row * weights), axis=1)
    material_ranking.sort_values(by=['Ranking'], inplace=True)
    material_ranking.drop(material_ranking.columns.difference(['Ranking']), 1, inplace=True)
    
    return list(material_ranking.index)

# import excel file with material properties
material_df = pd.read_excel('material_properties.xlsx', sheet_name='Sheet2', header=0, index_col=0)

# create dictionary for each material containing all properties
material = {}
for i in range(len(material_df)):
    material[material_df.index[i]] = material_df.iloc[i]



