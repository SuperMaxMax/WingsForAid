import pandas as pd

# import excel file with material properties
material_excel = pd.read_excel('material_properties.xlsx', sheet_name='Sheet2', header=0, index_col=0)

# create dictionary of material properties
material = material_excel.to_dict(orient='index')

# # example of how to print a property
# print(material['6061 T6']['Density'])

# # all materials in dictionary
# print(material.keys())

# # all properties of a material
# print(material['6061 T6'].keys())



