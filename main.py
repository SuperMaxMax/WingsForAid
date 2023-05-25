from aerodynamics.example import print_name, overwrite_value
from parameters import UAV

#     # --- saving
#     # save all attributes of object to csv file
#     members = [attr for attr in dir(concept) if not callable(getattr(concept, attr)) and not attr.startswith("__")]
#     values = [getattr(concept, member) for member in members]

#     # remove brackets and round values
#     values = [value[0] if isinstance(value, np.ndarray) else value for value in values]
#     values = [round(value, 4) if isinstance(value, float) else value for value in values]

#     # add to dataframe
#     df[concept.name] = values

# # set index of dataframe
# df.index = members

# # export dataframe of current design to csv file
# df['DET_CON_2_braced'].to_csv('DET_CON_2_braced.csv', sep=';')

# # remove row in dataframe if all values in that row are the same
# if remove_duplicates == True:
#     for i in df.index:
#         if all(element == df.loc[i].values[0] for element in df.loc[i].values):
#             df.drop(i, inplace=True)
        
# # save dataframe to csv file
# df.to_csv('concept_comparison.csv', sep=';')