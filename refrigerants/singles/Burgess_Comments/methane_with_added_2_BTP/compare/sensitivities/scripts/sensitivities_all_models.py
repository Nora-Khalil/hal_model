import cantera as ct
from PIL import Image
from subprocess import run
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import re 
import os
import numpy as np
import csv


###################### get the model #######################################

model_number = sys.argv[1]

directory = f'{model_number}'

file_name = directory.split('/')




############### if wanted to loop over all of the models from (120-140, plus 145)

# directory = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/cantera/Nora/2_BTP/FFCM_seed/2_BTP_seed/chemkin/copies'

# models = [file for file in os.listdir(directory) if re.search('copy_chem0([0-9]+).cti', file)]


############### some values that are constant and will be used in the loop 

BTPmole_list = list(np.linspace(0.0, 0.16268, 10))
species_con =  {'CH4(3)': 0.5, 'O2(4)': 1, 'N2': 3.76, '2-BTP(1)': BTPmole_list[3]}



############### if using a loop, start the loop here

sensitivities_across_models = dict() # this dictionary will store all the data across all the models 

#for file in models: 

match = re.search('copy_chem0([0-9]+).cti', file_name[-1])
assert match.group(1)

########### make the gas, specify the TPX, create the FreeFlame object in Cantera 

gas = ct.Solution(directory)
gas.TPX = 298, ct.one_atm, species_con

########### solve the flame

f = ct.FreeFlame(gas)
f.solve()
    
########## save a flux diagram 
diagram = ct.ReactionPathDiagram(f.gas,'H')

diagram.title = 'Reaction path diagram following H'
diagram.label_threshold = 0.01

dot_file = f'rxnpath{match.group(1)}.dot'
img_file = f'rxnpath{match.group(1)}.png'
img_path = Path.cwd().joinpath(img_file)

diagram.write_dot(dot_file)
#print(diagram.get_data())

print("Wrote graphviz input file to '{0}'.".format(Path.cwd().joinpath(dot_file)))

run('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file).split())
print("Wrote graphviz output file to '{0}'.".format(img_path))

########## get sensitivities of each reaction rate on flamespeed
'''
.get_flame_speed_reaction_sensitivities(): Compute the normalized sensitivities of the laminar flame speed S with respect to the reaction rate constants k:

                    s_i =   k_i    dS_u
                           -----  -------
                            S_U    dk_i

'''
######## sort the sensitivities by magnitude #######################################
sens = f.get_flame_speed_reaction_sensitivities()

sensitivity = {}
for m in range(gas.n_reactions):
    sensitivity[m] = abs(sens[m])

sorted_sensitivity = dict(sorted(sensitivity.items(), key=lambda item: item[1], reverse=True))

######### revert the sensitivity values back to original sign  ####################

sorted_sensitivity_list = [[k,sens[k],gas.reaction(k).equation] for k,v in sorted_sensitivity.items() ]

######## only store the top 10 most sensitive reactions #######################

rxn_sensitivity_values = [triplet for triplet in sorted_sensitivity_list[0:20]]

print(rxn_sensitivity_values)

 #sensitivities_across_models[int(match.group(1))] = rxn_sensitivity_values
    
    
# ############# sort sensitivities for each model ############################

# sorted_sensitivities_across_models = dict(sorted(sensitivities_across_models.items(), key=lambda item: item[0], reverse=False))



# ############# plot the sensitivities for each model  ############################

# for k,v in sorted_sensitivities_across_models.items():

#     fig, ax = plt.subplots()
#     for row in v: 
#         ax.barh(row[2], row[1], color='b', label=row[0], align='center')
#         plt.xlabel('Sensitivity')
#         plt.title(f"Sensitivity for Model {k}")
    



with open(f'sensitivites_{file_name[-1]}.csv', 'w+') as g:
    g.write(directory)
    g.write('\n')
    writers = csv.writer(g)
    writers.writerow(rxn_sensitivity_values)


