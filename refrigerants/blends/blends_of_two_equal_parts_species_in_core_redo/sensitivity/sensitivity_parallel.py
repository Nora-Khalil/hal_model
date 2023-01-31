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
import pandas as pd
import sys


###################### get the model #######################################


blend_name = sys.argv[1]

directory = f'/work/westgroup/nora/Code/projects/halogens/refrigerants/blends/blends_of_two_equal_parts_species_in_core_redo/{blend_name}/chemkin/copies/copy_chem.cti'


##################### here's the mapping of the species names in each rmg-built blend model #########################

species_mapping_for_rmg_built_blends = {'C2H5F_CH2F2': ['C2H5F(1)', 'CH2F2(34)', 'O2(2)'], 'C2H5F_CH2FCH2F': ['C2H5F(1)', 'S(179)', 'O2(2)'], 'C2H5F_CH2FCHF2': ['C2H5F(1)', 'S(179)', 'O2(2)'], 'C2H5F_CH3CHF2': ['C2H5F(1)', 'S(179)', 'O2(2)'], 'C2H5F_CH3F': ['C2H5F(1)', 'CH3F(34)', 'O2(2)'], 'CH2F2_CH3CF3': ['CH2F2(1)', 'CH3CF3(87)', 'O2(2)'], 'CH2F2_CH3F': ['CH2F2(1)', 'CH3F(22)', 'O2(2)'], 'CH2FCH2F_CH2F2': ['C2H4F2(1)', 'CH2F2(40)', 'O2(2)'], 'CH3CF3_C2H5F': ['CH3CF3(1)', 'C2H5F(188)', 'O2(2)'], 'CH3CF3_CH2FCH2F': ['CH3CF3(1)', 'S(188)', 'O2(2)'], 'CH3CF3_CH2FCHF2': ['CH3CF3(1)', 'S(188)', 'O2(2)'], 'CH3CF3_CH3CHF2': ['CH3CF3(1)', 'S(188)', 'O2(2)'], 'CH3CHF2_CH2F2': ['CH3CHF2(1)', 'CH2F2(40)', 'O2(2)'], 'CH3CHF2_CH2FCHF2': ['CH3CHF2(1)', 'S(180)', 'O2(2)'], 'CH3F_CH3CF3': ['CH3F(1)', 'CH3CF3(81)', 'O2(2)']}


vol_frac_of_max = {'C2H5F_CH2F2': 0.10258620689655173, 'C2H5F_CH2FCH2F': 0.0793103448275862, 'C2H5F_CH3CHF2': 0.0793103448275862, 'C2H5F_CH3F': 0.09482758620689655, 'CH2F2_CH3CF3': 0.12586206896551724, 'CH2F2_CH3F': 0.14913793103448275, 'CH2FCH2F_CH2F2': 0.1103448275862069, 'CH3CF3_C2H5F': 0.0793103448275862, 'CH3CF3_CH2FCH2F': 0.08706896551724139, 'CH3CF3_CH2FCHF2': 0.09482758620689655, 'CH3CF3_CH3CHF2': 0.08706896551724139, 'CH3CHF2_CH2F2': 0.1103448275862069, 'CH3CHF2_CH2FCHF2': 0.08706896551724139, 'CH3F_CH3CF3': 0.1103448275862069}

# {'C2H5F_CH2F2': 0.10258620689655173, 'C2H5F_CH2FCH2F': 0.0793103448275862, 'C2H5F_CH2FCHF2': 0.08706896551724139, 'C2H5F_CH3CHF2': 0.0793103448275862, 'C2H5F_CH3F': 0.09482758620689655, 'CH2F2_CH3CF3': 0.12586206896551724, 'CH2F2_CH3F': 0.14913793103448275, 'CH2FCH2F_CH2F2': 0.1103448275862069, 'CH2FCH2F_CH2FCHF2': 0.08706896551724139, 'CH2FCH2F_CH3CHF2': 0.0793103448275862, 'CH2FCH2F_CH3F': 0.10258620689655173, 'CH2FCHF2_CH2F2': 0.12586206896551724, 'CH2FCHF2_CH3F': 0.1103448275862069, 'CH3CF3_C2H5F': 0.09482758620689655, 'CH3CF3_CH2FCH2F': 0.08706896551724139, 'CH3CF3_CH2FCHF2': 0.09482758620689655, 'CH3CF3_CH3CHF2': 0.08706896551724139, 'CH3CHF2_CH2F2': 0.1103448275862069, 'CH3CHF2_CH2FCHF2': 0.08706896551724139, 'CH3F_CH3CF3': 0.1103448275862069, 'CH3F_CH3CHF2': 0.10258620689655173}


############### some values that are constant and will be used in the loop 

#get the labels needs
sp1, sp2, oxygen = species_mapping_for_rmg_built_blends[blend_name]

#create the species concentrations
x = vol_frac_of_max[blend_name]
norm_ox = (1-x)*.21
species_con =  {f'{sp1}': (x/2/norm_ox), f'{sp2}': (x/2/norm_ox), f'{oxygen}':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}



########### make the gas, specify the TPX

gas = ct.Solution(directory)
gas.TPX = 298, ct.one_atm, species_con

########### create the FreeFlame object in Cantera, solve the flame

f = ct.FreeFlame(gas)
f.solve()
    
########## save a flux diagram 
diagram = ct.ReactionPathDiagram(f.gas,'H')
diagram.show_details = True
diagram.flow_type = 'NetFlow'
diagram.scale=-1

diagram.title = 'Reaction path diagram following H'
diagram.label_threshold = 0.05

dot_file = f'rxnpath_{blend_name}_rmg_details_2nd_time.dot'
img_file = f'rxnpath_{blend_name}_rmg_details_2nd_time.png'
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
# ######## sort the sensitivities by magnitude #######################################

sens = f.get_flame_speed_reaction_sensitivities()

sensitivity = {}
for m in range(gas.n_reactions):
    sensitivity[m] = abs(sens[m])

sorted_sensitivity = dict(sorted(sensitivity.items(), key=lambda item: item[1], reverse=True)) #sort with highest magnitude first

######### revert the sensitivity values back to original sign  ####################

#sorted_sensitivity_list = [[k,sens[k],gas.reaction(k)] for k,v in sorted_sensitivity.items() ]

data = {
    'k_s': [k for k,v in sorted_sensitivity.items()], #this is number of reaction in gas.reactions list
    'sensitivity': [sens[k] for k,v in sorted_sensitivity.items()], #sensitivity
    'cantera equation': [gas.reaction(k).equation for k,v in sorted_sensitivity.items()],
    'cantera products': [gas.reaction(k).products for k,v in sorted_sensitivity.items()],
    'cantera reactants': [gas.reaction(k).reactants for k,v in sorted_sensitivity.items()],
}

df = pd.DataFrame(data)
df.to_csv(f'{blend_name}_RMG_sensitivities_2nd_time.csv', index=False)




