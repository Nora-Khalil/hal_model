# One_flip for RMG reactions 
# This script switches RMG reactions from one model 
#to the kinetics of another, then checks the flame speed 


####D is Davids (was previously rmg), N is Noras (was previously Nist)

import cantera as ct
import cantera.ck2cti
import rmgpy.chemkin
#import numpy as np
import subprocess
import csv
#import scipy
import copy
import os

################################################## functions #############################################

def N2D(N_reaction):
    '''
    Convert the NIST species in the reactions to RMG species, but keep the NIST kinetics
    '''
    
    D_reaction = copy.deepcopy(N_reaction)
    reactants = []
    for reactant in N_reaction.reactants:
        try:
            N_species_index = N_species_list.index(reactant)
            reactants.append(D_species_list[N2D_mapping[N_species_index]])
        except ValueError:
            if reactant in D_species_list:
                reactants.append(reactant)
        
    D_reaction.reactants = reactants
    
    products = []
    for product in N_reaction.products:
        try:
            N_species_index = N_species_list.index(product)
            products.append(D_species_list[N2D_mapping[N_species_index]])
        except ValueError:
            if product in D_species_list:
                products.append(product)
    D_reaction.products = products
    
    return D_reaction

def D2N(D_reaction):
    # takes in the D_reaction object to convert
    D_index = D_reaction_list.index(D_reaction)
    if D_index not in D2N_rxn_mapping.keys():
        # this reaction does not exist in Noras, so it will be deleted. return None
        return
    N_index = D2N_rxn_mapping[D_index]
    N_reaction = N_reaction_list[N_index]
    
    # convert the N model species in the N_reaction to D model species
    return N2D(N_reaction)



############################################### Load the models###########################################

#where to save the models
model_dir = 'flip_models'
os.makedirs(model_dir, exist_ok=True)


#Davids
full_path_D = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/compare/models/David/2-BTP/chemkin/'


D_chemkin_path = full_path_D + 'copies/copy_chem.inp'
D_dictionary_path = full_path_D + 'species_dictionary.txt'
D_transport_path = full_path_D + 'tran.dat'
D_cti_path = full_path_D + 'copies/copy_chem.cti'

D_species_list, D_reaction_list = rmgpy.chemkin.load_chemkin_file(D_chemkin_path, dictionary_path=D_dictionary_path, transport_path=D_transport_path)
D_gas = ct.Solution(D_cti_path)


# Noras
full_path_N = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/compare/models/Nora/2_BTP_seed/chemkin/'
N_cti_path = full_path_N + 'copies/copy_chem_130.cti'
N_chemkin_path = full_path_N + 'copies/copy_chem_130.inp'
N_dictionary_path = full_path_N + 'species_dictionary.txt'
N_transport_path = full_path_N + 'tran.dat'

N_gas = ct.Solution(N_cti_path)
N_dict = rmgpy.chemkin.load_species_dictionary(N_dictionary_path)
N_species_list, N_reaction_list = rmgpy.chemkin.load_chemkin_file(N_chemkin_path, dictionary_path=N_dictionary_path, transport_path=N_transport_path)


####################################### get the mapping between RMG and NIST models ####################################
# Species Diff
common_species = []
D2N_mapping = {}
N2D_mapping = {}
for i, D_sp in enumerate(D_species_list):
    for j, N_sp in enumerate(N_species_list):
        if D_sp.is_isomorphic(N_sp):
            D2N_mapping[i] = j
            N2D_mapping[j] = i
            common_species.append([D_sp, N_sp])
            break

# Reaction Diff
common_reactions = []
D2N_rxn_mapping = {}
N2D_rxn_mapping = {}
for i, D_rxn in enumerate(D_reaction_list):
    for j, N_rxn in enumerate(N_reaction_list):
        if D_rxn.is_isomorphic(N_rxn):
            D2N_rxn_mapping[i] = j
            N2D_rxn_mapping[j] = i
            common_reactions.append([D_rxn, N_rxn])
            break
print(f'{len(common_species)} common species')
print(f'{len(common_reactions)} common reactions')


common_reaction_index_D = []
common_reaction_index_N = []
for rD, rN in common_reactions: 
    common_reaction_index_D.append(rD.index)
    common_reaction_index_N.append(rN.index)
    
    

############################### convert the indicated reactions to use Noras kinetics ##########################################
print(len(D_reaction_list))
for reaction_to_flip in range(1732, len(D_reaction_list)):
    new_reaction_list = []
    deleted_duplicates = []
    for i in range(0, len(D_reaction_list)):
        if i == reaction_to_flip: 
            new_reaction = D2N(D_reaction_list[i]) #in this loop, will delete if not in Noras model 
            if new_reaction:
                new_reaction_list.append(new_reaction) 
            elif D_reaction_list[i].duplicate:
                deleted_duplicates.append(D_reaction_list[i])
        else:
            new_reaction_list.append(D_reaction_list[i])


    # get rid of duplicates
    for i, rxn in enumerate(new_reaction_list):
        if rxn.duplicate:
            duplicate_still_exists = False
            for j, rxn2 in enumerate(new_reaction_list):
                if rxn.is_isomorphic(rxn2) and rxn != rxn2:
                    duplicate_still_exists = True
                    break
            if not duplicate_still_exists:
                rxn.duplicate = False

    # mark reactions that are duplicates
    for i, rxn in enumerate(new_reaction_list):
        if not rxn.duplicate:
            duplicate = False
            for j, rxn2 in enumerate(new_reaction_list):
                if rxn.is_isomorphic(rxn2) and rxn != rxn2:
                    duplicate = True
                    break
            if duplicate:
                rxn.duplicate = True

    ########################################### save chemkin file ########################################################

    chemkin_file = os.path.join(model_dir, f'chem_{reaction_to_flip:04}.inp')
    rmgpy.chemkin.save_chemkin_file(chemkin_file, D_species_list, new_reaction_list, verbose=True, check_for_duplicates=True)
    subprocess.run(['ck2cti', f'--input={chemkin_file}', f'--transport={D_transport_path}', f'--output={model_dir}/chem_{reaction_to_flip:04}.cti'])

    #delete the input file because replaced with cti
    os.remove(chemkin_file)



