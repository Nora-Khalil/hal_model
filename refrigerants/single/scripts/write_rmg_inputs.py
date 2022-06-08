import os
from rmgpy.molecule import Molecule

directory = '../'

species = {
  
  'C2H6' : 'CC',
  'CH4': 'C',
  'CH2F2': 'FCF',
  'CH3F': 'CF',
  'C2H5F': 'CCF',
  'CH2FCH2F': 'FCCF',
  'CH2FCHF2' : 'FCC(F)F',
  'CH3CF3': 'CC(F)(F)F',
  'CH3CHF2': 'CC(F)F',
}

input_file= '''
thermolibs = [
'primaryThermoLibrary',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'Fluorine',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]

thermolibs_Creg = [
'primaryThermoLibrary',
'FFCM1(-)',
'DFT_QCI_thermo',
]





'''

#def add_database(kind):
def add_database(label): 
    if (label == 'CH4') or (label == 'C2H6'):
       block = '''
database(
thermoLibraries = thermolibs_Creg,
reactionLibraries = ['FFCM1(-)'],
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = 'default',
kineticsEstimator = 'rate rules',
)
'''
    else:    
       block = '''
database(
thermoLibraries = thermolibs,
reactionLibraries = ['halogens_pdep'],
seedMechanisms = ['FFCM1(-)'],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)

'''
    return block

def make_species_block(label,smiles,reactive=True):
    block = \
    f'''
species(
    label = '{label}',
    reactive = {reactive},
    structure = SMILES('{smiles}')
)
    '''
    return block

def make_reactor(label,mol):

   # Cs = mol.get_num_atoms('C')
   # low = str(round(0.15/Cs,4))
   # high = str(round(0.30/Cs,4))

    reactor = \
    '''
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure= (1.0,'bar'),
        nSims=12,
        initialMoleFractions={
        "halogen": [0.5,1.0],
        "O2": 1,
        "N2": 3.76,
        },
        # terminationConversion={
        # 'halogen': 0.999,
        # },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        #sensitivity=['halogen','OH'],
        #sensitivityThreshold=0.001,
        )
        '''.replace('halogen',label) #.replace("low",low).replace("high",high)
    return reactor

def write_settings(mol=None):
    if mol:
        C = mol.get_num_atoms('C')
        if C == 1:
            Cs = 4
            Os = 4
        elif C == 2:
            Cs = 6
            Os = 4
        elif C == 3:
            Cs = 6
            Os = 4
        else:
            Cs = 2*C
            Os = 2*C
    else:
        Cs = 8
        Os = 6

    settings = \
    f"""
simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms={Cs},
    maximumOxygenAtoms={Os},
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)

options(
    units = "si",
    generateSeedEachIteration = True,
    generateOutputHTML = True,
    generatePlots = True,
    saveSimulationProfiles = True,
    saveEdgeSpecies = False,
    keepIrreversible = True,
    verboseComments = False,
)
    """
    return settings


def make_bash_script(label):
    bash_script = """#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=55:00:00
#SBATCH --job-name=label
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py
""".replace('label',label)

    return bash_script

model_loose = \
"""
model(
    toleranceMoveToCore = 0.5,
    toleranceInterruptSimulation = 0.5,
    toleranceKeepInEdge = 0.01,
    maximumEdgeSpecies = 5e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 100,
    minSpeciesExistIterationsForPrune = 2,
)
"""

model_medium = \
"""
model(
    toleranceMoveToCore = 0.2,
    toleranceInterruptSimulation = 0.2,
    maximumEdgeSpecies = 5e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)
"""

model_tight = \
"""
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)
"""

pdep = \
"""
pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)
"""


if __name__ == "__main__":

    # write inputs for individual refrigerants
    #for kind,species_dict in species.items():
        for label,smiles in species.items():
            d = os.path.join(directory,label,'input.py')
            if os.path.exists(d):
                continue
            print(label)
            input0 = input_file
            input0 += add_database(label)
            species_block = ""
            species_block += make_species_block(label,smiles)
            species_block += make_species_block('O2','[O][O]')
            species_block += make_species_block('H2O','O')
            species_block += make_species_block('N2','N#N',False)
            #adding in extra species that we know we'll need below
            if (label == 'C2H6') or (label == 'CH4'):
               species_block += make_species_block('C:','[C]') 
               species_block += make_species_block('CH','[CH]')
            else: 
               species_block += make_species_block('CH4','C')
            input0 += species_block
            mol = Molecule(smiles=smiles)
            input0 += make_reactor(label,mol)
            input0 += model_tight
            input0 += pdep
            input0 += write_settings(mol)
            bash_script = make_bash_script(label)
            d = os.path.join(directory,label)#'input.py')
            os.makedirs(d,exist_ok=True)
            with open(os.path.join(d,'input.py'),'w') as f:
                for l in input0:
                    f.write(l)
            with open(os.path.join(d,'run.sh'),'w') as f:
                for l in bash_script:
                    f.write(l)





