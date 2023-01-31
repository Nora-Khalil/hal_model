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


database(
thermoLibraries = thermolibs,
reactionLibraries = ['FFCM1(-)','halogens_pdep'],
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)


species(
    label = 'CH2F2',
    reactive = True,
    structure = SMILES('FCF')
)
    
species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)
    
    
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)

species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)

species(
    label = 'OH',
    reactive = True,
    structure = SMILES('[OH]')
)
  
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure=[(1.0,'bar'),(10.0,'bar')],
        nSims=12,
        initialMoleFractions={
        "CH2F2": [0.5,1.0],
        "O2": 1,
        "N2": 3.76,
        },
        terminationConversion={
        'CH2F2': 0.999,
        },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        sensitivity=['CH2F2','OH'],
        sensitivityThreshold=0.001,
        )
        
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=4,
    maximumOxygenAtoms=4,
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
    

