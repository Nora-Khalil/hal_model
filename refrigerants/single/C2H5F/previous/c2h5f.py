database(
    thermoLibraries = [
    'primaryThermoLibrary','FFCM1(-)',
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
    ],
    
    reactionLibraries = ['FFCM1(-)','halogens_pdep'],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['default','halogens','Disproportionation-Y'],
    frequenciesLibraries = ['halogens_G4'],
    kineticsEstimator='rate rules',
)

species(
    label='C2H5F',
    reactive=True,
    structure = SMILES('CCF'),
)

species(
    label='N2',
    reactive=False,
    structure=adjacencyList("""
    1 N u0 p1 c0 {2,T}
    2 N u0 p1 c0 {1,T}
    """),
)

species(
    label='O2',
    structure=SMILES("[O][O]"),
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

species(
    label = 'CH4',
    reactive = True,
    structure = SMILES('C')
)
simpleReactor(
    temperature=[(1000, 'K'), (2000, 'K')],
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        "C2H5F": 1,
        "N2": 11.289,
        "O2": 3,
    },
    #terminationConversion={
      #  'C2H6': .99,
   # },
    #terminationRateRatio=1e-8,
    terminationTime=(.1, 's'), 
)

simulator(
    atol=1e-16,
    rtol=1e-8,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=6,
    maximumOxygenAtoms=6,
    #maximumHeavyAtoms=24,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)


model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=5e5,
    filterReactions=True,
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

options(
    units='si', 
    generateOutputHTML=True,
    generateSeedEachIteration = True,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=True,
    keepIrreversible = True,
)

