

database(
    thermoLibraries = ['primaryThermoLibrary', 'DFT_QCI_thermo','FFCM1(-)'],
    reactionLibraries = ['FFCM1(-)'],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

species(
    label='C2H6',
    reactive=True,
    structure=adjacencyList("""
    1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
    2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}
    6 H u0 p0 c0 {2,S}
    7 H u0 p0 c0 {2,S}
    8 H u0 p0 c0 {2,S}
    """),
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
    reactive=True,
    structure=SMILES("[O][O]"),
)

species(
    label='C:',
    reactive=True,
    structure=SMILES("[C]"),
)

species(
    label='CH',
    reactive=True,
    structure=SMILES("[CH]"),
)

simpleReactor(
    temperature=[(1000, 'K'), (2000, 'K')],
    # specifies reaction pressure with units
    pressure=(1.0, 'bar'),
    # list initial mole fractions of compounds using the label from the 'species' label.
    # RMG will normalize if sum/=1
    nSims = 12,
    initialMoleFractions={
        "C2H6": 1,
        "N2": 13.16,
        "O2": 3.5,
    },
   # terminationConversion={
   #             'C2H6': 0.99,
   # },
    terminationRateRatio=1e-8,
    terminationTime=(1e6, 's'), 
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    #maximumCarbonAtoms=7,
    #maximumOxygenAtoms=6,
    #maximumHeavyAtoms=24,
    #maximumRadicalElectrons=2,
    #maximumSingletCarbenes=1,
    #maximumCarbeneRadicals=0,
    #allowSingletO2 = True,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    #maximumAtoms=16,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=5e5,
    filterReactions=True,
    filterThreshold = 5e8, 
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)
options(
    units='si', 
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)




