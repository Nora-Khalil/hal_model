species(
    label = '[CH2]C(C)(F)C(=C)OF(10094)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {12,S} {13,S}
8  C u0 p0 c0 {6,D} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-110.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,223.833,230.676],'cm^-1')),
        HinderedRotor(inertia=(0.225848,'amu*angstrom^2'), symmetry=1, barrier=(8.10266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224394,'amu*angstrom^2'), symmetry=1, barrier=(8.05506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594646,'amu*angstrom^2'), symmetry=1, barrier=(22.6184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808368,'amu*angstrom^2'), symmetry=1, barrier=(29.2849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273072,0.0891216,-0.000126534,1.06148e-07,-3.67202e-11,-13175.9,27.5221], Tmin=(100,'K'), Tmax=(761.079,'K')), NASAPolynomial(coeffs=[8.69956,0.0398358,-1.95457e-05,3.8019e-09,-2.66619e-13,-14313.8,-9.87916], Tmin=(761.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P)"""),
)

species(
    label = 'HF(38)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-281.113,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(20.0062,'amu')),
        LinearRotor(inertia=(0.809097,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([4113.43],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.0064,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2743.78,'J/mol'), sigma=(3.148,'angstroms'), dipoleMoment=(1.92,'De'), polarizability=(2.46,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43657,0.000486021,-1.2524e-06,1.36475e-09,-4.09574e-13,-33800.1,1.20682], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[2.7813,0.00103959,-2.41735e-07,2.68416e-11,-1.09766e-15,-33504.2,5.0197], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-281.113,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""HF""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C(=C)F(1340)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-49.3492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0540617,'amu*angstrom^2'), symmetry=1, barrier=(25.2827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2755.66,'J/mol'), sigma=(4.80499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=430.43 K, Pc=56.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93971,0.00348861,8.08985e-05,-1.51352e-07,8.61066e-11,-5933.05,7.77268], Tmin=(10,'K'), Tmax=(574.141,'K')), NASAPolynomial(coeffs=[2.51346,0.0273693,-1.79226e-05,5.69534e-09,-6.965e-13,-5999.1,11.8607], Tmin=(574.141,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-49.3492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[CH2]C(DC)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45064e-06,-3.57997e-09,7.56176e-13,50400.6,-0.411759], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62651,0.00394757,-1.49921e-06,2.54533e-10,-1.62951e-14,50691.7,6.78356], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(F)C(=C)OF(593)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-46.8193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,174,267,591,721,1107,1278,1348,3273,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0928072,'amu*angstrom^2'), symmetry=1, barrier=(2.13382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0931087,'amu*angstrom^2'), symmetry=1, barrier=(2.14075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44631,'amu*angstrom^2'), symmetry=1, barrier=(33.2535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15071,0.0684572,-9.10202e-05,7.08301e-08,-2.30717e-11,-5533.8,23.7922], Tmin=(100,'K'), Tmax=(740.653,'K')), NASAPolynomial(coeffs=[8.02831,0.0313164,-1.58062e-05,3.13404e-09,-2.23148e-13,-6552.65,-7.32385], Tmin=(740.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.8193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(H)(H))"""),
)

species(
    label = 'C=C(OF)[C](F)CC(10098)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {4,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-139.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812673,0.0756435,-7.67897e-05,4.39633e-08,-1.07065e-11,-16667.4,25.6582], Tmin=(100,'K'), Tmax=(965.498,'K')), NASAPolynomial(coeffs=[9.89249,0.0380263,-1.83475e-05,3.6095e-09,-2.57531e-13,-18420.8,-17.8278], Tmin=(965.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCdCsF1s)"""),
)

species(
    label = 'C=C(C[C](C)F)OF(10099)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-118.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876164,0.0775058,-8.43431e-05,3.95022e-08,1.73413e-12,-14192,28.0373], Tmin=(100,'K'), Tmax=(571.072,'K')), NASAPolynomial(coeffs=[8.1751,0.0395344,-1.91546e-05,3.73328e-09,-2.63329e-13,-15240.1,-4.96454], Tmin=(571.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH2]C(C)(F)C(=O)CF(10100)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  C u1 p0 c0 {4,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-459.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0568554,0.0949838,-0.000143902,1.24312e-07,-4.31815e-11,-55127.3,27.4153], Tmin=(100,'K'), Tmax=(799.198,'K')), NASAPolynomial(coeffs=[9.07882,0.0392074,-1.92815e-05,3.72778e-09,-2.59411e-13,-56230.1,-11.966], Tmin=(799.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C)=C(CF)OF(10101)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u1 p0 c0 {6,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-129.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703017,0.0785518,-8.78324e-05,5.72418e-08,-1.58436e-11,-15518.2,25.5262], Tmin=(100,'K'), Tmax=(859.371,'K')), NASAPolynomial(coeffs=[9.19963,0.0390015,-1.87949e-05,3.6821e-09,-2.61653e-13,-16978.5,-14.1769], Tmin=(859.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(F)=C(CC)OF(10102)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {1,S} {6,D} {8,S}
8  C u1 p0 c0 {7,S} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-145.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770547,0.0775522,-8.84495e-05,6.0845e-08,-1.80273e-11,-17337.7,27.2959], Tmin=(100,'K'), Tmax=(801.551,'K')), NASAPolynomial(coeffs=[8.09711,0.0409916,-2.00335e-05,3.94397e-09,-2.80794e-13,-18512.3,-6.43008], Tmin=(801.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(OF)[C](C)F(592)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {6,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-115.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,346,659,817,1284,350,440,435,1725,2950,3100,1380,975,1025,1650,279.648,279.65],'cm^-1')),
        HinderedRotor(inertia=(0.113959,'amu*angstrom^2'), symmetry=1, barrier=(6.32444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720731,'amu*angstrom^2'), symmetry=1, barrier=(39.9972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303492,'amu*angstrom^2'), symmetry=1, barrier=(16.8424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42814,0.0612443,-6.43496e-05,3.77745e-08,-9.37615e-12,-13828.8,21.2021], Tmin=(100,'K'), Tmax=(951.225,'K')), NASAPolynomial(coeffs=[8.93482,0.0296774,-1.45704e-05,2.88609e-09,-2.06644e-13,-15256.9,-14.6378], Tmin=(951.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCdCsF1s)"""),
)

species(
    label = 'H(5)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,25474.2,-0.444972], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(C)(F)C(=C)OF-2(10103)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  C u2 p0 c0 {4,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (132.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.158464,0.0917439,-0.000136652,1.14875e-07,-3.91462e-11,16050.5,27.1504], Tmin=(100,'K'), Tmax=(774.042,'K')), NASAPolynomial(coeffs=[9.93125,0.0358249,-1.77911e-05,3.4618e-09,-2.42229e-13,14699.8,-16.4464], Tmin=(774.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC1(F)CC[C]1OF(10104)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-112.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33341,0.0605618,-4.17801e-05,1.44534e-08,-2.09273e-12,-13413,23.9792], Tmin=(100,'K'), Tmax=(1524.26,'K')), NASAPolynomial(coeffs=[11.5321,0.0337984,-1.54429e-05,2.93449e-09,-2.03483e-13,-16522.2,-29.5228], Tmin=(1524.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs-Cs) + radical(C2CsJO)"""),
)

species(
    label = '[CH2]C1(OF)CC1(C)F(10105)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {5,S} {14,S} {15,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-40.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218889,0.079334,-7.6344e-05,3.77615e-08,-7.46964e-12,-4789.06,24.5228], Tmin=(100,'K'), Tmax=(1217.71,'K')), NASAPolynomial(coeffs=[16.3242,0.0264301,-1.11757e-05,2.0833e-09,-1.44761e-13,-8711.37,-56.3483], Tmin=(1217.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsOs) + group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cs-Cs(F)(C)-Cs) + radical(CJC(C)2O)"""),
)

species(
    label = 'F(37)',
    structure = adjacencyList("""multiplicity 2
1 F u1 p3 c0
"""),
    E0 = (72.8916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.9984,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.158,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(C)C(=C)OF(10106)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (0.243838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,265.315,265.755],'cm^-1')),
        HinderedRotor(inertia=(0.304353,'amu*angstrom^2'), symmetry=1, barrier=(15.181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303242,'amu*angstrom^2'), symmetry=1, barrier=(15.181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304671,'amu*angstrom^2'), symmetry=1, barrier=(15.1755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850476,0.0692255,-6.88915e-05,3.68633e-08,-8.03185e-12,142.817,21.614], Tmin=(100,'K'), Tmax=(1100.49,'K')), NASAPolynomial(coeffs=[12.4103,0.027209,-1.16226e-05,2.17075e-09,-1.5081e-13,-2401.52,-35.2627], Tmin=(1100.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.243838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH3(19)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C(F)C(=C)OF(3267)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-157.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,280,518,736,852,873,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,281.521,282.621],'cm^-1')),
        HinderedRotor(inertia=(0.150947,'amu*angstrom^2'), symmetry=1, barrier=(8.33212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420798,'amu*angstrom^2'), symmetry=1, barrier=(24.0743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39908,0.0604488,-6.75445e-05,4.02685e-08,-9.82151e-12,-18822.6,19.6222], Tmin=(100,'K'), Tmax=(983.658,'K')), NASAPolynomial(coeffs=[10.4433,0.0236712,-1.14622e-05,2.25945e-09,-1.61445e-13,-20601.9,-23.8619], Tmin=(983.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]OF(202)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (235.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2950,3100,1380,975,1025,1650,1685,370,1858.56],'cm^-1')),
        HinderedRotor(inertia=(0.190472,'amu*angstrom^2'), symmetry=1, barrier=(8.45471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27376,0.0120477,-4.21846e-06,6.03973e-11,9.31488e-14,28267.1,10.4556], Tmin=(100,'K'), Tmax=(2667.03,'K')), NASAPolynomial(coeffs=[19.9694,-0.00396407,5.52785e-07,-7.38522e-11,6.52159e-15,17217.5,-85.6832], Tmin=(2667.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=CJO)"""),
)

species(
    label = 'C=C(C)F(497)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-203.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,323,467,575,827,1418,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.291017,'amu*angstrom^2'), symmetry=1, barrier=(6.69106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2755.66,'J/mol'), sigma=(4.80499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=430.43 K, Pc=56.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96409,0.00231263,9.07676e-05,-1.79095e-07,1.14016e-10,-24461.7,8.35177], Tmin=(10,'K'), Tmax=(461.451,'K')), NASAPolynomial(coeffs=[1.43034,0.0300497,-1.81635e-05,5.39498e-09,-6.25069e-13,-24289.3,17.95], Tmin=(461.451,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-203.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CDC(C)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(C)=C([CH2])OF(10107)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {4,D} {7,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (184.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180.002],'cm^-1')),
        HinderedRotor(inertia=(0.642516,'amu*angstrom^2'), symmetry=1, barrier=(16.2163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835099,'amu*angstrom^2'), symmetry=1, barrier=(16.3463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706256,'amu*angstrom^2'), symmetry=1, barrier=(16.2613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259294,'amu*angstrom^2'), symmetry=1, barrier=(74.036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839188,0.070149,-7.18968e-05,4.01709e-08,-9.1619e-12,22328.4,24.5731], Tmin=(100,'K'), Tmax=(1053.25,'K')), NASAPolynomial(coeffs=[11.9253,0.0280467,-1.19366e-05,2.21861e-09,-1.53569e-13,19993.1,-29.4862], Tmin=(1053.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(C)(F)C(=C)[O](10108)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-130.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.682247,'amu*angstrom^2'), symmetry=1, barrier=(15.6862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682286,'amu*angstrom^2'), symmetry=1, barrier=(15.6871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682317,'amu*angstrom^2'), symmetry=1, barrier=(15.6878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938203,0.0711793,-7.97749e-05,5.08169e-08,-1.34387e-11,-15602,24.1803], Tmin=(100,'K'), Tmax=(907.395,'K')), NASAPolynomial(coeffs=[9.91908,0.0315899,-1.43307e-05,2.73523e-09,-1.91642e-13,-17231.8,-18.2746], Tmin=(907.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_P)"""),
)

species(
    label = '[O]F(126)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(C)(F)[C]=C(10109)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {6,D} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (164.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.676955,'amu*angstrom^2'), symmetry=1, barrier=(15.5645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677009,'amu*angstrom^2'), symmetry=1, barrier=(15.5658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.676796,'amu*angstrom^2'), symmetry=1, barrier=(15.5609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9877,0.0526007,-2.94259e-05,-3.85828e-08,5.22665e-11,19817.8,19.9582], Tmin=(100,'K'), Tmax=(483.492,'K')), NASAPolynomial(coeffs=[5.87425,0.0344177,-1.63581e-05,3.16349e-09,-2.22116e-13,19278.7,2.34355], Tmin=(483.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(F)=C([CH2])OF(3234)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (36.4894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,271,519,563,612,1379,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00367833,'amu*angstrom^2'), symmetry=1, barrier=(10.8951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84841,'amu*angstrom^2'), symmetry=1, barrier=(42.4985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0143197,'amu*angstrom^2'), symmetry=1, barrier=(42.4894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05261,0.0704445,-0.000105556,8.91382e-08,-3.03705e-11,4489.41,24.0394], Tmin=(100,'K'), Tmax=(789.832,'K')), NASAPolynomial(coeffs=[8.43861,0.0275159,-1.35394e-05,2.61708e-09,-1.82208e-13,3494.95,-8.76053], Tmin=(789.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.4894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(O)CJ) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = '[CH2][C](C)F(9483)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (68.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.222019,'amu*angstrom^2'), symmetry=1, barrier=(5.10466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0018878,'amu*angstrom^2'), symmetry=1, barrier=(9.28569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69728,0.022342,-9.99926e-06,1.80524e-09,-1.13534e-13,8240.19,13.0709], Tmin=(100,'K'), Tmax=(3036.29,'K')), NASAPolynomial(coeffs=[27.4193,-0.00471389,1.29434e-06,-2.19373e-10,1.5696e-14,-8099.13,-130.906], Tmin=(3036.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2]C([CH2])(F)C(=C)OF(10110)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  C u0 p0 c0 {5,D} {13,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (94.6363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,210.72,1608.13],'cm^-1')),
        HinderedRotor(inertia=(0.292504,'amu*angstrom^2'), symmetry=1, barrier=(9.20668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292425,'amu*angstrom^2'), symmetry=1, barrier=(9.20767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96588,'amu*angstrom^2'), symmetry=1, barrier=(30.4405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96698,'amu*angstrom^2'), symmetry=1, barrier=(30.4408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.256525,0.0905313,-0.000139528,1.21891e-07,-4.26904e-11,11509,28.6719], Tmin=(100,'K'), Tmax=(798.922,'K')), NASAPolynomial(coeffs=[8.70616,0.0371613,-1.85502e-05,3.60426e-09,-2.51357e-13,10512,-7.9859], Tmin=(798.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.6363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cs_P)"""),
)

species(
    label = '[CH]=C(OF)C([CH2])(C)F(10111)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {12,S} {13,S}
8  C u1 p0 c0 {6,D} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (136.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.82315,'amu*angstrom^2'), symmetry=1, barrier=(18.9258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.446149,'amu*angstrom^2'), symmetry=1, barrier=(10.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49792,'amu*angstrom^2'), symmetry=1, barrier=(34.4402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.446159,'amu*angstrom^2'), symmetry=1, barrier=(10.2581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135223,0.0934388,-0.000146118,1.27437e-07,-4.42895e-11,16546.5,28.4171], Tmin=(100,'K'), Tmax=(806.305,'K')), NASAPolynomial(coeffs=[9.24117,0.0364705,-1.8196e-05,3.52705e-09,-2.45344e-13,15461.5,-11.1759], Tmin=(806.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=C)C(=C)OF(9294)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {12,S} {13,S}
7  C u0 p0 c0 {4,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (151.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,360.91,360.914],'cm^-1')),
        HinderedRotor(inertia=(0.209799,'amu*angstrom^2'), symmetry=1, barrier=(19.3923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20979,'amu*angstrom^2'), symmetry=1, barrier=(19.3923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.342605,'amu*angstrom^2'), symmetry=1, barrier=(31.6673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622777,0.0694365,-6.98601e-05,3.58556e-08,-7.26651e-12,18376.2,22.0159], Tmin=(100,'K'), Tmax=(1201.06,'K')), NASAPolynomial(coeffs=[15.5426,0.0197471,-7.80246e-06,1.40923e-09,-9.64465e-14,14792.3,-52.697], Tmin=(1201.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(OF)C(C)(C)F(10112)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u1 p0 c0 {7,D} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-68.7603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.768032,'amu*angstrom^2'), symmetry=1, barrier=(17.6586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767091,'amu*angstrom^2'), symmetry=1, barrier=(17.6369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768738,'amu*angstrom^2'), symmetry=1, barrier=(17.6748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767806,'amu*angstrom^2'), symmetry=1, barrier=(17.6534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146007,0.0921046,-0.000133424,1.12134e-07,-3.85267e-11,-8138.18,25.9012], Tmin=(100,'K'), Tmax=(771.878,'K')), NASAPolynomial(coeffs=[9.25289,0.0391115,-1.91712e-05,3.71971e-09,-2.60182e-13,-9371.29,-14.557], Tmin=(771.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.7603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C(C)(F)CF(10113)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-510.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530394,0.0800385,-8.76794e-05,5.28643e-08,-1.3123e-11,-61219.5,25.166], Tmin=(100,'K'), Tmax=(966.181,'K')), NASAPolynomial(coeffs=[11.7856,0.0334425,-1.53397e-05,2.95032e-09,-2.07911e-13,-63394.4,-28.7467], Tmin=(966.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(OF)[C](C)CF(10114)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u1 p0 c0 {4,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-141.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763059,0.0682407,-5.57946e-05,2.34715e-08,-4.03496e-12,-16874.3,27.3538], Tmin=(100,'K'), Tmax=(1362.12,'K')), NASAPolynomial(coeffs=[14.0815,0.0291294,-1.27241e-05,2.39121e-09,-1.6592e-13,-20502.5,-41.016], Tmin=(1362.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJ(C)CO)"""),
)

species(
    label = 'N2',
    structure = adjacencyList("""1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
"""),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0137,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = adjacencyList("""1 Ne u0 p4 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1801,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (-93.0747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (353.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (27.8636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (27.8636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (33.3601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (169.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (243.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (246.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (325.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-3.93479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-59.8217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (118.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-39.3859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (37.0002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (238.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-41.9396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (248.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (168.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (284.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (287.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (329.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (97.4504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-43.2964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-74.1381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (36.5911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['HF(38)', 'CH2CO(28)', '[CH2]C(=C)F(1340)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(36.3801,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[CH2]C(F)C(=C)OF(593)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.16448e+06,'m^3/(mol*s)'), n=-0.144618, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['C=C(OF)[C](F)CC(10098)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 2 used for cCs(-R!HR!H)CJ;CsJ-HH;CH3
Exact match found for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['C=C(C[C](C)F)OF(10099)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['[CH2]C(C)(F)C(=O)CF(10100)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(162.815,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['[CH2]C(C)=C(CF)OF(10101)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(298.698,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['[CH2]C(F)=C(CC)OF(10102)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3431.13,'s^-1'), n=2.77015, Ea=(372.584,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_2R!H->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(18)', 'C=C(OF)[C](C)F(592)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH]C(C)(F)C(=C)OF-2(10103)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['CC1(F)CC[C]1OF(10104)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.82e+08,'s^-1'), n=0.91, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['[CH2]C1(OF)CC1(C)F(10105)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(69.6331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 68.7 to 69.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C=C(C)C(=C)OF(10106)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(63.8763,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(19)', 'C=C(F)C(=C)OF(3267)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.97336e-14,'m^3/(mol*s)'), n=5.33869, Ea=(0.527887,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.259264086083826, var=5.650144167053507, Tref=1000.0, N=67, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]OF(202)', 'C=C(C)F(497)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1107.8,'m^3/(mol*s)'), n=0.867872, Ea=(23.975,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.042098679164349666, var=0.7527123938779839, Tref=1000.0, N=88, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[CH2]C(C)=C([CH2])OF(10107)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', '[CH2]C(C)(F)C(=C)[O](10108)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(34.625,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]F(126)', '[CH2]C(C)(F)[C]=C(10109)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(19)', '[CH2]C(F)=C([CH2])OF(3234)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.29446e+18,'m^3/(mol*s)'), n=-4.19701, Ea=(14.7964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OF(202)', '[CH2][C](C)F(9483)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -20.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[CH2]C([CH2])(F)C(=C)OF(10110)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.3422e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C
Multiplied by reaction path degeneracy 2.0
Ea raised from -0.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[CH]=C(OF)C([CH2])(C)F(10111)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[CH2]C(=C)C(=C)OF(9294)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(50.0704,'m^3/(mol*s)'), n=1.25316, Ea=(245.666,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(OF)C(C)(C)F(10112)'],
    products = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['C=C([O])C(C)(F)CF(10113)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(55.3167,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C)(F)C(=C)OF(10094)'],
    products = ['C=C(OF)[C](C)CF(10114)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(166.046,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3523',
    isomers = [
        '[CH2]C(C)(F)C(=C)OF(10094)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', '[CH2]C(=C)F(1340)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3523',
    Tmin = (300,'K'),
    Tmax = (2500,'K'),
    Tcount = 8,
    Tlist = ([302.558,324.028,372.925,464.512,632.697,950.724,1545.17,2335.46],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.0125282,0.0667467,1,14.982,79.8202],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

