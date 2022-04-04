species(
    label = 'F[CH]C(CCF)OF(1953)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-365.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,299.722,299.726,299.726,299.728],'cm^-1')),
        HinderedRotor(inertia=(0.215942,'amu*angstrom^2'), symmetry=1, barrier=(13.7659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134242,'amu*angstrom^2'), symmetry=1, barrier=(8.55772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13424,'amu*angstrom^2'), symmetry=1, barrier=(8.55774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516179,'amu*angstrom^2'), symmetry=1, barrier=(32.9064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152927,0.0920766,-0.000138265,1.1707e-07,-4.01812e-11,-43813.9,28.1452], Tmin=(100,'K'), Tmax=(770.696,'K')), NASAPolynomial(coeffs=[9.86136,0.0360054,-1.80721e-05,3.53291e-09,-2.47847e-13,-45141.5,-15.0687], Tmin=(770.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
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
    label = 'O=C[CH]F(338)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-191.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.21522,'amu*angstrom^2'), symmetry=1, barrier=(35.2035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2974.27,'J/mol'), sigma=(4.91649,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.57 K, Pc=56.79 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96722,0.00209998,6.04576e-05,-1.2409e-07,8.01311e-11,-23018.4,8.77531], Tmin=(10,'K'), Tmax=(478.869,'K')), NASAPolynomial(coeffs=[2.63637,0.0192853,-1.23829e-05,3.78104e-09,-4.41625e-13,-22960.5,13.4894], Tmin=(478.869,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-191.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CHF(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[CH]C(CF)OF(1957)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {5,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-338.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.817023,'amu*angstrom^2'), symmetry=1, barrier=(18.785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816807,'amu*angstrom^2'), symmetry=1, barrier=(18.78,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817003,'amu*angstrom^2'), symmetry=1, barrier=(18.7845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59147,0.0322212,0.00016733,-6.25008e-07,5.9444e-10,-40687.9,12.3899], Tmin=(10,'K'), Tmax=(396.037,'K')), NASAPolynomial(coeffs=[8.19887,0.0343097,-2.4743e-05,8.32438e-09,-1.04771e-12,-41434.2,-10.3848], Tmin=(396.037,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-338.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""F[CH]C(CF)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF(40)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (138.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'CC([CH]F)OF(1958)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {4,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-161.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.200453,'amu*angstrom^2'), symmetry=1, barrier=(4.6088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200156,'amu*angstrom^2'), symmetry=1, barrier=(4.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22416,'amu*angstrom^2'), symmetry=1, barrier=(28.1459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.72824,0.0242375,0.000178231,-7.45183e-07,8.64754e-10,-19443.9,11.181], Tmin=(10,'K'), Tmax=(309.414,'K')), NASAPolynomial(coeffs=[5.24309,0.0351916,-2.2915e-05,7.18093e-09,-8.62324e-13,-19683.8,3.28758], Tmin=(309.414,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-161.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""CC([CH]F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FCCC(F)[CH]OF(1959)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
8  C u1 p0 c0 {4,S} {6,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-387.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293168,0.0871823,-0.000110912,7.544e-08,-2.07867e-11,-46419.3,26.8405], Tmin=(100,'K'), Tmax=(880.12,'K')), NASAPolynomial(coeffs=[12.6366,0.0310832,-1.5302e-05,3.01754e-09,-2.14926e-13,-48592,-31.1333], Tmin=(880.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFHH) + radical(CCsJO)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(CCF)OF-2(1960)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
7  C u2 p0 c0 {5,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (77.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936087,0.0701555,-7.60301e-05,4.40555e-08,-1.0427e-11,9382.62,25.0835], Tmin=(100,'K'), Tmax=(1013.7,'K')), NASAPolynomial(coeffs=[11.7614,0.027441,-1.28265e-05,2.49076e-09,-1.76603e-13,7187.81,-27.2901], Tmin=(1013.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(230)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'FCC[CH]OF(514)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {4,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-167.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.774123,'amu*angstrom^2'), symmetry=1, barrier=(17.7986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773701,'amu*angstrom^2'), symmetry=1, barrier=(17.7889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773398,'amu*angstrom^2'), symmetry=1, barrier=(17.7819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13181,0.0676732,-8.98214e-05,6.31273e-08,-1.78485e-11,-20079.4,19.8973], Tmin=(100,'K'), Tmax=(860.822,'K')), NASAPolynomial(coeffs=[10.8446,0.0225393,-1.11728e-05,2.21612e-09,-1.58289e-13,-21751.6,-25.5056], Tmin=(860.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
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
    label = 'F[C]C(CCF)OF(1961)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
8  C u2 p0 c0 {2,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-189.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,90,150,250,247,323,377,431,433,1065,1274,197.715,450.628,702.659,2417.91,2804.05,3460.57],'cm^-1')),
        HinderedRotor(inertia=(0.0162379,'amu*angstrom^2'), symmetry=1, barrier=(0.384494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162379,'amu*angstrom^2'), symmetry=1, barrier=(0.384494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162379,'amu*angstrom^2'), symmetry=1, barrier=(0.384494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162379,'amu*angstrom^2'), symmetry=1, barrier=(0.384494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.551761,0.0820903,-0.000109307,7.94499e-08,-2.35569e-11,-22656.6,26.0208], Tmin=(100,'K'), Tmax=(817.933,'K')), NASAPolynomial(coeffs=[10.9973,0.0310087,-1.56306e-05,3.09929e-09,-2.20875e-13,-24365.4,-22.2738], Tmin=(817.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(CsCCl_triplet)"""),
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
    label = 'FC=CCCF(232)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {2,S} {5,D} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-397.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.27409,'amu*angstrom^2'), symmetry=1, barrier=(6.30186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407872,'amu*angstrom^2'), symmetry=1, barrier=(9.37777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0872,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51497,0.0460759,-8.16879e-05,1.46099e-07,-1.06698e-10,-47781.5,10.8748], Tmin=(10,'K'), Tmax=(460.173,'K')), NASAPolynomial(coeffs=[2.81637,0.0403103,-2.43065e-05,7.06539e-09,-7.94085e-13,-47591.9,15.0648], Tmin=(460.173,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-397.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FCDCCCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=COF(1541)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-182.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.052636,'amu*angstrom^2'), symmetry=1, barrier=(19.8022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.901,0.00613816,9.01112e-05,-2.02909e-07,1.32775e-10,-21918.9,10.3718], Tmin=(10,'K'), Tmax=(529.501,'K')), NASAPolynomial(coeffs=[4.44857,0.0228458,-1.62676e-05,5.37178e-09,-6.64336e-13,-22269.1,5.31897], Tmin=(529.501,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-182.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""FCDCOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2F-CH2(62)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-75.6489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0750522,'amu*angstrom^2'), symmetry=1, barrier=(7.14957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0516,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(2,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99845,0.00714341,1.60323e-05,-2.0007e-08,6.54682e-12,-9101.39,8.71379], Tmin=(10,'K'), Tmax=(1078.03,'K')), NASAPolynomial(coeffs=[3.84792,0.0159784,-7.77703e-06,1.83859e-09,-1.70814e-13,-9549.86,7.22073], Tmin=(1078.03,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-75.6489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[CH2]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C(CCF)OF(1962)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {2,S} {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-427.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,194,682,905,1196,1383,3221,292.468,292.522,1991.51],'cm^-1')),
        HinderedRotor(inertia=(0.131732,'amu*angstrom^2'), symmetry=1, barrier=(7.99951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131772,'amu*angstrom^2'), symmetry=1, barrier=(7.99988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547057,'amu*angstrom^2'), symmetry=1, barrier=(33.2173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32921,0.0677295,-6.90114e-05,2.04327e-08,1.36453e-11,-51281.7,24.6748], Tmin=(100,'K'), Tmax=(532.72,'K')), NASAPolynomial(coeffs=[7.22121,0.0366515,-1.8567e-05,3.68684e-09,-2.6273e-13,-52096.2,-1.79318], Tmin=(532.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-427.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[CH2]CC([CH]F)OF(1963)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {4,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (28.7003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,301.942,301.949,301.951],'cm^-1')),
        HinderedRotor(inertia=(0.108622,'amu*angstrom^2'), symmetry=1, barrier=(7.02766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.527381,'amu*angstrom^2'), symmetry=1, barrier=(34.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108621,'amu*angstrom^2'), symmetry=1, barrier=(7.02764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00337749,'amu*angstrom^2'), symmetry=1, barrier=(34.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666225,0.0797046,-0.000117358,9.90078e-08,-3.39424e-11,3565.79,27.151], Tmin=(100,'K'), Tmax=(779.439,'K')), NASAPolynomial(coeffs=[8.73298,0.0326865,-1.60572e-05,3.11254e-09,-2.17376e-13,2479.01,-8.66119], Tmin=(779.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.7003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(Cs-CsHHH) + radical(Csj(Cs-O2sCsH)(F1s)(H)) + radical(RCCJ)"""),
)

species(
    label = '[O]C([CH]F)CCF(1964)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-269.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,319.764,322.095,325.323,327.07],'cm^-1')),
        HinderedRotor(inertia=(0.0952345,'amu*angstrom^2'), symmetry=1, barrier=(7.12549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889042,'amu*angstrom^2'), symmetry=1, barrier=(7.10071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312529,'amu*angstrom^2'), symmetry=1, barrier=(23.9918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09759,0.0713119,-8.53659e-05,5.08545e-08,-7.8279e-12,-32319.3,24.4373], Tmin=(100,'K'), Tmax=(595.603,'K')), NASAPolynomial(coeffs=[8.29409,0.0334578,-1.64172e-05,3.21267e-09,-2.26961e-13,-33362.4,-8.11259], Tmin=(595.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]CCF(236)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {6,S} {11,S}
6  C u1 p0 c0 {2,S} {5,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-132.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,3025,407.5,1350,352.5,334,575,1197,1424,3202,238.226,1895.59,1895.59],'cm^-1')),
        HinderedRotor(inertia=(0.216241,'amu*angstrom^2'), symmetry=1, barrier=(8.70859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216242,'amu*angstrom^2'), symmetry=1, barrier=(8.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00885,'amu*angstrom^2'), symmetry=1, barrier=(40.628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0872,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70912,0.0563524,-7.72297e-05,7.12206e-08,-2.72522e-11,-15860.3,24.9519], Tmin=(100,'K'), Tmax=(780.439,'K')), NASAPolynomial(coeffs=[3.92183,0.0350026,-1.69582e-05,3.30268e-09,-2.31964e-13,-15900.9,16.7782], Tmin=(780.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]OF(1538)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u1 p0 c0 {3,S} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (63.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.363305,'amu*angstrom^2'), symmetry=1, barrier=(8.35309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164954,'amu*angstrom^2'), symmetry=1, barrier=(13.8793,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89237,0.0524723,-9.6079e-05,8.76882e-08,-3.04288e-11,7710.07,18.7902], Tmin=(100,'K'), Tmax=(846.972,'K')), NASAPolynomial(coeffs=[7.45195,0.0143542,-7.56376e-06,1.48098e-09,-1.02319e-13,7193.77,-4.59654], Tmin=(846.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = 'CH2F(46)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-42.5685,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(33.0141,'amu')),
        NonlinearRotor(inertia=([1.91548,16.2277,17.9803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([576.418,1180.5,1217.62,1485.55,3118.23,3268.88],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03338,-0.00262849,2.74227e-05,-3.89096e-08,1.85259e-11,-5119.82,5.20374], Tmin=(10,'K'), Tmax=(594.366,'K')), NASAPolynomial(coeffs=[2.59024,0.00857266,-4.60348e-06,1.22743e-09,-1.29255e-13,-4974.57,11.194], Tmin=(594.366,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-42.5685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[CH2]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C([CH]F)OF(1965)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5  C u1 p0 c0 {1,S} {4,S} {8,S}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (48.6069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.202817,'amu*angstrom^2'), symmetry=1, barrier=(4.66317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203698,'amu*angstrom^2'), symmetry=1, barrier=(4.68341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204037,'amu*angstrom^2'), symmetry=1, barrier=(4.69122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1511,0.0734523,-0.000136616,1.30618e-07,-4.7012e-11,5938.29,21.7074], Tmin=(100,'K'), Tmax=(856.021,'K')), NASAPolynomial(coeffs=[5.95303,0.027251,-1.4018e-05,2.70997e-09,-1.85897e-13,5986.81,4.37276], Tmin=(856.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.6069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-O2sCsH)(F1s)(H)) + radical(CJCO)"""),
)

species(
    label = 'F[CH]C([CH]CF)OF(1966)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {6,S} {12,S}
8  C u1 p0 c0 {2,S} {5,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-165.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,334,575,1197,1424,3202,264.596,264.635,264.635,264.645],'cm^-1')),
        HinderedRotor(inertia=(0.00240669,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240834,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266392,'amu*angstrom^2'), symmetry=1, barrier=(13.2385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680451,'amu*angstrom^2'), symmetry=1, barrier=(33.8202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.601172,0.0831422,-0.000112244,7.61195e-08,-1.71394e-11,-19788.1,29.0741], Tmin=(100,'K'), Tmax=(613.429,'K')), NASAPolynomial(coeffs=[10.2485,0.0320575,-1.62377e-05,3.20042e-09,-2.26262e-13,-21194.1,-14.5671], Tmin=(613.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(CCJCO) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH][C](CCF)OF(1967)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {8,S}
8  C u1 p0 c0 {2,S} {7,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-185.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,360,370,350,334,575,1197,1424,3202,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.175392,'amu*angstrom^2'), symmetry=1, barrier=(4.03261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17486,'amu*angstrom^2'), symmetry=1, barrier=(4.02039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174451,'amu*angstrom^2'), symmetry=1, barrier=(4.01098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15207,'amu*angstrom^2'), symmetry=1, barrier=(26.4884,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0792863,0.0990475,-0.000176188,1.64719e-07,-5.88768e-11,-22154.8,28.8266], Tmin=(100,'K'), Tmax=(845.62,'K')), NASAPolynomial(coeffs=[7.49915,0.0372781,-1.93079e-05,3.73961e-09,-2.57517e-13,-22456.1,-0.0872646], Tmin=(845.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(C2CsJO) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]CC([CH]F)OF(1968)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {5,S} {12,S}
8  C u1 p0 c0 {2,S} {6,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-173.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,262,406,528,622,1148,1246,1368,1480,3164,3240,248.015,248.016,248.017,2905.83],'cm^-1')),
        HinderedRotor(inertia=(0.203653,'amu*angstrom^2'), symmetry=1, barrier=(8.88949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709297,'amu*angstrom^2'), symmetry=1, barrier=(30.961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203653,'amu*angstrom^2'), symmetry=1, barrier=(8.88949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709303,'amu*angstrom^2'), symmetry=1, barrier=(30.961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0331956,0.0973807,-0.000164115,1.47233e-07,-5.14569e-11,-20707.1,29.7124], Tmin=(100,'K'), Tmax=(828.283,'K')), NASAPolynomial(coeffs=[9.43089,0.0342372,-1.76018e-05,3.42088e-09,-2.36991e-13,-21654.7,-10.178], Tmin=(828.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(Csj(Cs-O2sCsH)(F1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'C=CC([CH]F)OF(1969)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u1 p0 c0 {1,S} {4,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-48.4099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,180,516.561,516.575],'cm^-1')),
        HinderedRotor(inertia=(0.681275,'amu*angstrom^2'), symmetry=1, barrier=(15.6639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0827197,'amu*angstrom^2'), symmetry=1, barrier=(15.6638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681278,'amu*angstrom^2'), symmetry=1, barrier=(15.6639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18166,0.0636125,-6.65412e-05,3.61049e-08,-7.93911e-12,-5722.19,25.2555], Tmin=(100,'K'), Tmax=(1089.77,'K')), NASAPolynomial(coeffs=[12.1589,0.0233199,-1.10798e-05,2.17571e-09,-1.55384e-13,-8114.66,-28.6465], Tmin=(1089.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.4099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CsCsF1sH)"""),
)

species(
    label = 'O=C([CH]F)CCF(1970)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {4,S} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-456.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,375,552.5,462.5,1710,235,1215,1347,1486,3221,323.822,323.823,323.823],'cm^-1')),
        HinderedRotor(inertia=(0.0759378,'amu*angstrom^2'), symmetry=1, barrier=(5.65066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0759378,'amu*angstrom^2'), symmetry=1, barrier=(5.65066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533031,'amu*angstrom^2'), symmetry=1, barrier=(39.6637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50127,0.0598531,-7.20124e-05,5.41144e-08,-1.75545e-11,-54822.1,21.6635], Tmin=(100,'K'), Tmax=(736.27,'K')), NASAPolynomial(coeffs=[6.59367,0.0321836,-1.56343e-05,3.05945e-09,-2.16653e-13,-55571.9,-1.34463], Tmin=(736.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + radical(Csj(CO-CsO2d)(F1s)(H))"""),
)

species(
    label = 'F2(78)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C2H4(30)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[C]C(CCF)OF-2(1971)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
8  C u0 p1 c0 {2,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-193.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,617,898,1187,286.172,286.454,287.126],'cm^-1')),
        HinderedRotor(inertia=(0.613586,'amu*angstrom^2'), symmetry=1, barrier=(35.5993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187364,'amu*angstrom^2'), symmetry=1, barrier=(10.9691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189689,'amu*angstrom^2'), symmetry=1, barrier=(10.9714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188292,'amu*angstrom^2'), symmetry=1, barrier=(10.9602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723772,0.0770515,-9.50093e-05,6.3672e-08,-1.74209e-11,-23152.9,27.0959], Tmin=(100,'K'), Tmax=(883.392,'K')), NASAPolynomial(coeffs=[11.2126,0.0295601,-1.43718e-05,2.81977e-09,-2.00355e-13,-25006.1,-22.2063], Tmin=(883.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-CsCsOsH) + group(CsCsFHH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'FCC[C](CF)OF(1177)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
8  C u1 p0 c0 {4,S} {5,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-379.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,551,1088,1226,1380,1420,1481,3057,3119,360,370,350,274.504,274.642,274.717,274.778],'cm^-1')),
        HinderedRotor(inertia=(0.168046,'amu*angstrom^2'), symmetry=1, barrier=(8.98385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549649,'amu*angstrom^2'), symmetry=1, barrier=(29.4224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167843,'amu*angstrom^2'), symmetry=1, barrier=(8.98438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167875,'amu*angstrom^2'), symmetry=1, barrier=(8.97959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28302,0.089128,-0.000133353,1.13699e-07,-3.92418e-11,-45556.6,27.7725], Tmin=(100,'K'), Tmax=(783.751,'K')), NASAPolynomial(coeffs=[9.19806,0.0362285,-1.79471e-05,3.48679e-09,-2.43634e-13,-46726.7,-11.6152], Tmin=(783.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(C2CsJO)"""),
)

species(
    label = 'FC[CH]C(CF)OF(1972)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
8  C u1 p0 c0 {5,S} {7,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-360.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,307.404,308.291,308.817,308.967],'cm^-1')),
        HinderedRotor(inertia=(0.243366,'amu*angstrom^2'), symmetry=1, barrier=(16.4991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244174,'amu*angstrom^2'), symmetry=1, barrier=(16.5122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.371677,'amu*angstrom^2'), symmetry=1, barrier=(25.1447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00176764,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.746328,0.0745277,-7.75658e-05,4.27963e-08,-9.67727e-12,-43187.5,28.195], Tmin=(100,'K'), Tmax=(1055.14,'K')), NASAPolynomial(coeffs=[12.4947,0.0299902,-1.42512e-05,2.79274e-09,-1.99078e-13,-45666.8,-29.115], Tmin=(1055.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(CCJCO)"""),
)

species(
    label = 'F[CH]CC(CF)OF(1973)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-367.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414851,0.0851149,-0.000111772,8.16949e-08,-2.44755e-11,-44116.4,28.0353], Tmin=(100,'K'), Tmax=(808.797,'K')), NASAPolynomial(coeffs=[10.7892,0.0338069,-1.66154e-05,3.25955e-09,-2.30891e-13,-45794.6,-19.8134], Tmin=(808.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O]C(CCF)C(F)F(1974)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-685.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798395,0.0741444,-7.74236e-05,4.36043e-08,-1.01476e-11,-82345.2,25.9377], Tmin=(100,'K'), Tmax=(1022.76,'K')), NASAPolynomial(coeffs=[11.6719,0.0316183,-1.5054e-05,2.94985e-09,-2.10181e-13,-84569.4,-26.7654], Tmin=(1022.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-685.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(OF)C(F)F(1975)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-387.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653789,0.078285,-8.94129e-05,5.55734e-08,-1.42011e-11,-46472.2,27.6724], Tmin=(100,'K'), Tmax=(939.599,'K')), NASAPolynomial(coeffs=[11.5542,0.031879,-1.53269e-05,3.006e-09,-2.13995e-13,-48520.6,-24.2362], Tmin=(939.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(RCCJ)"""),
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
    E0 = (-227.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (331.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (100.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-108.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (247.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (144.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (119.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-137.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-142.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-115.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (198.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-92.1523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (67.3061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (85.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (103.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (143.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (125.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (135.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-5.56328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-152.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-3.54291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (68.1083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (119.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-101.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-122.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-198.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-169.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-37.6928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['HF(38)', 'O=C[CH]F(338)', 'CH2CHF(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(40.7622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', 'F[CH]C(CF)OF(1957)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(153.348,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'CC([CH]F)OF(1958)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.39871e-05,'m^3/(mol*s)'), n=3.10993, Ea=(26.0301,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['FCCC(F)[CH]OF(1959)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', '[CH]C(CCF)OF-2(1960)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(230)', 'FCC[CH]OF(514)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'F[C]C(CCF)OF(1961)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]F(126)', 'FC=CCCF(232)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.20421e+10,'m^3/(mol*s)'), n=-1.2407, Ea=(60.2313,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08524545916686703, var=37.982834433071, Tref=1000.0, N=20, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC=COF(1541)', 'CH2F-CH2(62)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(18.4505,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'FC=C(CCF)OF(1962)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(3.13517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2]CC([CH]F)OF(1963)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -34.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C([CH]F)CCF(1964)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(7.3568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]F(126)', 'F[CH][CH]CCF(236)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[CH][CH]OF(1538)', 'CH2F-CH2(62)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -3.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2F(46)', '[CH2]C([CH]F)OF(1965)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -9.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'F[CH]C([CH]CF)OF(1966)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'F[CH][C](CCF)OF(1967)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(1.35813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'F[CH]CC([CH]F)OF(1968)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.0689648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'C=CC([CH]F)OF(1969)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(51.7395,'m^3/(mol*s)'), n=1.07494, Ea=(226.832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07038556604361025, var=0.028808843118544165, Tref=1000.0, N=2, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'O=C([CH]F)CCF(1970)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(488.354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['F2(78)', 'O=C[CH]F(338)', 'C2H4(30)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(264.711,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CHF(40)', 'FCC[CH]OF(514)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'F[C]C(CCF)OF-2(1971)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(4.33194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['FCC[C](CF)OF(1177)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.14713e+11,'s^-1'), n=0.563333, Ea=(166.523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC[CH]C(CF)OF(1972)'],
    products = ['F[CH]C(CCF)OF(1953)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.73898e+09,'s^-1'), n=0.87, Ea=(139.955,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['F[CH]CC(CF)OF(1973)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.67145,'s^-1'), n=3.19634, Ea=(69.5908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_1H;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['[O]C(CCF)C(F)F(1974)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(98.7168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(CCF)OF(1953)'],
    products = ['[CH2]CC(OF)C(F)F(1975)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(230.561,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #644',
    isomers = [
        'F[CH]C(CCF)OF(1953)',
    ],
    reactants = [
        ('HF(38)', 'O=C[CH]F(338)', 'CH2CHF(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #644',
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

