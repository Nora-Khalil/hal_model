species(
    label = 'CC(F)(CC([O])F)OF(7520)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {13,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-585.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,391,562,707,872,1109,1210,1289,3137,183.9,184.225,184.367],'cm^-1')),
        HinderedRotor(inertia=(0.48799,'amu*angstrom^2'), symmetry=1, barrier=(11.7012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489209,'amu*angstrom^2'), symmetry=1, barrier=(11.7066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.487959,'amu*angstrom^2'), symmetry=1, barrier=(11.7007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994524,'amu*angstrom^2'), symmetry=1, barrier=(23.745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132153,0.0911912,-0.000109383,7.05464e-08,-1.86071e-11,-70234.5,28.2101], Tmin=(100,'K'), Tmax=(913.431,'K')), NASAPolynomial(coeffs=[12.9125,0.0352247,-1.74774e-05,3.46876e-09,-2.48349e-13,-72569.3,-32.2904], Tmin=(913.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-585.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH))"""),
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
    label = 'CH2CFO(78)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-261.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,611,648,830,1210,1753],'cm^-1')),
        HinderedRotor(inertia=(0.112588,'amu*angstrom^2'), symmetry=1, barrier=(22.3442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(3.302,'De'), polarizability=(3.442,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95484,0.00266306,5.93929e-05,-1.151e-07,6.74789e-11,-31406.5,8.78105], Tmin=(10,'K'), Tmax=(558.359,'K')), NASAPolynomial(coeffs=[3.00171,0.0196207,-1.33752e-05,4.27391e-09,-5.17274e-13,-31457.9,11.4099], Tmin=(558.359,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-261.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[CH2]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC(=O)F(499)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-453.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926],'cm^-1')),
        HinderedRotor(inertia=(0.163766,'amu*angstrom^2'), symmetry=1, barrier=(3.76529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3044.49,'J/mol'), sigma=(4.92747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.54 K, Pc=57.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95746,0.00861016,1.67316e-05,-2.39571e-08,8.75423e-12,-54563.9,8.39893], Tmin=(10,'K'), Tmax=(963.489,'K')), NASAPolynomial(coeffs=[3.63904,0.0176469,-9.3479e-06,2.39871e-09,-2.40789e-13,-54860.6,8.06499], Tmin=(963.489,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-453.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C[C]F(124)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (80.9091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1049.2],'cm^-1')),
        HinderedRotor(inertia=(0.0264144,'amu*angstrom^2'), symmetry=1, barrier=(3.87782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1971.36,'J/mol'), sigma=(5.118e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50526,0.00692148,1.56045e-05,-2.26115e-08,8.49098e-12,9752.45,9.13904], Tmin=(100,'K'), Tmax=(981.304,'K')), NASAPolynomial(coeffs=[5.2739,0.00958149,-3.54759e-06,6.48822e-10,-4.59773e-14,8930.14,-1.78146], Tmin=(981.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]C(F)COF(3105)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-292.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.939135,'amu*angstrom^2'), symmetry=1, barrier=(21.5926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942161,'amu*angstrom^2'), symmetry=1, barrier=(21.6621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.80935,0.012114,0.000141719,-3.48996e-07,2.41923e-10,-35222.1,11.9311], Tmin=(10,'K'), Tmax=(515.887,'K')), NASAPolynomial(coeffs=[6.45981,0.0300484,-2.23267e-05,7.60143e-09,-9.605e-13,-36007.7,-4.06507], Tmin=(515.887,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-292.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(F)COF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC(F)(OF)C([O])F(5845)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-546.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,391,562,707,872,1109,1210,1289,3137,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.820782,'amu*angstrom^2'), symmetry=1, barrier=(18.8714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825007,'amu*angstrom^2'), symmetry=1, barrier=(18.9685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826616,'amu*angstrom^2'), symmetry=1, barrier=(19.0055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.805921,0.0739227,-8.66232e-05,5.24944e-08,-1.28279e-11,-65631.8,23.9091], Tmin=(100,'K'), Tmax=(988.275,'K')), NASAPolynomial(coeffs=[12.8718,0.0250867,-1.25002e-05,2.49278e-09,-1.79246e-13,-68016.6,-34.1594], Tmin=(988.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-546.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(F)CC(F)OF(7595)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-530.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,391,562,707,872,1109,1210,1289,3137,277.162,277.171,277.413],'cm^-1')),
        HinderedRotor(inertia=(0.147657,'amu*angstrom^2'), symmetry=1, barrier=(8.06977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14824,'amu*angstrom^2'), symmetry=1, barrier=(8.07154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499397,'amu*angstrom^2'), symmetry=1, barrier=(27.1671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96038,0.0715217,-8.47872e-05,5.34194e-08,-1.3733e-11,-63700.7,25.107], Tmin=(100,'K'), Tmax=(936.56,'K')), NASAPolynomial(coeffs=[11.3774,0.0270313,-1.35314e-05,2.698e-09,-1.93756e-13,-65652,-24.4664], Tmin=(936.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'OF(169)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC(F)=CC([O])F(7596)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-396.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,361,769,965,1078,1132,1246,3247,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.358966,'amu*angstrom^2'), symmetry=1, barrier=(8.25333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359492,'amu*angstrom^2'), symmetry=1, barrier=(8.26543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65071,0.0546032,-5.01868e-05,2.53709e-08,-5.43485e-12,-47632.2,21.2248], Tmin=(100,'K'), Tmax=(1087.22,'K')), NASAPolynomial(coeffs=[8.99585,0.0275797,-1.29034e-05,2.50933e-09,-1.77969e-13,-49229.4,-14.8254], Tmin=(1087.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFHO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(O2sj(Cs-F1sCdH))"""),
)

species(
    label = 'C=C(F)CC([O])F(7597)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-383.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.220608,'amu*angstrom^2'), symmetry=1, barrier=(5.0722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220591,'amu*angstrom^2'), symmetry=1, barrier=(5.07183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51791,0.0553401,-5.08439e-05,2.50394e-08,-5.10156e-12,-46093.8,22.8637], Tmin=(100,'K'), Tmax=(1157.43,'K')), NASAPolynomial(coeffs=[10.4078,0.0246172,-1.10275e-05,2.10554e-09,-1.47913e-13,-48151.7,-21.3245], Tmin=(1157.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsHH) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC(F)(C[CH]F)OF(2940)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-422.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,334,575,1197,1424,3202,236.517,240.919,243.916],'cm^-1')),
        HinderedRotor(inertia=(0.288487,'amu*angstrom^2'), symmetry=1, barrier=(11.6103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281486,'amu*angstrom^2'), symmetry=1, barrier=(11.5625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413134,'amu*angstrom^2'), symmetry=1, barrier=(17.4972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678434,'amu*angstrom^2'), symmetry=1, barrier=(29.4221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390238,0.0860839,-0.000109963,7.72222e-08,-2.22792e-11,-50732.3,25.586], Tmin=(100,'K'), Tmax=(837.139,'K')), NASAPolynomial(coeffs=[11.2109,0.0343804,-1.73192e-05,3.4437e-09,-2.46189e-13,-52544,-24.6939], Tmin=(837.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    label = 'CC(F)(CC=O)OF(5730)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {4,D} {6,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-508.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.524087,'amu*angstrom^2'), symmetry=1, barrier=(12.0498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52367,'amu*angstrom^2'), symmetry=1, barrier=(12.0402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524472,'amu*angstrom^2'), symmetry=1, barrier=(12.0586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3157,'amu*angstrom^2'), symmetry=1, barrier=(30.2505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3506.73,'J/mol'), sigma=(5.8765,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=547.74 K, Pc=39.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674212,0.079931,-9.91672e-05,7.0864e-08,-2.12598e-11,-61094.5,24.5199], Tmin=(100,'K'), Tmax=(799.784,'K')), NASAPolynomial(coeffs=[9.26587,0.0369625,-1.8582e-05,3.69355e-09,-2.6403e-13,-62468.8,-15.0106], Tmin=(799.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-508.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CHFO(47)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(C)(F)OF(2912)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-198.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.423989,'amu*angstrom^2'), symmetry=1, barrier=(9.74833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423828,'amu*angstrom^2'), symmetry=1, barrier=(9.74464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978508,'amu*angstrom^2'), symmetry=1, barrier=(22.4978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73228,0.0203706,0.000170993,-5.4485e-07,4.85075e-10,-23912.6,11.0548], Tmin=(10,'K'), Tmax=(401.706,'K')), NASAPolynomial(coeffs=[5.86655,0.0359386,-2.4628e-05,7.97562e-09,-9.77322e-13,-24381.2,-0.993067], Tmin=(401.706,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-198.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""[CH2]C(C)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,25474.2,-0.444973], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC(F)(CC(=O)F)OF(7598)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-762.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,180,182.251,193.823],'cm^-1')),
        HinderedRotor(inertia=(0.413986,'amu*angstrom^2'), symmetry=1, barrier=(9.51835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682958,'amu*angstrom^2'), symmetry=1, barrier=(15.7026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42909,'amu*angstrom^2'), symmetry=1, barrier=(32.8576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416088,'amu*angstrom^2'), symmetry=1, barrier=(9.56668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.174986,0.0917843,-0.000128769,1.01397e-07,-3.28226e-11,-91625.8,27.2295], Tmin=(100,'K'), Tmax=(749.642,'K')), NASAPolynomial(coeffs=[10.4549,0.0369315,-1.90103e-05,3.78615e-09,-2.69783e-13,-93167.1,-19.4027], Tmin=(749.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-762.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(COCsFO)"""),
)

species(
    label = 'C[C](CC([O])F)OF(7599)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
7  C u0 p0 c0 {8,S} {12,S} {13,S} {14,S}
8  C u1 p0 c0 {3,S} {5,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-161.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,229.212,229.212,1871.33],'cm^-1')),
        HinderedRotor(inertia=(0.218196,'amu*angstrom^2'), symmetry=1, barrier=(8.13481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218196,'amu*angstrom^2'), symmetry=1, barrier=(8.13481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218196,'amu*angstrom^2'), symmetry=1, barrier=(8.13481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218196,'amu*angstrom^2'), symmetry=1, barrier=(8.13481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375351,0.086457,-0.000127176,1.06672e-07,-3.61436e-11,-19254.8,29.3413], Tmin=(100,'K'), Tmax=(797.727,'K')), NASAPolynomial(coeffs=[9.26096,0.0348054,-1.67085e-05,3.20073e-09,-2.21801e-13,-20446.6,-10.1032], Tmin=(797.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJO)"""),
)

species(
    label = 'CC(F)(C[CH][O])OF(7600)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {4,S} {6,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-179.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,222.881,222.881,222.882,222.882],'cm^-1')),
        HinderedRotor(inertia=(0.238665,'amu*angstrom^2'), symmetry=1, barrier=(8.41324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238665,'amu*angstrom^2'), symmetry=1, barrier=(8.41324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00339354,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05381,'amu*angstrom^2'), symmetry=1, barrier=(37.1481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.011533,0.0978438,-0.000159153,1.43403e-07,-5.0988e-11,-21507.8,27.9564], Tmin=(100,'K'), Tmax=(813.245,'K')), NASAPolynomial(coeffs=[8.48243,0.0393182,-2.01055e-05,3.92335e-09,-2.73395e-13,-22328,-7.73167], Tmin=(813.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CC([O])(F)CC([O])F(7601)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-485.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,334,437,422,648,895,1187,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,391,562,707,872,1109,1210,1289,3137,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.21357,'amu*angstrom^2'), symmetry=1, barrier=(4.9104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213636,'amu*angstrom^2'), symmetry=1, barrier=(4.91192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213417,'amu*angstrom^2'), symmetry=1, barrier=(4.90687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.720474,0.0777634,-9.5888e-05,6.76395e-08,-1.98228e-11,-58295.2,27.0645], Tmin=(100,'K'), Tmax=(822.041,'K')), NASAPolynomial(coeffs=[9.66918,0.0342198,-1.64335e-05,3.20319e-09,-2.26479e-13,-59766.5,-14.3542], Tmin=(822.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sCsCs)) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]F(160)',
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
    label = 'C[C](F)CC([O])F(7602)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-328.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,180,2934.77],'cm^-1')),
        HinderedRotor(inertia=(0.353301,'amu*angstrom^2'), symmetry=1, barrier=(8.12309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353621,'amu*angstrom^2'), symmetry=1, barrier=(8.13045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353826,'amu*angstrom^2'), symmetry=1, barrier=(8.13516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13412,0.0689583,-9.63882e-05,8.27128e-08,-2.91638e-11,-39471.7,25.2849], Tmin=(100,'K'), Tmax=(795.053,'K')), NASAPolynomial(coeffs=[6.43023,0.034427,-1.63609e-05,3.1328e-09,-2.17451e-13,-40064.6,2.5163], Tmin=(795.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsHH)(Cs-HHH)(F1s))"""),
)

species(
    label = '[CH2]C([O])F(527)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-49.3749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00206763,'amu*angstrom^2'), symmetry=1, barrier=(12.9133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3215.69,'J/mol'), sigma=(5.42929,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.28 K, Pc=45.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95775,0.0246236,-2.50869e-05,1.62783e-08,-4.59836e-12,-5902.4,14.7453], Tmin=(100,'K'), Tmax=(836.479,'K')), NASAPolynomial(coeffs=[5.16834,0.0140521,-6.12857e-06,1.16774e-09,-8.19841e-14,-6272.21,4.47531], Tmin=(836.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.3749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = 'C[C](F)OF(5141)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-167.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,180],'cm^-1')),
        HinderedRotor(inertia=(0.288089,'amu*angstrom^2'), symmetry=1, barrier=(6.62374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289015,'amu*angstrom^2'), symmetry=1, barrier=(6.64503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75318,0.033278,-2.07715e-05,-3.05874e-08,4.32155e-11,-20144.3,14.9774], Tmin=(100,'K'), Tmax=(468.497,'K')), NASAPolynomial(coeffs=[5.68823,0.0189402,-9.19318e-06,1.78427e-09,-1.2484e-13,-20537,1.78713], Tmin=(468.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s)"""),
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
    label = '[O]C(F)C[C](F)OF(7603)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
8  C u1 p0 c0 {2,S} {4,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-335.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,395,473,707,1436,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.267687,'amu*angstrom^2'), symmetry=1, barrier=(6.15466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267769,'amu*angstrom^2'), symmetry=1, barrier=(6.15653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.509148,'amu*angstrom^2'), symmetry=1, barrier=(11.7063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550115,0.0835699,-0.000138851,1.21591e-07,-4.17446e-11,-40290,27.3645], Tmin=(100,'K'), Tmax=(816.363,'K')), NASAPolynomial(coeffs=[9.85756,0.0267573,-1.38688e-05,2.70853e-09,-1.88397e-13,-41436.2,-13.3626], Tmin=(816.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][CH]F(1657)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53592,0.00799825,-5.22271e-07,-4.56612e-09,2.08752e-12,-2348.33,9.31392], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97187,0.0038822,-1.62979e-06,3.36437e-10,-2.54188e-14,-3160.18,-3.94923], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'CC(F)([CH]C([O])F)OF(7604)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-385.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,391,562,707,872,1109,1210,1289,3137,3025,407.5,1350,352.5,220.526,220.568,220.607],'cm^-1')),
        HinderedRotor(inertia=(0.110426,'amu*angstrom^2'), symmetry=1, barrier=(3.81516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00346527,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305228,'amu*angstrom^2'), symmetry=1, barrier=(10.5377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685978,'amu*angstrom^2'), symmetry=1, barrier=(23.6746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292987,0.0864335,-0.000102615,6.38004e-08,-1.6061e-11,-46196.5,30.1253], Tmin=(100,'K'), Tmax=(959.016,'K')), NASAPolynomial(coeffs=[13.694,0.0305376,-1.51866e-05,3.02301e-09,-2.17022e-13,-48766.8,-33.9659], Tmin=(959.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(F)(CC([O])F)OF(7605)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {7,S} {12,S}
9  C u1 p0 c0 {6,S} {13,S} {14,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-374.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,316,385,515,654,689,1295,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.328715,'amu*angstrom^2'), symmetry=1, barrier=(7.5578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329812,'amu*angstrom^2'), symmetry=1, barrier=(7.58302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.551,'amu*angstrom^2'), symmetry=1, barrier=(35.6605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329362,'amu*angstrom^2'), symmetry=1, barrier=(7.57268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399342,0.104979,-0.000164347,1.3828e-07,-4.65181e-11,-44858.5,30.3481], Tmin=(100,'K'), Tmax=(774.078,'K')), NASAPolynomial(coeffs=[12.2883,0.0343434,-1.76404e-05,3.46487e-09,-2.4315e-13,-46670.7,-26.6318], Tmin=(774.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'CC(F)(C[C]([O])F)OF(7606)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {2,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-399.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,300.553,300.672,300.857,301.248],'cm^-1')),
        HinderedRotor(inertia=(0.35274,'amu*angstrom^2'), symmetry=1, barrier=(22.5843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35149,'amu*angstrom^2'), symmetry=1, barrier=(22.5845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352154,'amu*angstrom^2'), symmetry=1, barrier=(22.5554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13381,'amu*angstrom^2'), symmetry=1, barrier=(8.59314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320164,0.0870287,-0.000100396,6.0517e-08,-1.49053e-11,-47865.2,28.4835], Tmin=(100,'K'), Tmax=(973.294,'K')), NASAPolynomial(coeffs=[13.4716,0.0329805,-1.71007e-05,3.46444e-09,-2.51095e-13,-50425.3,-34.6089], Tmin=(973.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsHH)(F1s)(O2s-H))"""),
)

species(
    label = 'F2(77)',
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
    label = 'CC(=O)CC([O])F(7607)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-417.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180,2826.82],'cm^-1')),
        HinderedRotor(inertia=(0.327899,'amu*angstrom^2'), symmetry=1, barrier=(7.53903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327815,'amu*angstrom^2'), symmetry=1, barrier=(7.53712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327797,'amu*angstrom^2'), symmetry=1, barrier=(7.5367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70116,0.0569756,-5.23662e-05,1.77848e-08,5.65626e-12,-50095.7,23.0672], Tmin=(100,'K'), Tmax=(564.904,'K')), NASAPolynomial(coeffs=[6.13202,0.0345049,-1.63411e-05,3.17082e-09,-2.23813e-13,-50738.4,2.96387], Tmin=(564.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCFHO) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'CC(=CC([O])F)OF(7608)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-244.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,361,769,965,1078,1132,1246,3247,350,440,435,1725,3010,987.5,1337.5,450,1655,338.577,338.58,338.582],'cm^-1')),
        HinderedRotor(inertia=(0.100036,'amu*angstrom^2'), symmetry=1, barrier=(8.13768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100033,'amu*angstrom^2'), symmetry=1, barrier=(8.13761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100037,'amu*angstrom^2'), symmetry=1, barrier=(8.1378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01313,0.0708091,-8.25347e-05,5.39535e-08,-1.46765e-11,-29360.4,25.9505], Tmin=(100,'K'), Tmax=(880.856,'K')), NASAPolynomial(coeffs=[9.72416,0.0312558,-1.51863e-05,2.98643e-09,-2.12666e-13,-30895.1,-14.9709], Tmin=(880.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCFHO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(O2sj(Cs-F1sCdH))"""),
)

species(
    label = 'CC(F)(C=C[O])OF(5878)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {4,S} {7,D} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-366.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3161,'amu*angstrom^2'), symmetry=1, barrier=(30.2597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31601,'amu*angstrom^2'), symmetry=1, barrier=(30.2576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31566,'amu*angstrom^2'), symmetry=1, barrier=(30.2497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0953729,0.0786668,-8.16995e-05,4.11183e-08,-8.0415e-12,-43968.8,25.6234], Tmin=(100,'K'), Tmax=(1251.57,'K')), NASAPolynomial(coeffs=[19.5953,0.0163451,-7.00674e-06,1.33189e-09,-9.41503e-14,-48849.9,-72.8279], Tmin=(1251.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C=C(CC([O])F)OF(7609)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-232.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,350,440,435,1725,2950,3100,1380,975,1025,1650,295.115,295.283,295.316],'cm^-1')),
        HinderedRotor(inertia=(0.00193273,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106372,'amu*angstrom^2'), symmetry=1, barrier=(6.58107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106316,'amu*angstrom^2'), symmetry=1, barrier=(6.58149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966928,0.0705272,-7.96566e-05,4.90811e-08,-1.24211e-11,-27825.7,27.2788], Tmin=(100,'K'), Tmax=(949.165,'K')), NASAPolynomial(coeffs=[10.8957,0.0286837,-1.35275e-05,2.63246e-09,-1.86649e-13,-29710.4,-20.1034], Tmin=(949.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'CH2CHO(35)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,530.062,530.609,975.335,1640.89,1642.3],'cm^-1')),
        HinderedRotor(inertia=(0.00115918,'amu*angstrom^2'), symmetry=1, barrier=(2.21299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CC(F)(C[C](O)F)OF(7610)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {9,S} {15,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {2,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-625.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.108626,0.0906128,-9.87771e-05,5.51125e-08,-1.25061e-11,-75107.1,28.3345], Tmin=(100,'K'), Tmax=(1053.23,'K')), NASAPolynomial(coeffs=[15.1852,0.0333538,-1.72286e-05,3.494e-09,-2.53592e-13,-78282.9,-45.1831], Tmin=(1053.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(Csj(Cs-CsHH)(F1s)(O2s-H))"""),
)

species(
    label = 'CC(F)([CH]C(O)F)OF(7611)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {8,S} {15,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-611.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0794747,0.0900195,-0.000100916,5.81965e-08,-1.3546e-11,-73438.3,29.985], Tmin=(100,'K'), Tmax=(1034.21,'K')), NASAPolynomial(coeffs=[15.3235,0.0310608,-1.54039e-05,3.07423e-09,-2.2135e-13,-76591.4,-44.0712], Tmin=(1034.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(F)(CC(O)F)OF(7612)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,S} {15,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {12,S}
9  C u1 p0 c0 {7,S} {13,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-600.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.475787,0.106671,-0.00015466,1.20232e-07,-3.7687e-11,-72106,29.7339], Tmin=(100,'K'), Tmax=(778.521,'K')), NASAPolynomial(coeffs=[13.2722,0.0360345,-1.85618e-05,3.68784e-09,-2.62072e-13,-74246.6,-33.1502], Tmin=(778.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-600.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'CC([O])(F)CC(F)OF(7613)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {12,S}
9  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-577.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,286,334,437,422,648,895,1187,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,205.027,205.68,206.046],'cm^-1')),
        HinderedRotor(inertia=(0.00746213,'amu*angstrom^2'), symmetry=1, barrier=(11.4762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387396,'amu*angstrom^2'), symmetry=1, barrier=(11.4498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.749116,'amu*angstrom^2'), symmetry=1, barrier=(21.9256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378726,'amu*angstrom^2'), symmetry=1, barrier=(11.4711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.114224,0.0925678,-0.000116162,7.97636e-08,-2.2502e-11,-69347.6,28.5157], Tmin=(100,'K'), Tmax=(855.132,'K')), NASAPolynomial(coeffs=[12.0043,0.0369494,-1.85997e-05,3.7026e-09,-2.65067e-13,-71381,-26.9862], Tmin=(855.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-577.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sCsCs))"""),
)

species(
    label = 'C[C](CC(F)OF)OF(7614)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {3,S} {9,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
8  C u0 p0 c0 {9,S} {13,S} {14,S} {15,S}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-253.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([504,610,1067,1155,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0367733,0.0983257,-0.000133299,9.25877e-08,-2.25601e-11,-30315.4,30.1322], Tmin=(100,'K'), Tmax=(623.819,'K')), NASAPolynomial(coeffs=[11.3977,0.0379243,-1.91223e-05,3.76284e-09,-2.65861e-13,-31993.3,-21.6511], Tmin=(623.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + radical(C2CsJO)"""),
)

species(
    label = 'CC(F)(C[CH]OF)OF(7615)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {2,S} {6,S}
5  O u0 p2 c0 {3,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {5,S} {7,S} {15,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-267.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([504,610,1067,1155,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.427153,0.10352,-0.000134928,9.07015e-08,-2.43774e-11,-32044.1,29.5176], Tmin=(100,'K'), Tmax=(906.002,'K')), NASAPolynomial(coeffs=[15.8412,0.0316925,-1.60035e-05,3.18941e-09,-2.28563e-13,-34991.8,-47.3611], Tmin=(906.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2sCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJO)"""),
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
    E0 = (-363.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (174.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (231.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (97.4824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-51.2842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-67.4175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (29.0972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-206.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-340.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-305.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (120.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (101.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-203.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-17.4036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-8.27983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (9.1331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-9.60213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (35.5362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (46.4852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (21.6782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-141.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-82.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-144.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-101.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-166.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-250.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-297.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-344.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-174.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (27.9741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (17.1296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(F)(CC([O])F)OF(7520)'],
    products = ['HF(38)', 'CH2CFO(78)', 'CC(=O)F(499)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(12.4956,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[C]F(124)', '[O]C(F)COF(3105)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(178.033,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', 'CC(F)(OF)C([O])F(5845)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(149.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(S)(25)', '[O]C(F)CC(F)OF(7595)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(486642,'m^3/(mol*s)'), n=0.225671, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(169)', 'CC(F)=CC([O])F(7596)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OF(169)', 'C=C(F)CC([O])F(7597)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1250,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', 'CC(F)(C[CH]F)OF(2940)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'CC(F)(CC=O)OF(5730)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(20.2035,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CHFO(47)', '[CH2]C(C)(F)OF(2912)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.43214e-05,'m^3/(mol*s)'), n=3.00879, Ea=(44.0701,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Sp-5R!H-4R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Sp-5R!H-4R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'CC(F)(CC(=O)F)OF(7598)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(36.3925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C[C](CC([O])F)OF(7599)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'CC(F)(C[CH][O])OF(7600)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'CC([O])(F)CC([O])F(7601)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(160)', 'C[C](F)CC([O])F(7602)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])F(527)', 'C[C](F)OF(5141)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH3(19)', '[O]C(F)C[C](F)OF(7603)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][CH]F(1657)', '[CH2]C(C)(F)OF(2912)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'CC(F)([CH]C([O])F)OF(7604)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[CH2]C(F)(CC([O])F)OF(7605)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'CC(F)(C[C]([O])F)OF(7606)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F2(77)', 'CC(=O)CC([O])F(7607)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(75.9776,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'CC(=CC([O])F)OF(7608)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(234.741,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'CC(F)(C=C[O])OF(5878)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(294.996,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'C=C(CC([O])F)OF(7609)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(202.493,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(F)(CC([O])F)OF(7520)'],
    products = ['F2(77)', 'CH2CHO(35)', 'CC(=O)F(499)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(209.215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC(F)(CC([O])F)OF(7520)'],
    products = ['CC(F)(C[C](O)F)OF(7610)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC(F)([CH]C(O)F)OF(7611)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC(F)(CC([O])F)OF(7520)'],
    products = ['[CH2]C(F)(CC(O)F)OF(7612)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC([O])(F)CC(F)OF(7613)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(194.503,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[C](CC(F)OF)OF(7614)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(72.2579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC(F)(C[CH]OF)OF(7615)'],
    products = ['CC(F)(CC([O])F)OF(7520)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(75.9345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #2858',
    isomers = [
        'CC(F)(CC([O])F)OF(7520)',
    ],
    reactants = [
        ('HF(38)', 'CH2CFO(78)', 'CC(=O)F(499)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2858',
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

