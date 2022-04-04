species(
    label = 'OC(F)(CF)C(F)C(F)F(3507)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {15,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1322.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,235,523,627,1123,1142,1372,1406,3097,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0047,'amu*angstrom^2'), symmetry=1, barrier=(23.1,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00439,'amu*angstrom^2'), symmetry=1, barrier=(23.0929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00477,'amu*angstrom^2'), symmetry=1, barrier=(23.1017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0042,'amu*angstrom^2'), symmetry=1, barrier=(23.0885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.284395,0.0989688,-0.00011891,7.32524e-08,-1.80953e-11,-158951,29.5223], Tmin=(100,'K'), Tmax=(980.767,'K')), NASAPolynomial(coeffs=[16.2503,0.0315347,-1.57776e-05,3.15088e-09,-2.26723e-13,-162195,-49.9272], Tmin=(980.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1322.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F))"""),
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
    label = 'CHFCHF[Z](59)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=C(F)CF(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OCF(781)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-436.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,548,1085,1183,1302,1466,1520,3060,3119],'cm^-1')),
        HinderedRotor(inertia=(0.193143,'amu*angstrom^2'), symmetry=1, barrier=(16.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0324,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03993,-0.00338798,6.00172e-05,-9.49101e-08,4.84175e-11,-52506.5,7.78779], Tmin=(10,'K'), Tmax=(607.141,'K')), NASAPolynomial(coeffs=[1.98113,0.0176998,-1.06708e-05,3.11934e-09,-3.52169e-13,-52395.2,15.551], Tmin=(607.141,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-436.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C(F)C(F)F(793)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
7 C u0 p1 c0 {4,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-532.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([353,444,1253,3145,235,523,627,1123,1142,1372,1406,3097,617,898,1187,353.605,356.141,356.23,358.402],'cm^-1')),
        HinderedRotor(inertia=(0.082132,'amu*angstrom^2'), symmetry=1, barrier=(7.3375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109694,'amu*angstrom^2'), symmetry=1, barrier=(34.096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78352,0.0528863,-6.48827e-05,4.24744e-08,-1.13901e-11,-63946.6,20.1219], Tmin=(100,'K'), Tmax=(897.526,'K')), NASAPolynomial(coeffs=[9.0773,0.0203829,-1.05657e-05,2.13221e-09,-1.53958e-13,-65256,-14.2784], Tmin=(897.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-532.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'F[C]CF(126)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p1 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-105.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,734,1109,1255,1358,2983,3011,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.0460479,'amu*angstrom^2'), symmetry=1, barrier=(35.5232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4418.31,'J/mol'), sigma=(4.687e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8659,0.0107238,1.79401e-05,-3.81644e-08,1.95617e-11,-12729.5,8.48308], Tmin=(10,'K'), Tmax=(672.698,'K')), NASAPolynomial(coeffs=[3.4045,0.0188633,-1.22417e-05,3.67097e-09,-4.17395e-13,-12789.5,9.6187], Tmin=(672.698,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-105.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)C(F)F(778)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {9,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-880.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,180],'cm^-1')),
        HinderedRotor(inertia=(1.07194,'amu*angstrom^2'), symmetry=1, barrier=(24.6459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0722,'amu*angstrom^2'), symmetry=1, barrier=(24.652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84462,0.00955472,0.000130159,-2.95984e-07,1.92886e-10,-105845,11.1136], Tmin=(10,'K'), Tmax=(538.658,'K')), NASAPolynomial(coeffs=[5.31308,0.0316785,-2.34231e-05,7.92737e-09,-9.9477e-13,-106482,0.489879], Tmin=(538.658,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-880.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2F2(41)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-461.906,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(52.0125,'amu')),
        NonlinearRotor(inertia=([10.1284,47.6239,54.4207],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([531.881,1130.54,1136.09,1194.11,1273.86,1488.46,1538.71,3014.54,3074.22],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0234,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06909,-0.00454455,4.03365e-05,-4.87809e-08,1.9023e-11,-55555.1,6.40435], Tmin=(10,'K'), Tmax=(788.976,'K')), NASAPolynomial(coeffs=[1.30655,0.0154345,-9.00408e-06,2.50681e-09,-2.68973e-13,-55305.1,17.899], Tmin=(788.976,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-461.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)([C]F)CF(4240)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7  C u0 p1 c0 {3,S} {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-533.352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.26549,'amu*angstrom^2'), symmetry=1, barrier=(29.096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26493,'amu*angstrom^2'), symmetry=1, barrier=(29.0832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26483,'amu*angstrom^2'), symmetry=1, barrier=(29.081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33255,0.0619959,-7.69453e-05,4.91174e-08,-1.25652e-11,-64054.2,21.2564], Tmin=(100,'K'), Tmax=(948.856,'K')), NASAPolynomial(coeffs=[11.4028,0.019543,-9.83278e-06,1.96354e-09,-1.41131e-13,-65965.2,-26.7981], Tmin=(948.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-533.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'F[C]C(F)F(125)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p1 c0 {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-295.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([257,409,1143,1361,2944,617,898,1187,180,1242.56,1243.85],'cm^-1')),
        HinderedRotor(inertia=(0.0989364,'amu*angstrom^2'), symmetry=1, barrier=(2.27474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89361,0.016532,-3.68777e-06,-7.04441e-09,3.72589e-12,-35574,10.6199], Tmin=(10,'K'), Tmax=(1003.12,'K')), NASAPolynomial(coeffs=[6.83938,0.0116737,-6.72299e-06,1.81811e-09,-1.88858e-13,-36511.5,-5.32834], Tmin=(1003.12,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-295.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)CF(773)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {9,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-662.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,180],'cm^-1')),
        HinderedRotor(inertia=(1.27049,'amu*angstrom^2'), symmetry=1, barrier=(29.211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27447,'amu*angstrom^2'), symmetry=1, barrier=(29.3025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0494,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89252,0.00611482,0.000108385,-2.11378e-07,1.20595e-10,-79630.2,9.64641], Tmin=(10,'K'), Tmax=(599.394,'K')), NASAPolynomial(coeffs=[3.69664,0.03329,-2.43572e-05,8.26298e-09,-1.04496e-12,-80071.4,6.61484], Tmin=(599.394,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-662.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(CF)C(F)F(780)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1108.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,235,523,627,1123,1142,1372,1406,3097,528,1116,1182,1331,1402,1494,3075,3110,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3187,'amu*angstrom^2'), symmetry=1, barrier=(30.3195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31813,'amu*angstrom^2'), symmetry=1, barrier=(30.3065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31967,'amu*angstrom^2'), symmetry=1, barrier=(30.3418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73978,0.015966,0.000200994,-4.66073e-07,3.05235e-10,-133256,13.0707], Tmin=(10,'K'), Tmax=(542.837,'K')), NASAPolynomial(coeffs=[7.13315,0.0462887,-3.56795e-05,1.23478e-08,-1.56876e-12,-134440,-8.73559], Tmin=(542.837,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1108.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(CF)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'H2(8)',
    structure = adjacencyList("""1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-8.60349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3765.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (2.01594,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(496.376,'J/mol'), sigma=(2.8327,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212708,-2.78619e-07,3.40262e-10,-7.7602e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.09,'K')), NASAPolynomial(coeffs=[2.78813,0.000587684,1.5899e-07,-5.52698e-11,4.34281e-15,-596.125,0.112931], Tmin=(1959.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(F)([C]F)C(F)C(F)F(5080)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
10 C u0 p1 c0 {5,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-971.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,417,511,1155,1315,1456,3119,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,617,898,1187,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12597,'amu*angstrom^2'), symmetry=1, barrier=(25.8884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12578,'amu*angstrom^2'), symmetry=1, barrier=(25.8838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12581,'amu*angstrom^2'), symmetry=1, barrier=(25.8846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12663,'amu*angstrom^2'), symmetry=1, barrier=(25.9035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254009,0.10029,-0.000142745,1.03498e-07,-2.97835e-11,-116641,30.3542], Tmin=(100,'K'), Tmax=(850.398,'K')), NASAPolynomial(coeffs=[15.2157,0.0275272,-1.44035e-05,2.88719e-09,-2.06661e-13,-119272,-41.7715], Tmin=(850.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-971.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)C(F)C(F)F(384)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1100.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,417,511,1155,1315,1456,3119,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.38268,'amu*angstrom^2'), symmetry=1, barrier=(31.7905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05499,'amu*angstrom^2'), symmetry=1, barrier=(24.2562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05433,'amu*angstrom^2'), symmetry=1, barrier=(24.2412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3403.26,'J/mol'), sigma=(5.47543,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=531.58 K, Pc=47.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51381,0.0450837,9.32577e-06,-7.7418e-08,4.98377e-11,-132374,12.4997], Tmin=(10,'K'), Tmax=(669.25,'K')), NASAPolynomial(coeffs=[8.46531,0.0381082,-2.5736e-05,8.00866e-09,-9.38069e-13,-133543,-13.1843], Tmin=(669.25,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1100.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)C(F)C(F)F(2802)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1332.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.921855,'amu*angstrom^2'), symmetry=1, barrier=(21.1953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.921854,'amu*angstrom^2'), symmetry=1, barrier=(21.1952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92307,'amu*angstrom^2'), symmetry=1, barrier=(21.2232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.047,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3397,'J/mol'), sigma=(5.67045,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.60 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.49088,0.0409947,0.000125294,-4.391e-07,3.64164e-10,-160305,13.3177], Tmin=(10,'K'), Tmax=(459.778,'K')), NASAPolynomial(coeffs=[8.77559,0.0420379,-3.15065e-05,1.06796e-08,-1.33929e-12,-161288,-13.4763], Tmin=(459.778,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1332.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'OC(F)(CF)CF(4239)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-888.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,442,614,1089,1143,1126,1238,1254,1408,1339,1465,1463,1525,3023,3127,3088,3132,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.34942,'amu*angstrom^2'), symmetry=1, barrier=(31.0259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35007,'amu*angstrom^2'), symmetry=1, barrier=(31.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35016,'amu*angstrom^2'), symmetry=1, barrier=(31.0429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69082,0.0293284,6.18299e-05,-1.49179e-07,8.60911e-11,-106898,11.8627], Tmin=(10,'K'), Tmax=(636.626,'K')), NASAPolynomial(coeffs=[5.90292,0.0413771,-2.76957e-05,8.59329e-09,-1.00627e-12,-107705,-1.94025], Tmin=(636.626,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-888.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(CF)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(CF)C(F)[C]F(5081)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
9  C u0 p1 c0 {4,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-741.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,353,444,1253,3145,528,1116,1182,1331,1402,1494,3075,3110,617,898,1187,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390344,0.0851958,-0.00010661,6.94189e-08,-1.82638e-11,-89051.5,27.6123], Tmin=(100,'K'), Tmax=(919.785,'K')), NASAPolynomial(coeffs=[13.3056,0.0290292,-1.50126e-05,3.02823e-09,-2.18661e-13,-91427.3,-33.6165], Tmin=(919.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-741.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(C(F)F)C(F)(F)CF(5082)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
10 C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1340.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506071,0.106635,-0.000150697,1.12505e-07,-3.36831e-11,-161046,30.9398], Tmin=(100,'K'), Tmax=(815.608,'K')), NASAPolynomial(coeffs=[14.2922,0.0340598,-1.7222e-05,3.40426e-09,-2.41612e-13,-163460,-37.437], Tmin=(815.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1340.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + group(CsCsFHH)"""),
)

species(
    label = 'OCC(F)(F)C(F)C(F)F(5083)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {15,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1331.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.601759,0.108942,-0.000160928,1.26624e-07,-3.98983e-11,-159950,31.0481], Tmin=(100,'K'), Tmax=(776.814,'K')), NASAPolynomial(coeffs=[13.8383,0.0345842,-1.7343e-05,3.39537e-09,-2.38822e-13,-162194,-34.9696], Tmin=(776.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1331.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-CsOsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787932,2.90876e-06,-1.47517e-09,2.12837e-13,-30281.6,-0.311363], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81029e-07,1.07243e-10,-5.79388e-15,-29958.6,5.91041], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FCC(F)=C(F)C(F)F(4565)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {7,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-963.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,171,205,598,1104,1143,1317,1411,3153,260,386,409,525,515,635,761,893,1354,1482,180,1267.42],'cm^-1')),
        HinderedRotor(inertia=(0.245007,'amu*angstrom^2'), symmetry=1, barrier=(5.63319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245119,'amu*angstrom^2'), symmetry=1, barrier=(5.63578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26777,0.0766109,-0.000216489,4.56294e-07,-3.86658e-10,-115901,14.2371], Tmin=(10,'K'), Tmax=(349.172,'K')), NASAPolynomial(coeffs=[5.08067,0.0446257,-3.08965e-05,9.94161e-09,-1.20474e-12,-115959,8.37763], Tmin=(349.172,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-963.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FCC(F)DC(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OF(482)',
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
    label = 'FCC(F)=CC(F)F(4123)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-807.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.368662,'amu*angstrom^2'), symmetry=1, barrier=(8.47627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369055,'amu*angstrom^2'), symmetry=1, barrier=(8.4853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19176,0.0618505,-8.24397e-05,6.65356e-08,-2.26698e-11,-97140.8,13.6216], Tmin=(10,'K'), Tmax=(759.39,'K')), NASAPolynomial(coeffs=[7.8049,0.0339801,-2.0334e-05,5.82046e-09,-6.42929e-13,-97738.5,-6.68639], Tmin=(759.39,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-807.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FCC(F)DCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C(F)C(F)C(F)F(3694)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-975.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,235,523,627,1123,1142,1372,1406,3097,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.421415,'amu*angstrom^2'), symmetry=1, barrier=(9.68915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423013,'amu*angstrom^2'), symmetry=1, barrier=(9.7259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.35802,0.0566747,-1.34892e-05,-1.05504e-07,1.09249e-10,-117268,13.9595], Tmin=(10,'K'), Tmax=(473.071,'K')), NASAPolynomial(coeffs=[7.81826,0.0400714,-2.77779e-05,8.95749e-09,-1.08698e-12,-117927,-6.71663], Tmin=(473.071,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-975.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FCDC(F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(F)C(F)C(F)F(4566)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-807.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,235,523,627,1123,1142,1372,1406,3097,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.278975,'amu*angstrom^2'), symmetry=1, barrier=(6.41417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639653,'amu*angstrom^2'), symmetry=1, barrier=(14.7069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43175,0.051312,-4.0944e-05,1.46939e-08,-1.48123e-12,-97087.7,13.2318], Tmin=(10,'K'), Tmax=(988.05,'K')), NASAPolynomial(coeffs=[10.2965,0.0295136,-1.69492e-05,4.64271e-09,-4.91304e-13,-98736.7,-21.2847], Tmin=(988.05,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-807.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""CDC(F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)F(368)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-687.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,519,593,830,1115,1166,1389,1413,3114],'cm^-1')),
        HinderedRotor(inertia=(0.761553,'amu*angstrom^2'), symmetry=1, barrier=(17.5096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0228,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93609,0.00366784,6.01076e-05,-1.20272e-07,7.01069e-11,-82625.1,8.12343], Tmin=(10,'K'), Tmax=(593.851,'K')), NASAPolynomial(coeffs=[4.25321,0.0168559,-1.19109e-05,4.02989e-09,-5.14929e-13,-82933,4.48375], Tmin=(593.851,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-687.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C(F)CF(486)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u0 p0 c0 {3,S} {5,D} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-541.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,323,467,575,827,1418,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.00655221,'amu*angstrom^2'), symmetry=1, barrier=(8.81279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75121,0.0252454,2.08009e-07,-1.91229e-08,9.65324e-12,-65113.8,11.461], Tmin=(10,'K'), Tmax=(869.181,'K')), NASAPolynomial(coeffs=[5.83446,0.0246751,-1.43688e-05,3.99307e-09,-4.28487e-13,-65816.6,-0.256776], Tmin=(869.181,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-541.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FCDC(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O[C](CF)C(F)C(F)F(5084)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {14,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-931.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,551,1088,1226,1380,1420,1481,3057,3119,360,370,350,192.302,192.303,192.304],'cm^-1')),
        HinderedRotor(inertia=(0.39354,'amu*angstrom^2'), symmetry=1, barrier=(10.3273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984551,'amu*angstrom^2'), symmetry=1, barrier=(25.8369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393536,'amu*angstrom^2'), symmetry=1, barrier=(10.3273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(10.3273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.330955,0.102983,-0.00016149,1.34936e-07,-4.47489e-11,-111911,29.9612], Tmin=(100,'K'), Tmax=(798.337,'K')), NASAPolynomial(coeffs=[12.4765,0.0321964,-1.60563e-05,3.10751e-09,-2.15905e-13,-113745,-27.6224], Tmin=(798.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-931.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(C2CsJOH)"""),
)

species(
    label = 'OC(F)([CH]C(F)F)CF(5085)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-945.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,528,1116,1182,1331,1402,1494,3075,3110,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,245.711,247.641,248.963],'cm^-1')),
        HinderedRotor(inertia=(0.522247,'amu*angstrom^2'), symmetry=1, barrier=(22.8052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297208,'amu*angstrom^2'), symmetry=1, barrier=(13.0003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522012,'amu*angstrom^2'), symmetry=1, barrier=(22.8172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202976,'amu*angstrom^2'), symmetry=1, barrier=(8.8944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.297736,0.0851548,-9.73405e-05,5.72871e-08,-1.3587e-11,-113554,29.1903], Tmin=(100,'K'), Tmax=(1016.94,'K')), NASAPolynomial(coeffs=[14.6105,0.0288583,-1.43038e-05,2.85232e-09,-2.05244e-13,-116465,-40.1011], Tmin=(1016.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-945.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + radical(CCJCO)"""),
)

species(
    label = 'OC(F)(CF)C(F)[CH]F(5086)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-914.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,259,529,569,1128,1321,1390,3140,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.974561,'amu*angstrom^2'), symmetry=1, barrier=(22.4071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97355,'amu*angstrom^2'), symmetry=1, barrier=(22.3838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975275,'amu*angstrom^2'), symmetry=1, barrier=(22.4235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973516,'amu*angstrom^2'), symmetry=1, barrier=(22.3831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0694023,0.0957645,-0.000125814,8.61177e-08,-2.36677e-11,-109824,28.3951], Tmin=(100,'K'), Tmax=(885.085,'K')), NASAPolynomial(coeffs=[14.3176,0.0307443,-1.56209e-05,3.11731e-09,-2.23392e-13,-112371,-39.2578], Tmin=(885.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-914.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(O)(F)C(F)C(F)F(5087)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {11,S}
9  C u1 p0 c0 {6,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-937.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,250,417,511,1155,1315,1456,3119,235,523,627,1123,1142,1372,1406,3097,3000,3100,440,815,1455,1000,180,183.502],'cm^-1')),
        HinderedRotor(inertia=(0.401823,'amu*angstrom^2'), symmetry=1, barrier=(9.23869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36477,'amu*angstrom^2'), symmetry=1, barrier=(9.22396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.921346,'amu*angstrom^2'), symmetry=1, barrier=(21.1835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60522,'amu*angstrom^2'), symmetry=1, barrier=(37.5961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.237792,0.100438,-0.000144561,1.1018e-07,-3.3674e-11,-112649,28.2948], Tmin=(100,'K'), Tmax=(799.562,'K')), NASAPolynomial(coeffs=[13.429,0.0320631,-1.62818e-05,3.21613e-09,-2.27876e-13,-114835,-34.5815], Tmin=(799.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-937.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC[C](F)C(F)C(F)F(4569)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-929.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.205249,'amu*angstrom^2'), symmetry=1, barrier=(4.71908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205499,'amu*angstrom^2'), symmetry=1, barrier=(4.72482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1224,'amu*angstrom^2'), symmetry=1, barrier=(25.8062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09951,0.0919707,-0.000258795,4.88045e-07,-3.69288e-10,-111836,15.5517], Tmin=(10,'K'), Tmax=(371.797,'K')), NASAPolynomial(coeffs=[6.57913,0.0449748,-3.06215e-05,9.74975e-09,-1.17324e-12,-112029,3.0959], Tmin=(371.797,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-929.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""FC[C](F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(CF)C(F)C(F)F(4095)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1088.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,334,437,422,648,895,1187,250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,235,523,627,1123,1142,1372,1406,3097,214.495,214.552,214.557],'cm^-1')),
        HinderedRotor(inertia=(0.394826,'amu*angstrom^2'), symmetry=1, barrier=(12.8961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.394842,'amu*angstrom^2'), symmetry=1, barrier=(12.8965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01664,'amu*angstrom^2'), symmetry=1, barrier=(33.2105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10154,0.0969379,-0.000128127,8.91619e-08,-2.50105e-11,-130822,30.0118], Tmin=(100,'K'), Tmax=(866.48,'K')), NASAPolynomial(coeffs=[13.863,0.0324717,-1.65261e-05,3.29586e-09,-2.35945e-13,-133242,-35.3577], Tmin=(866.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1088.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-F1sCsCs))"""),
)

species(
    label = 'CHF2-CHF(67)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-481.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.191999,'amu*angstrom^2'), symmetry=1, barrier=(4.41444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0324,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03194,0.0198962,-1.06494e-05,1.96763e-09,3.06708e-14,-57935.2,10.3725], Tmin=(10,'K'), Tmax=(1633.48,'K')), NASAPolynomial(coeffs=[15.0609,-0.000860555,2.67159e-06,-1.12649e-09,1.457e-13,-62372.2,-50.8003], Tmin=(1633.48,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-481.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[CH]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)CF(2940)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {8,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-471.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180],'cm^-1')),
        HinderedRotor(inertia=(0.625922,'amu*angstrom^2'), symmetry=1, barrier=(14.3912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628429,'amu*angstrom^2'), symmetry=1, barrier=(14.4488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78656,0.0175844,0.000121424,-4.80809e-07,5.06984e-10,-56694.5,10.4186], Tmin=(10,'K'), Tmax=(351.977,'K')), NASAPolynomial(coeffs=[6.00134,0.0208372,-1.35648e-05,4.29025e-09,-5.21054e-13,-57026.5,-0.454778], Tmin=(351.977,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-471.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""O[C](F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O[C](F)C(F)C(F)F(804)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {11,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-901.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.25657,'amu*angstrom^2'), symmetry=1, barrier=(5.89906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697821,'amu*angstrom^2'), symmetry=1, barrier=(16.0443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53517,'amu*angstrom^2'), symmetry=1, barrier=(35.2965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4724,0.0550086,-5.60381e-05,2.91701e-08,-6.17523e-12,-108475,14.0234], Tmin=(10,'K'), Tmax=(1075.24,'K')), NASAPolynomial(coeffs=[10.9758,0.0270952,-1.70977e-05,5.02636e-09,-5.61629e-13,-110089,-22.7203], Tmin=(1075.24,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-901.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""O[C](F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF2(82)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)([CH]F)CF(5088)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {3,S} {5,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-689.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17284,'amu*angstrom^2'), symmetry=1, barrier=(26.966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03354,'amu*angstrom^2'), symmetry=1, barrier=(23.7631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02415,'amu*angstrom^2'), symmetry=1, barrier=(23.5472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59965,0.0341675,5.98545e-05,-1.94464e-07,1.3783e-10,-82975.7,12.7457], Tmin=(10,'K'), Tmax=(539.264,'K')), NASAPolynomial(coeffs=[6.99729,0.0370429,-2.62432e-05,8.52514e-09,-1.03511e-12,-83750.4,-5.33307], Tmin=(539.264,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-689.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)([CH]F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(CF)[C](F)C(F)F(5089)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {5,S} {7,S} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1132.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,528,1116,1182,1331,1402,1494,3075,3110,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.664149,'amu*angstrom^2'), symmetry=1, barrier=(15.2701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664267,'amu*angstrom^2'), symmetry=1, barrier=(15.2728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66419,'amu*angstrom^2'), symmetry=1, barrier=(38.2631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664504,'amu*angstrom^2'), symmetry=1, barrier=(15.2782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.412424,0.10509,-0.000154745,1.20346e-07,-3.75416e-11,-136041,31.9895], Tmin=(100,'K'), Tmax=(783.504,'K')), NASAPolynomial(coeffs=[13.6414,0.0333409,-1.73803e-05,3.4639e-09,-2.46415e-13,-138244,-32.3829], Tmin=(783.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1132.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(CsCsCsF1s)"""),
)

species(
    label = 'OC(F)([CH]F)C(F)C(F)F(5090)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {7,S} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1123.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,250,417,511,1155,1315,1456,3119,235,523,627,1123,1142,1372,1406,3097,334,575,1197,1424,3202,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.477269,'amu*angstrom^2'), symmetry=1, barrier=(10.9734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478535,'amu*angstrom^2'), symmetry=1, barrier=(11.0025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.695,'amu*angstrom^2'), symmetry=1, barrier=(38.9715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675824,'amu*angstrom^2'), symmetry=1, barrier=(15.5385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.947473,0.118964,-0.000198693,1.71367e-07,-5.79277e-11,-134928,31.5877], Tmin=(100,'K'), Tmax=(808.842,'K')), NASAPolynomial(coeffs=[13.8178,0.0336748,-1.77708e-05,3.49191e-09,-2.43645e-13,-136915,-34.0326], Tmin=(808.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1123.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sO2sCs)(F1s)(H))"""),
)

species(
    label = 'OC(F)(CF)C(F)[C](F)F(5091)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1121.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,259,529,569,1128,1321,1390,3140,528,1116,1182,1331,1402,1494,3075,3110,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.763439,'amu*angstrom^2'), symmetry=1, barrier=(17.553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42897,'amu*angstrom^2'), symmetry=1, barrier=(32.8549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76201,'amu*angstrom^2'), symmetry=1, barrier=(17.5201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.761482,'amu*angstrom^2'), symmetry=1, barrier=(17.508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.402882,0.104715,-0.000151775,1.15495e-07,-3.52177e-11,-134730,31.2508], Tmin=(100,'K'), Tmax=(801.239,'K')), NASAPolynomial(coeffs=[14.0094,0.0327614,-1.70662e-05,3.40666e-09,-2.42832e-13,-137040,-35.086], Tmin=(801.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1121.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(CF)C(F)C(F)F(4573)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {5,D} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-998.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,235,523,627,1123,1142,1372,1406,3097,244,1102,1254,1379,1468,1475,3037,3083,375,552.5,462.5,1710,353.985,354.29,354.297],'cm^-1')),
        HinderedRotor(inertia=(0.121835,'amu*angstrom^2'), symmetry=1, barrier=(10.8491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121934,'amu*angstrom^2'), symmetry=1, barrier=(10.8486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345304,'amu*angstrom^2'), symmetry=1, barrier=(30.7562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59292,0.081801,-0.000108321,7.98812e-08,-2.42974e-11,-119944,26.7705], Tmin=(100,'K'), Tmax=(794.61,'K')), NASAPolynomial(coeffs=[10.1788,0.0335432,-1.72178e-05,3.4416e-09,-2.46395e-13,-121467,-17.2714], Tmin=(794.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-998.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs)"""),
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
    label = 'OC(=CC(F)F)CF(5092)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {7,S} {13,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-794.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,212,931,1104,1251,1325,1474,3102,3155,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,3010,987.5,1337.5,450,1655,276.813,279.802],'cm^-1')),
        HinderedRotor(inertia=(0.307511,'amu*angstrom^2'), symmetry=1, barrier=(16.6417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29826,'amu*angstrom^2'), symmetry=1, barrier=(16.6469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305149,'amu*angstrom^2'), symmetry=1, barrier=(16.6183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817576,0.0710594,-7.44384e-05,4.00845e-08,-8.69738e-12,-95404.7,24.3585], Tmin=(100,'K'), Tmax=(1107.34,'K')), NASAPolynomial(coeffs=[13.6497,0.0247063,-1.16484e-05,2.28206e-09,-1.62854e-13,-98246.6,-38.8574], Tmin=(1107.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-794.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCFFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=C(O)C(F)C(F)F(5093)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {7,S} {13,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-801.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,284,328,853,1146,1135,1297,3239,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.890757,'amu*angstrom^2'), symmetry=1, barrier=(20.4803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.89033,'amu*angstrom^2'), symmetry=1, barrier=(20.4704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890464,'amu*angstrom^2'), symmetry=1, barrier=(20.4735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387467,0.0763784,-8.19687e-05,4.36726e-08,-9.14204e-12,-96315.9,24.4694], Tmin=(100,'K'), Tmax=(1164.2,'K')), NASAPolynomial(coeffs=[16.9676,0.0194113,-8.56959e-06,1.64105e-09,-1.16135e-13,-100176,-58.0409], Tmin=(1164.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-801.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'OC(CF)=C(F)C(F)F(5094)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {7,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-959.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,212,931,1104,1251,1325,1474,3102,3155,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,323,467,575,827,1418,292.626,293.737],'cm^-1')),
        HinderedRotor(inertia=(0.136743,'amu*angstrom^2'), symmetry=1, barrier=(8.5868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143426,'amu*angstrom^2'), symmetry=1, barrier=(8.61531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.448706,'amu*angstrom^2'), symmetry=1, barrier=(28.3811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.412674,0.0857905,-0.000118452,8.93959e-08,-2.75067e-11,-115249,26.5587], Tmin=(100,'K'), Tmax=(789.343,'K')), NASAPolynomial(coeffs=[10.9993,0.0321422,-1.65016e-05,3.28949e-09,-2.34786e-13,-116920,-22.011], Tmin=(789.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-959.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'OC(F)(C=CF)CF(5095)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {13,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {7,D} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-799.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1968,'amu*angstrom^2'), symmetry=1, barrier=(27.5167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19669,'amu*angstrom^2'), symmetry=1, barrier=(27.5142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19627,'amu*angstrom^2'), symmetry=1, barrier=(27.5046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835473,0.0733281,-8.12751e-05,4.76704e-08,-1.14133e-11,-96098.3,23.6987], Tmin=(100,'K'), Tmax=(1001.87,'K')), NASAPolynomial(coeffs=[12.1392,0.0281976,-1.37059e-05,2.7086e-09,-1.93872e-13,-98363.3,-30.8567], Tmin=(1001.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-799.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH)"""),
)

species(
    label = 'OC(F)(CF)C(F)=CF(3053)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-973.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17748,'amu*angstrom^2'), symmetry=1, barrier=(27.0725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18628,'amu*angstrom^2'), symmetry=1, barrier=(27.275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18086,'amu*angstrom^2'), symmetry=1, barrier=(27.1503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23905,0.0886561,-0.000119432,8.48404e-08,-2.42239e-11,-116917,26.2017], Tmin=(100,'K'), Tmax=(853.108,'K')), NASAPolynomial(coeffs=[12.9607,0.0290067,-1.45501e-05,2.87867e-09,-2.04952e-13,-119087,-33.1517], Tmin=(853.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-973.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OC(=CF)C(F)C(F)F(5096)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-966.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,284,328,853,1146,1135,1297,3239,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,194,682,905,1196,1383,3221,205.57,205.603],'cm^-1')),
        HinderedRotor(inertia=(0.444586,'amu*angstrom^2'), symmetry=1, barrier=(13.3342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44461,'amu*angstrom^2'), symmetry=1, barrier=(13.3293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03057,'amu*angstrom^2'), symmetry=1, barrier=(30.9365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.411159,0.0846116,-0.00010648,7.02297e-08,-1.87237e-11,-116169,25.863], Tmin=(100,'K'), Tmax=(908.367,'K')), NASAPolynomial(coeffs=[13.0136,0.0291171,-1.48419e-05,2.97507e-09,-2.1403e-13,-118458,-33.7254], Tmin=(908.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-966.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'OC(F)(C=C(F)F)CF(3081)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1000.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19316,'amu*angstrom^2'), symmetry=1, barrier=(27.4331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19868,'amu*angstrom^2'), symmetry=1, barrier=(27.5599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19874,'amu*angstrom^2'), symmetry=1, barrier=(27.5615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504253,0.0823626,-0.000101676,6.62963e-08,-1.75387e-11,-120266,25.9004], Tmin=(100,'K'), Tmax=(913.774,'K')), NASAPolynomial(coeffs=[12.5823,0.0294905,-1.48824e-05,2.9725e-09,-2.13496e-13,-122474,-31.2798], Tmin=(913.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1000.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF)"""),
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
    E0 = (-543.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-254.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-264.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-410.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-381.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-215.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-368.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-409.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-222.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-433.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-465.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-494.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-494.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-454.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-141.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-465.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-170.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-369.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-329.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-342.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-312.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-335.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-372.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-347.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-423.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-415.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-417.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-390.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-382.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-380.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-602.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-266.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-273.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-487.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-271.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-531.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-491.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-509.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)(CF)C(F)C(F)F(3507)'],
    products = ['HF(38)', 'CHFCHF[Z](59)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.97084e+10,'s^-1'), n=0.978818, Ea=(250.083,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['OCF(781)', 'F[C]C(F)C(F)F(793)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(185.469,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]CF(126)', 'OC(F)C(F)F(778)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(191.698,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2F2(41)', 'OC(F)([C]F)CF(4240)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.67744e-07,'m^3/(mol*s)'), n=3.38172, Ea=(55.0622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]C(F)F(125)', 'OC(F)CF(773)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.33872e-07,'m^3/(mol*s)'), n=3.38172, Ea=(46.6958,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHF(40)', 'OC(F)(CF)C(F)F(780)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.7066e-06,'m^3/(mol*s)'), n=3.19155, Ea=(224.993,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H2(8)', 'OC(F)([C]F)C(F)C(F)F(5080)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CHF(40)', 'OC(F)C(F)C(F)F(384)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.79957e-05,'m^3/(mol*s)'), n=3.10993, Ea=(22.9309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(S)(25)', 'OC(F)(F)C(F)C(F)F(2802)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(161.575,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2(43)', 'OC(F)(CF)CF(4239)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.67207e-10,'m^3/(mol*s)'), n=4.72997, Ea=(129.609,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['HF(38)', 'OC(F)(CF)C(F)[C]F(5081)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4487,'m^3/(mol*s)'), n=1.52711, Ea=(27.7169,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.02090651070256912, var=1.7269912187795557, Tref=1000.0, N=2, data_mean=0.0, correlation='HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R',), comment="""Estimated from node HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OC(F)(CF)C(F)C(F)F(3507)'],
    products = ['OC(C(F)F)C(F)(F)CF(5082)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OC(F)(CF)C(F)C(F)F(3507)'],
    products = ['OCC(F)(F)C(F)C(F)F(5083)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H2O(3)', 'FCC(F)=C(F)C(F)F(4565)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OF(482)', 'FCC(F)=CC(F)F(4123)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H2O(3)', 'FC=C(F)C(F)C(F)F(3694)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OF(482)', 'C=C(F)C(F)C(F)F(4566)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1250,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OC(F)F(368)', 'FC=C(F)CF(486)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/monosub_Cd/disub;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F(37)', 'O[C](CF)C(F)C(F)F(5084)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F(37)', 'OC(F)([CH]C(F)F)CF(5085)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F(37)', 'OC(F)(CF)C(F)[CH]F(5086)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['F(37)', '[CH2]C(O)(F)C(F)C(F)F(5087)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['OH(7)', 'FC[C](F)C(F)C(F)F(4569)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', '[O]C(F)(CF)C(F)C(F)F(4095)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CHF2-CHF(67)', 'O[C](F)CF(2940)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -8.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2F(46)', 'O[C](F)C(F)C(F)F(804)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CHF2(82)', 'OC(F)([CH]F)CF(5088)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -10.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(5)', 'OC(F)(CF)[C](F)C(F)F(5089)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.261562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(5)', 'OC(F)([CH]F)C(F)C(F)F(5090)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(5)', 'OC(F)(CF)C(F)[C](F)F(5091)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['HF(38)', 'O=C(CF)C(F)C(F)F(4573)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(147.999,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F2(78)', 'OC(=CC(F)F)CF(5092)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.06825,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F2(78)', 'C=C(O)C(F)C(F)F(5093)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.59832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction34',
    reactants = ['HF(38)', 'OC(CF)=C(F)C(F)F(5094)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(223.632,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F2(78)', 'OC(F)(C=CF)CF(5095)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.45545,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction36',
    reactants = ['HF(38)', 'OC(F)(CF)C(F)=CF(3053)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(193.66,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction37',
    reactants = ['HF(38)', 'OC(=CF)C(F)C(F)F(5096)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(227.026,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction38',
    reactants = ['HF(38)', 'OC(F)(C=C(F)F)CF(3081)'],
    products = ['OC(F)(CF)C(F)C(F)F(3507)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(243,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #639',
    isomers = [
        'OC(F)(CF)C(F)C(F)F(3507)',
    ],
    reactants = [
        ('HF(38)', 'CHFCHF[Z](59)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #639',
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

