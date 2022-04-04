species(
    label = 'C=C(O)C(C)(F)OF(5731)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,S} {14,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-494.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16938,'amu*angstrom^2'), symmetry=1, barrier=(26.8863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17456,'amu*angstrom^2'), symmetry=1, barrier=(27.0054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17384,'amu*angstrom^2'), symmetry=1, barrier=(26.9889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17334,'amu*angstrom^2'), symmetry=1, barrier=(26.9773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0610827,0.0869222,-9.74213e-05,5.40292e-08,-1.17438e-11,-59300.3,24.4619], Tmin=(100,'K'), Tmax=(1124.04,'K')), NASAPolynomial(coeffs=[18.4821,0.0209346,-9.36287e-06,1.80197e-09,-1.27882e-13,-63469,-67.1667], Tmin=(1124.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)OF(1265)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u0 p2 c0 {1,S} {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-200.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.3132,'amu*angstrom^2'), symmetry=1, barrier=(30.193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31177,'amu*angstrom^2'), symmetry=1, barrier=(30.1601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87494,0.00722414,0.000106356,-2.21274e-07,1.3152e-10,-24092.6,10.0683], Tmin=(10,'K'), Tmax=(591.099,'K')), NASAPolynomial(coeffs=[5.13645,0.0281847,-2.16884e-05,7.56312e-09,-9.70859e-13,-24757.1,0.286506], Tmin=(591.099,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-200.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CDC(O)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'COF(400)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-101.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.263855,'amu*angstrom^2'), symmetry=1, barrier=(6.06655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0324,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98349,0.00100505,4.5348e-05,-8.07633e-08,4.67681e-11,-12262,7.03919], Tmin=(10,'K'), Tmax=(445.12,'K')), NASAPolynomial(coeffs=[2.13239,0.0176377,-1.06953e-05,3.16422e-09,-3.63845e-13,-12097.1,14.4716], Tmin=(445.12,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-101.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""COF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(O)[C]F(2149)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u0 p1 c0 {1,S} {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-45.5134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,315,622,1128],'cm^-1')),
        HinderedRotor(inertia=(0.955095,'amu*angstrom^2'), symmetry=1, barrier=(21.9595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949253,'amu*angstrom^2'), symmetry=1, barrier=(21.8252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04376,0.0341448,-1.44809e-05,-1.42949e-08,1.01497e-11,-5395.46,15.8839], Tmin=(100,'K'), Tmax=(949.386,'K')), NASAPolynomial(coeffs=[14.4866,0.00352096,-5.411e-07,1.03954e-10,-1.11514e-14,-8740.58,-48.6739], Tmin=(949.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.5134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(CJ2_singlet-FC)"""),
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
    label = 'C=C(O)C(F)OF(2868)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-439.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,236,527,855,1015,1182,1348,3236,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06105,'amu*angstrom^2'), symmetry=1, barrier=(24.3956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06138,'amu*angstrom^2'), symmetry=1, barrier=(24.4031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06335,'amu*angstrom^2'), symmetry=1, barrier=(24.4484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3528.52,'J/mol'), sigma=(5.6427,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.15 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672922,0.0683159,-7.63262e-05,4.11652e-08,-8.58402e-12,-52762.4,21.7001], Tmin=(100,'K'), Tmax=(1179.64,'K')), NASAPolynomial(coeffs=[17.1591,0.0124136,-5.24208e-06,9.92413e-10,-7.02288e-14,-56652,-60.5596], Tmin=(1179.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-439.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)C(=C)F(1180)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,D} {9,S} {10,S}
6  C u0 p0 c0 {4,D} {7,S} {8,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-315.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,280,518,736,852,873,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.970345,'amu*angstrom^2'), symmetry=1, barrier=(22.3101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.971457,'amu*angstrom^2'), symmetry=1, barrier=(22.3357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19555,0.0527577,-3.70345e-05,1.69008e-09,5.46285e-12,-37784.3,17.8739], Tmin=(100,'K'), Tmax=(962.139,'K')), NASAPolynomial(coeffs=[16.3428,0.0104306,-3.23346e-06,5.72441e-10,-4.19315e-14,-41654.6,-59.5845], Tmin=(962.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(OF)=C(O)CF(5806)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,S} {14,S}
5  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {4,S} {6,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-456.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,203.036,203.039],'cm^-1')),
        HinderedRotor(inertia=(0.00408913,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543191,'amu*angstrom^2'), symmetry=1, barrier=(15.8901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543193,'amu*angstrom^2'), symmetry=1, barrier=(15.8901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543168,'amu*angstrom^2'), symmetry=1, barrier=(15.8901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406656,0.081574,-9.54141e-05,5.83507e-08,-1.43085e-11,-54833.7,25.5993], Tmin=(100,'K'), Tmax=(989.531,'K')), NASAPolynomial(coeffs=[13.9645,0.0267693,-1.23383e-05,2.38156e-09,-1.6832e-13,-57516.9,-39.6671], Tmin=(989.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs)"""),
)

species(
    label = 'CC(F)=C(O)COF(5807)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-444.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,323,467,575,827,1418,302.246,302.284],'cm^-1')),
        HinderedRotor(inertia=(0.15115,'amu*angstrom^2'), symmetry=1, barrier=(9.79723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151144,'amu*angstrom^2'), symmetry=1, barrier=(9.79755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151093,'amu*angstrom^2'), symmetry=1, barrier=(9.79776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151133,'amu*angstrom^2'), symmetry=1, barrier=(9.79764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538739,0.0810939,-0.000100549,6.81008e-08,-1.8796e-11,-53376,26.0989], Tmin=(100,'K'), Tmax=(877.19,'K')), NASAPolynomial(coeffs=[11.5479,0.0308923,-1.47044e-05,2.85907e-09,-2.02078e-13,-55307.4,-25.5713], Tmin=(877.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'CCC(O)=C(F)OF(5808)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {2,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {1,S} {4,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-440.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,326,540,652,719,1357,274.881,275.567],'cm^-1')),
        HinderedRotor(inertia=(0.00223109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127036,'amu*angstrom^2'), symmetry=1, barrier=(6.87316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127369,'amu*angstrom^2'), symmetry=1, barrier=(6.86353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341662,'amu*angstrom^2'), symmetry=1, barrier=(18.3107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.740119,0.076952,-8.93925e-05,5.76405e-08,-1.53926e-11,-52864.5,26.1969], Tmin=(100,'K'), Tmax=(898.305,'K')), NASAPolynomial(coeffs=[10.6199,0.0329607,-1.59383e-05,3.12954e-09,-2.22741e-13,-54639.6,-20.4083], Tmin=(898.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-440.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[CH2]C([O])C(C)(F)OF(5772)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-156.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,264.546,264.548,440.807],'cm^-1')),
        HinderedRotor(inertia=(0.397961,'amu*angstrom^2'), symmetry=1, barrier=(19.7641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397962,'amu*angstrom^2'), symmetry=1, barrier=(19.7641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210702,'amu*angstrom^2'), symmetry=1, barrier=(10.4639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397935,'amu*angstrom^2'), symmetry=1, barrier=(19.7641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209815,0.0872739,-0.000104633,6.50654e-08,-1.62441e-11,-18705.6,26.9483], Tmin=(100,'K'), Tmax=(971.067,'K')), NASAPolynomial(coeffs=[14.4796,0.0284929,-1.38329e-05,2.72708e-09,-1.94833e-13,-21476.9,-41.476], Tmin=(971.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C([CH2])(F)OF(5809)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {14,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {5,S} {10,S} {11,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-176.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,316,385,515,654,689,1295,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,256.92,256.92],'cm^-1')),
        HinderedRotor(inertia=(0.457004,'amu*angstrom^2'), symmetry=1, barrier=(21.4064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15426,'amu*angstrom^2'), symmetry=1, barrier=(7.22566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154261,'amu*angstrom^2'), symmetry=1, barrier=(7.22566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655972,'amu*angstrom^2'), symmetry=1, barrier=(30.726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154261,'amu*angstrom^2'), symmetry=1, barrier=(7.22566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.406022,0.103991,-0.000161069,1.30938e-07,-4.22339e-11,-21033.4,29.5898], Tmin=(100,'K'), Tmax=(779.735,'K')), NASAPolynomial(coeffs=[13.6927,0.0298194,-1.48329e-05,2.87208e-09,-1.99844e-13,-23175.9,-34.5607], Tmin=(779.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Csj(Cs-F1sO2sCs)(H)(H)) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(F)(OF)[C](C)O(5774)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,S} {14,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {5,S} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-211.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,222.063,222.071],'cm^-1')),
        HinderedRotor(inertia=(0.247595,'amu*angstrom^2'), symmetry=1, barrier=(8.66235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24756,'amu*angstrom^2'), symmetry=1, barrier=(8.66309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247675,'amu*angstrom^2'), symmetry=1, barrier=(8.66265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247597,'amu*angstrom^2'), symmetry=1, barrier=(8.66277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990792,'amu*angstrom^2'), symmetry=1, barrier=(34.6479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.465585,0.108282,-0.000181303,1.58287e-07,-5.37238e-11,-25239.1,28.5898], Tmin=(100,'K'), Tmax=(836.345,'K')), NASAPolynomial(coeffs=[11.712,0.0333023,-1.68062e-05,3.23331e-09,-2.224e-13,-26690.6,-24.484], Tmin=(836.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'CC(=O)C(C)(F)OF(5733)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {8,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {4,D} {5,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-526.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.89656,'amu*angstrom^2'), symmetry=1, barrier=(20.6137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8967,'amu*angstrom^2'), symmetry=1, barrier=(20.6169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896423,'amu*angstrom^2'), symmetry=1, barrier=(20.6105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.89577,'amu*angstrom^2'), symmetry=1, barrier=(20.5955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3575.07,'J/mol'), sigma=(5.88724,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.42 K, Pc=39.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.765051,0.0761366,-8.29549e-05,4.91947e-08,-1.21053e-11,-63239.2,23.6675], Tmin=(100,'K'), Tmax=(968.33,'K')), NASAPolynomial(coeffs=[11.2064,0.0330043,-1.61393e-05,3.19332e-09,-2.28561e-13,-65261.3,-26.3698], Tmin=(968.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-526.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs)"""),
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
    label = '[CH2]C(O)=C(C)OF(3657)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,S} {13,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {7,S}
7  C u1 p0 c0 {6,S} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-142.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,264.17],'cm^-1')),
        HinderedRotor(inertia=(0.286132,'amu*angstrom^2'), symmetry=1, barrier=(14.1689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00241572,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286118,'amu*angstrom^2'), symmetry=1, barrier=(14.1689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286131,'amu*angstrom^2'), symmetry=1, barrier=(14.169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367626,0.0751663,-8.61458e-05,4.95786e-08,-1.10625e-11,-16979.4,25.2794], Tmin=(100,'K'), Tmax=(1105.94,'K')), NASAPolynomial(coeffs=[16.8289,0.0156291,-5.39536e-06,9.02103e-10,-5.91727e-14,-20620.5,-55.7943], Tmin=(1105.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(O)C(C)([O])F(5810)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {6,S} {13,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-398.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23191,'amu*angstrom^2'), symmetry=1, barrier=(28.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2269,'amu*angstrom^2'), symmetry=1, barrier=(28.2088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23187,'amu*angstrom^2'), symmetry=1, barrier=(28.3232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.31792,0.0736401,-7.55931e-05,3.81638e-08,-7.49694e-12,-47781.5,22.7376], Tmin=(100,'K'), Tmax=(1247.15,'K')), NASAPolynomial(coeffs=[18.2588,0.016098,-6.38482e-06,1.16844e-09,-8.09549e-14,-52256.5,-67.7794], Tmin=(1247.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-398.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
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
    label = 'C=C(O)[C](C)F(5811)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {12,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u1 p0 c0 {1,S} {3,S} {5,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-273.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,346,659,817,1284,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.975415,'amu*angstrom^2'), symmetry=1, barrier=(22.4267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976253,'amu*angstrom^2'), symmetry=1, barrier=(22.446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97738,'amu*angstrom^2'), symmetry=1, barrier=(22.4719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0882,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1122,0.0549338,-3.88985e-05,6.02664e-09,2.89989e-12,-32785.7,19.8531], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[15.3686,0.0155158,-5.80513e-06,1.07136e-09,-7.6474e-14,-36530,-53.3593], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCdCsF1s)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C(C)(F)OF(5812)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-61.6253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14423,'amu*angstrom^2'), symmetry=1, barrier=(26.308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1435,'amu*angstrom^2'), symmetry=1, barrier=(26.2914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14076,'amu*angstrom^2'), symmetry=1, barrier=(26.2284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09166,0.069531,-8.73609e-05,6.12084e-08,-1.77444e-11,-7311.99,21.2661], Tmin=(100,'K'), Tmax=(830.914,'K')), NASAPolynomial(coeffs=[9.45145,0.0292885,-1.4716e-05,2.92533e-09,-2.09164e-13,-8701.29,-17.5169], Tmin=(830.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.6253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'C=C([O])C(C)(F)OF(5778)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-356.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15652,'amu*angstrom^2'), symmetry=1, barrier=(26.5907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1571,'amu*angstrom^2'), symmetry=1, barrier=(26.6041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15708,'amu*angstrom^2'), symmetry=1, barrier=(26.6036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467866,0.0809127,-9.73482e-05,6.07617e-08,-1.519e-11,-42750.1,24.0735], Tmin=(100,'K'), Tmax=(971.291,'K')), NASAPolynomial(coeffs=[13.8342,0.0258669,-1.23385e-05,2.41306e-09,-1.71643e-13,-45346.6,-40.0217], Tmin=(971.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)=C(F)OF(3081)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {2,S} {6,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u1 p0 c0 {5,S} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-258.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,350,440,435,1725,326,540,652,719,1357,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.147273,'amu*angstrom^2'), symmetry=1, barrier=(3.38611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14721,'amu*angstrom^2'), symmetry=1, barrier=(3.38465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03963,0.0695247,-0.000104803,8.28382e-08,-2.5979e-11,-31038,22.8848], Tmin=(100,'K'), Tmax=(783.408,'K')), NASAPolynomial(coeffs=[10.7637,0.0198586,-9.67615e-06,1.86041e-09,-1.29139e-13,-32561.1,-21.6514], Tmin=(783.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=[C]O(142)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989113,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16243,0.0134242,5.5648e-06,-1.95529e-08,9.36453e-12,12455.2,10.1543], Tmin=(100,'K'), Tmax=(925.605,'K')), NASAPolynomial(coeffs=[8.19866,0.00453477,-8.93536e-07,1.26104e-10,-9.46689e-15,10971.4,-16.7325], Tmin=(925.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(F)(OF)C(=C)O(5813)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-282.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19546,'amu*angstrom^2'), symmetry=1, barrier=(27.486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19624,'amu*angstrom^2'), symmetry=1, barrier=(27.5039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19459,'amu*angstrom^2'), symmetry=1, barrier=(27.4661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19525,'amu*angstrom^2'), symmetry=1, barrier=(27.4811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.174129,0.091615,-0.00011532,7.071e-08,-1.68419e-11,-33849.8,25.5043], Tmin=(100,'K'), Tmax=(1032.81,'K')), NASAPolynomial(coeffs=[18.8589,0.017901,-8.26125e-06,1.6046e-09,-1.14274e-13,-37781.3,-66.9336], Tmin=(1032.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CJCO)"""),
)

species(
    label = '[CH]=C(O)C(C)(F)OF(5814)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-247.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.18173,'amu*angstrom^2'), symmetry=1, barrier=(27.1704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18394,'amu*angstrom^2'), symmetry=1, barrier=(27.2212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18644,'amu*angstrom^2'), symmetry=1, barrier=(27.2785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18587,'amu*angstrom^2'), symmetry=1, barrier=(27.2655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0753705,0.0897022,-0.000111329,6.77073e-08,-1.60498e-11,-29583,24.9197], Tmin=(100,'K'), Tmax=(1035.88,'K')), NASAPolynomial(coeffs=[18.2633,0.0188875,-8.78438e-06,1.71125e-09,-1.22006e-13,-33382.3,-64.1997], Tmin=(1035.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(O)C(C)(F)OF(5815)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {14,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u0 p1 c0 {6,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-147.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,333,384,448,608,1254,1480,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.029961,0.0943592,-0.000134423,1.02895e-07,-3.17668e-11,-17601.8,27.2931], Tmin=(100,'K'), Tmax=(790.112,'K')), NASAPolynomial(coeffs=[12.3204,0.0321372,-1.62949e-05,3.22176e-09,-2.28516e-13,-19544,-29.1054], Tmin=(790.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=C(O)C(C)=O(1008)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {3,S} {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-333.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,350,440,435,1725,2950,3100,1380,975,1025,1650,336.014],'cm^-1')),
        HinderedRotor(inertia=(0.167823,'amu*angstrom^2'), symmetry=1, barrier=(13.5284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169232,'amu*angstrom^2'), symmetry=1, barrier=(13.5218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167866,'amu*angstrom^2'), symmetry=1, barrier=(13.533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34047,0.0507425,-3.5558e-05,7.1326e-09,1.49076e-12,-39954.3,19.6241], Tmin=(100,'K'), Tmax=(1048.53,'K')), NASAPolynomial(coeffs=[13.8723,0.0167055,-6.56506e-06,1.22369e-09,-8.67526e-14,-43339.3,-45.0383], Tmin=(1048.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)C(=C)OF(2861)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {12,S}
3  O u0 p2 c0 {1,S} {5,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-176.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.875599,'amu*angstrom^2'), symmetry=1, barrier=(20.1317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87561,'amu*angstrom^2'), symmetry=1, barrier=(20.132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875844,'amu*angstrom^2'), symmetry=1, barrier=(20.1374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3673.09,'J/mol'), sigma=(5.85032,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.73 K, Pc=41.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305392,0.076102,-9.13118e-05,5.35484e-08,-1.20534e-11,-21037.4,22.2077], Tmin=(100,'K'), Tmax=(1101.23,'K')), NASAPolynomial(coeffs=[17.9715,0.0119336,-3.9078e-06,6.3573e-10,-4.12645e-14,-24928.4,-64.7248], Tmin=(1101.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    E0 = (-315.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (203.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (188.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (126.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-59.9878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-19.0393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-90.5922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (71.3122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (13.6425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (34.6715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-29.1719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-141.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (82.8963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-178.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-23.4205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (114.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (2.74133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (24.6871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (82.8737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (76.5256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (112.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (10.8627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-115.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-91.7879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C(C)(F)OF(5731)'],
    products = ['HF(38)', 'CH2CO(28)', 'CC(=O)F(499)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(31.4635,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[C]F(124)', 'C=C(O)OF(1265)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(175.332,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['COF(400)', 'C=C(O)[C]F(2149)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(188.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(S)(25)', 'C=C(O)C(F)OF(2868)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.16448e+06,'m^3/(mol*s)'), n=-0.144618, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(169)', 'C=C(O)C(=C)F(1180)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1250,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC(OF)=C(O)CF(5806)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(290.516,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(F)=C(O)COF(5807)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(206.785,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCC(O)=C(F)OF(5808)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3431.13,'s^-1'), n=2.77015, Ea=(364.375,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_2R!H->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([O])C(C)(F)OF(5772)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)C([CH2])(F)OF(5809)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(F)(OF)[C](C)O(5774)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C(C)(F)OF(5731)'],
    products = ['CC(=O)C(C)(F)OF(5733)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.913,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[CH2]C(O)=C(C)OF(3657)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(4.88476,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'C=C(O)C(C)([O])F(5810)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]F(160)', 'C=C(O)[C](C)F(5811)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(7)', 'C=[C]C(C)(F)OF(5812)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'C=C([O])C(C)(F)OF(5778)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(19)', '[CH2]C(O)=C(F)OF(3081)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.29446e+18,'m^3/(mol*s)'), n=-4.19701, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R
Ea raised from -7.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]O(142)', 'C[C](F)OF(5141)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[CH2]C(F)(OF)C(=C)O(5813)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[CH]=C(O)C(C)(F)OF(5814)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -5.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(O)C(C)(F)OF(5815)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(10.9316,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F2(77)', 'C=C(O)C(C)=O(1008)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(78.4819,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'C=C(O)C(=C)OF(2861)'],
    products = ['C=C(O)C(C)(F)OF(5731)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(217.968,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

network(
    label = 'PDepNetwork #2043',
    isomers = [
        'C=C(O)C(C)(F)OF(5731)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', 'CC(=O)F(499)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2043',
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

