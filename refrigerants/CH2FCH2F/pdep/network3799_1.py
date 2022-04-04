species(
    label = 'C=C(F)OC(O)(F)CF(11594)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-979.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.934961,'amu*angstrom^2'), symmetry=1, barrier=(21.4966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934394,'amu*angstrom^2'), symmetry=1, barrier=(21.4836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935124,'amu*angstrom^2'), symmetry=1, barrier=(21.5003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934792,'amu*angstrom^2'), symmetry=1, barrier=(21.4927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.112632,0.0909268,-0.00010506,6.02519e-08,-1.36297e-11,-117603,27.0853], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[17.9789,0.0237083,-1.14044e-05,2.25608e-09,-1.62084e-13,-121499,-61.5308], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-979.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH)"""),
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
    label = 'O=C(F)CF(476)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212711,-2.78626e-07,3.40268e-10,-7.76035e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.07,'K')), NASAPolynomial(coeffs=[2.78817,0.000587632,1.59015e-07,-5.52748e-11,4.34317e-15,-596.149,0.11269], Tmin=(1959.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(F)OC(O)(F)[C]F(11771)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {10,S} {11,S}
9  C u0 p1 c0 {3,S} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-627.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,326,540,652,719,1357,2950,3100,1380,975,1025,1650,617,898,1187,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08902,'amu*angstrom^2'), symmetry=1, barrier=(25.0388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08888,'amu*angstrom^2'), symmetry=1, barrier=(25.0355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08939,'amu*angstrom^2'), symmetry=1, barrier=(25.0471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08887,'amu*angstrom^2'), symmetry=1, barrier=(25.0352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0615704,0.0904798,-0.000122398,8.16438e-08,-2.13541e-11,-75298.9,27.4065], Tmin=(100,'K'), Tmax=(937.871,'K')), NASAPolynomial(coeffs=[16.464,0.0205216,-1.05051e-05,2.10476e-09,-1.51353e-13,-78375.5,-50.673], Tmin=(937.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-627.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(CdCFO) + group(Cds-CdsHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'C=C(F)OC(O)F(4123)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-760.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,326,540,652,719,1357,2950,3100,1380,975,1025,1650,191.121,191.121,191.121],'cm^-1')),
        HinderedRotor(inertia=(0.724308,'amu*angstrom^2'), symmetry=1, barrier=(18.7745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724309,'amu*angstrom^2'), symmetry=1, barrier=(18.7745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724308,'amu*angstrom^2'), symmetry=1, barrier=(18.7745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3528.52,'J/mol'), sigma=(5.6427,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.15 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.859677,0.0643527,-6.9315e-05,3.63684e-08,-7.41528e-12,-91404.9,21.9654], Tmin=(100,'K'), Tmax=(1202.29,'K')), NASAPolynomial(coeffs=[16.1953,0.013331,-5.65892e-06,1.071e-09,-7.55965e-14,-95092.5,-54.8451], Tmin=(1202.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-760.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFHOO) + group(CdCFO) + group(Cds-CdsHH)"""),
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
    label = 'C=C(F)OC(O)(F)F(3586)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-993.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180,282.928],'cm^-1')),
        HinderedRotor(inertia=(0.725805,'amu*angstrom^2'), symmetry=1, barrier=(16.6877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889931,'amu*angstrom^2'), symmetry=1, barrier=(20.4613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890036,'amu*angstrom^2'), symmetry=1, barrier=(20.4637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3527.82,'J/mol'), sigma=(5.84146,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.04 K, Pc=40.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66652,0.0703466,-7.92621e-05,4.3119e-08,-9.1237e-12,-119317,23.3172], Tmin=(100,'K'), Tmax=(1157.51,'K')), NASAPolynomial(coeffs=[16.8749,0.0143352,-6.6774e-06,1.31373e-09,-9.45053e-14,-123070,-57.2494], Tmin=(1157.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-993.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(F)OCC(O)(F)F(11772)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {14,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-990.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0509258,0.0886127,-0.000103516,6.10063e-08,-1.42684e-11,-119044,27.3585], Tmin=(100,'K'), Tmax=(1039.57,'K')), NASAPolynomial(coeffs=[16.5164,0.0252582,-1.21018e-05,2.38378e-09,-1.70643e-13,-122467,-52.7168], Tmin=(1039.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-990.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFFO) + group(CdCFO) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(F)OC(F)(F)CO(11773)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {14,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-990.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0509258,0.0886127,-0.000103516,6.10063e-08,-1.42684e-11,-119044,27.3585], Tmin=(100,'K'), Tmax=(1039.57,'K')), NASAPolynomial(coeffs=[16.5164,0.0252582,-1.21018e-05,2.38378e-09,-1.70643e-13,-122467,-52.7168], Tmin=(1039.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-990.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFFO) + group(Cs-CsOsHH) + group(CdCFO) + group(Cds-CdsHH)"""),
)

species(
    label = 'OC(F)=CF(3655)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {2,S} {4,D} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-508.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.537415,'amu*angstrom^2'), symmetry=1, barrier=(12.3562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3292.1,'J/mol'), sigma=(5.16209,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.22 K, Pc=54.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86684,0.00827899,9.49292e-05,-2.30798e-07,1.56874e-10,-61152.3,9.15162], Tmin=(10,'K'), Tmax=(532.988,'K')), NASAPolynomial(coeffs=[6.38193,0.0179589,-1.26771e-05,4.31723e-09,-5.56992e-13,-61826,-5.20451], Tmin=(532.988,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-508.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(O)F(213)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-348.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.540544,'amu*angstrom^2'), symmetry=1, barrier=(12.4282,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3213.72,'J/mol'), sigma=(5.20847,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.98 K, Pc=51.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93506,0.00386927,7.15726e-05,-1.45436e-07,8.82281e-11,-41920.4,8.43739], Tmin=(10,'K'), Tmax=(555.806,'K')), NASAPolynomial(coeffs=[3.62097,0.0208063,-1.37458e-05,4.40936e-09,-5.41481e-13,-42112.2,7.7289], Tmin=(555.806,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-348.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CDC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(F)OF(479)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-208.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,326,540,652,719,1357,2950,3100,1380,975,1025,1650,490.191],'cm^-1')),
        HinderedRotor(inertia=(0.13294,'amu*angstrom^2'), symmetry=1, barrier=(22.668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90632,0.00544022,8.64998e-05,-1.77457e-07,1.05368e-10,-25015.8,10.2412], Tmin=(10,'K'), Tmax=(583.792,'K')), NASAPolynomial(coeffs=[4.36855,0.0241832,-1.79543e-05,6.11266e-09,-7.72989e-13,-25443.1,5.06228], Tmin=(583.792,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-208.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CDC(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787929,2.90875e-06,-1.47516e-09,2.12833e-13,-30281.6,-0.311362], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81028e-07,1.07243e-10,-5.79385e-15,-29958.6,5.9104], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(F)OC(F)=CF(10044)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {6,D} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-567.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,390,483,597,572,732,631,807,1275,1439,2950,3100,1380,975,1025,1650,194,682,905,1196,1383,3221,180,180,1402.48],'cm^-1')),
        HinderedRotor(inertia=(0.195916,'amu*angstrom^2'), symmetry=1, barrier=(4.50449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196507,'amu*angstrom^2'), symmetry=1, barrier=(4.51808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17587,0.0703703,-0.000113535,1.06708e-07,-3.98682e-11,-68157,24.7626], Tmin=(100,'K'), Tmax=(798.05,'K')), NASAPolynomial(coeffs=[5.39939,0.0343564,-1.79426e-05,3.54497e-09,-2.49293e-13,-68358.4,8.30096], Tmin=(798.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-567.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CdCFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OF(140)',
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
    label = 'C=C(F)OC(=C)F(1563)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-400.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,390,483,597,572,732,631,807,1275,1439,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,1014.38],'cm^-1')),
        HinderedRotor(inertia=(0.188897,'amu*angstrom^2'), symmetry=1, barrier=(4.34311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187015,'amu*angstrom^2'), symmetry=1, barrier=(4.29985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68833,0.0563022,-7.66948e-05,6.89057e-08,-2.62612e-11,-48102.8,21.4617], Tmin=(100,'K'), Tmax=(745.408,'K')), NASAPolynomial(coeffs=[4.74728,0.0332118,-1.67964e-05,3.32045e-09,-2.35401e-13,-48373.3,8.84683], Tmin=(745.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-400.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CdCFO) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=C(F)CC(O)(F)CF(11596)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {14,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (-1071.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.278181,0.0889902,-0.000120053,8.98647e-08,-2.75953e-11,-128753,27.9999], Tmin=(100,'K'), Tmax=(789.233,'K')), NASAPolynomial(coeffs=[10.8858,0.0352273,-1.78697e-05,3.54866e-09,-2.52912e-13,-130428,-20.6648], Tmin=(789.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1071.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO)"""),
)

species(
    label = '[CH2]C(F)OC([O])(F)CF(11774)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
9  C u1 p0 c0 {7,S} {13,S} {14,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-636.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,487,638,688,1119,1325,1387,3149,528,1116,1182,1331,1402,1494,3075,3110,3000,3100,440,815,1455,1000,235.899,235.899,235.899,235.899],'cm^-1')),
        HinderedRotor(inertia=(0.241896,'amu*angstrom^2'), symmetry=1, barrier=(9.55231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879451,'amu*angstrom^2'), symmetry=1, barrier=(34.7292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241895,'amu*angstrom^2'), symmetry=1, barrier=(9.55232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879455,'amu*angstrom^2'), symmetry=1, barrier=(34.7292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00736788,0.0957316,-0.000135922,1.06357e-07,-3.39869e-11,-76360.8,30.2393], Tmin=(100,'K'), Tmax=(760.645,'K')), NASAPolynomial(coeffs=[11.2946,0.0363753,-1.88702e-05,3.76709e-09,-2.68693e-13,-78077.9,-21.127], Tmin=(760.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-636.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = '[CH2]C(F)OC(O)(F)[CH]F(11775)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-693.699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.659932,0.111508,-0.000181228,1.54391e-07,-5.19852e-11,-83273.5,32.3587], Tmin=(100,'K'), Tmax=(793.627,'K')), NASAPolynomial(coeffs=[13.23,0.0331355,-1.72903e-05,3.39845e-09,-2.37815e-13,-85214.8,-29.7823], Tmin=(793.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-693.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = 'C[C](F)OC([O])(F)CF(11776)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-654.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,202.711,205.468,205.952,209.707],'cm^-1')),
        HinderedRotor(inertia=(0.274672,'amu*angstrom^2'), symmetry=1, barrier=(8.00623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268993,'amu*angstrom^2'), symmetry=1, barrier=(8.01382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15098,'amu*angstrom^2'), symmetry=1, barrier=(34.3661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626131,'amu*angstrom^2'), symmetry=1, barrier=(19.3103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0662495,0.0976171,-0.000148757,1.26849e-07,-4.39446e-11,-78544.6,29.9063], Tmin=(100,'K'), Tmax=(751.507,'K')), NASAPolynomial(coeffs=[10.3476,0.03751,-1.9446e-05,3.85339e-09,-2.72487e-13,-79977.7,-16.4804], Tmin=(751.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-654.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs))"""),
)

species(
    label = 'C[C](F)OC(O)(F)[CH]F(11777)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u1 p0 c0 {3,S} {6,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-711.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,334,575,1197,1424,3202,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612399,0.111756,-0.000187213,1.6405e-07,-5.62572e-11,-85462.5,31.6036], Tmin=(100,'K'), Tmax=(816.631,'K')), NASAPolynomial(coeffs=[12.069,0.0346658,-1.8108e-05,3.54437e-09,-2.4671e-13,-87034.4,-23.9509], Tmin=(816.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
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
    label = 'C=C(F)O[C](O)CF(11778)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {6,S} {13,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-554.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,360,370,350,326,540,652,719,1357,2950,3100,1380,975,1025,1650,373.536,373.541,373.544,373.545],'cm^-1')),
        HinderedRotor(inertia=(0.134793,'amu*angstrom^2'), symmetry=1, barrier=(13.3458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134761,'amu*angstrom^2'), symmetry=1, barrier=(13.3455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134782,'amu*angstrom^2'), symmetry=1, barrier=(13.3457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20816,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524442,0.0766457,-8.47136e-05,4.69295e-08,-1.03165e-11,-66569.8,27.3596], Tmin=(100,'K'), Tmax=(1103.01,'K')), NASAPolynomial(coeffs=[15.5758,0.022063,-1.04858e-05,2.06572e-09,-1.48056e-13,-69890.2,-46.73], Tmin=(1103.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(CdCFO) + group(Cds-CdsHH) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C(O)(F)OC(=C)F(11779)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {5,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-590.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,326,540,652,719,1357,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,216.937,221.8,222.266],'cm^-1')),
        HinderedRotor(inertia=(0.243477,'amu*angstrom^2'), symmetry=1, barrier=(8.83924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.621201,'amu*angstrom^2'), symmetry=1, barrier=(21.2818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.606939,'amu*angstrom^2'), symmetry=1, barrier=(21.285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668281,'amu*angstrom^2'), symmetry=1, barrier=(21.2856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.166692,0.0870555,-0.000109006,6.85299e-08,-1.7015e-11,-70861.5,25.2266], Tmin=(100,'K'), Tmax=(983.711,'K')), NASAPolynomial(coeffs=[15.9523,0.0228674,-1.11287e-05,2.19791e-09,-1.57291e-13,-73967.2,-50.6703], Tmin=(983.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-590.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cs-F1sO2sO2s)(H)(H))"""),
)

species(
    label = 'C=[C]OC(O)(F)CF(11780)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {8,S}
4  O u0 p2 c0 {5,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {8,D} {11,S} {12,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-554.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.947465,'amu*angstrom^2'), symmetry=1, barrier=(21.7841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94869,'amu*angstrom^2'), symmetry=1, barrier=(21.8123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.947575,'amu*angstrom^2'), symmetry=1, barrier=(21.7866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949147,'amu*angstrom^2'), symmetry=1, barrier=(21.8228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.000774826,0.0825284,-9.34812e-05,5.13645e-08,-1.09121e-11,-66492.8,28.391], Tmin=(100,'K'), Tmax=(1158.87,'K')), NASAPolynomial(coeffs=[19.5353,0.0151026,-6.20794e-06,1.1587e-09,-8.1373e-14,-71020.4,-68.732], Tmin=(1158.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'CH2CFO(79)',
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
    label = 'O[C](F)CF(9946)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78656,0.0175844,0.000121424,-4.80809e-07,5.06984e-10,-56694.5,10.4186], Tmin=(10,'K'), Tmax=(351.977,'K')), NASAPolynomial(coeffs=[6.00134,0.0208372,-1.35648e-05,4.29025e-09,-5.21054e-13,-57026.5,-0.454778], Tmin=(351.977,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-471.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""O[C](F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CF(70)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (96.6638,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(45.014,'amu')),
        NonlinearRotor(inertia=([4.28106,49.5096,53.7907],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([443.335,623.171,819.82,951.292,1166.86,1400.31,1729.5,3121.37,3250.42],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0356,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06008,-0.00578847,7.22738e-05,-1.3105e-07,7.80971e-11,11626.8,7.18025], Tmin=(10,'K'), Tmax=(523.261,'K')), NASAPolynomial(coeffs=[2.54,0.0141591,-8.78068e-06,2.63244e-09,-3.03932e-13,11671.9,12.4399], Tmin=(523.261,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(96.6638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)(F)CF(9948)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6 C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-628.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,180],'cm^-1')),
        HinderedRotor(inertia=(1.00515,'amu*angstrom^2'), symmetry=1, barrier=(23.1104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00489,'amu*angstrom^2'), symmetry=1, barrier=(23.1043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78076,0.0150815,0.000141095,-3.93332e-07,3.04195e-10,-75641.9,11.4432], Tmin=(10,'K'), Tmax=(466.103,'K')), NASAPolynomial(coeffs=[6.42151,0.0289801,-2.12929e-05,7.21919e-09,-9.08945e-13,-76285.2,-3.54139], Tmin=(466.103,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-628.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(O)(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(F)O[C](F)CF(10005)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-576.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.975331,'amu*angstrom^2'), symmetry=1, barrier=(22.4248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975856,'amu*angstrom^2'), symmetry=1, barrier=(22.4368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975442,'amu*angstrom^2'), symmetry=1, barrier=(22.4273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767538,0.0734787,-8.38888e-05,4.87544e-08,-1.13481e-11,-69166,23.987], Tmin=(100,'K'), Tmax=(1038.81,'K')), NASAPolynomial(coeffs=[13.7913,0.0233297,-1.14754e-05,2.28226e-09,-1.6408e-13,-71871.8,-39.341], Tmin=(1038.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-576.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(CsCsF1sO2s)"""),
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
    label = 'C=C(F)OC([O])(F)CF(10014)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-720.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02228,'amu*angstrom^2'), symmetry=1, barrier=(23.5043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02339,'amu*angstrom^2'), symmetry=1, barrier=(23.5297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02276,'amu*angstrom^2'), symmetry=1, barrier=(23.5152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305121,0.0840394,-9.76465e-05,5.73228e-08,-1.34323e-11,-86509.1,26.407], Tmin=(100,'K'), Tmax=(1033.73,'K')), NASAPolynomial(coeffs=[15.4342,0.0254973,-1.26981e-05,2.53808e-09,-1.82958e-13,-89636.9,-47.0838], Tmin=(1033.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(O2sj(Cs-F1sO2sCs))"""),
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
    label = 'C=C(F)O[C](O)F(3813)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {6,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-554.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,482,586,761,1411,2950,3100,1380,975,1025,1650,180,192.149,1222.94],'cm^-1')),
        HinderedRotor(inertia=(1.11312,'amu*angstrom^2'), symmetry=1, barrier=(25.653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.548691,'amu*angstrom^2'), symmetry=1, barrier=(14.7073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00279,'amu*angstrom^2'), symmetry=1, barrier=(23.0589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13823,0.0579512,-6.1194e-05,3.12925e-08,-6.21325e-12,-66616.2,23.1484], Tmin=(100,'K'), Tmax=(1233.96,'K')), NASAPolynomial(coeffs=[15.3858,0.0117659,-5.05101e-06,9.60276e-10,-6.79169e-14,-70132.4,-48.5833], Tmin=(1233.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFHOO) + group(CdCFO) + group(Cds-CdsHH) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'C=C(F)OC(O)(F)[CH]F(11781)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {2,S} {6,S} {10,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-778.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,326,540,652,719,1357,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,240.31,240.45,240.615,240.673],'cm^-1')),
        HinderedRotor(inertia=(0.326858,'amu*angstrom^2'), symmetry=1, barrier=(13.4306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.326981,'amu*angstrom^2'), symmetry=1, barrier=(13.4283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713184,'amu*angstrom^2'), symmetry=1, barrier=(29.2833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500916,'amu*angstrom^2'), symmetry=1, barrier=(20.6108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.23372,0.0979995,-0.00013515,9.31176e-08,-2.52376e-11,-93427.1,28.0851], Tmin=(100,'K'), Tmax=(904.953,'K')), NASAPolynomial(coeffs=[16.575,0.0236995,-1.19894e-05,2.38235e-09,-1.70194e-13,-96469.1,-51.3278], Tmin=(904.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-778.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = '[CH]=C(F)OC(O)(F)CF(11107)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-716.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,293,496,537,1218,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.975846,'amu*angstrom^2'), symmetry=1, barrier=(22.4366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97587,'amu*angstrom^2'), symmetry=1, barrier=(22.4372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976461,'amu*angstrom^2'), symmetry=1, barrier=(22.4508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976468,'amu*angstrom^2'), symmetry=1, barrier=(22.4509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.269216,0.0960023,-0.000123457,7.75047e-08,-1.90079e-11,-86033.5,27.7629], Tmin=(100,'K'), Tmax=(1000.2,'K')), NASAPolynomial(coeffs=[18.6055,0.0205182,-1.02525e-05,2.04992e-09,-1.47834e-13,-89809.1,-63.3002], Tmin=(1000.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'CC(=O)F(214)',
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
    label = 'OC(F)(CF)[O+]=[C-]CF(11782)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p1 c+1 {6,S} {9,D}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
9  C u0 p1 c-1 {4,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-614.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.409528,0.108576,-0.000181764,1.64834e-07,-5.859e-11,-73773.5,40.3569], Tmin=(100,'K'), Tmax=(816.455,'K')), NASAPolynomial(coeffs=[9.29372,0.041195,-2.15147e-05,4.21758e-09,-2.94123e-13,-74696.6,-0.437788], Tmin=(816.455,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(F)OC(O)(F)CF(11783)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u0 p1 c0 {8,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-656.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.309571,0.103337,-0.000156141,1.26582e-07,-4.1407e-11,-78788.4,28.9995], Tmin=(100,'K'), Tmax=(746.583,'K')), NASAPolynomial(coeffs=[12.4695,0.034866,-1.85658e-05,3.72673e-09,-2.65859e-13,-80696.5,-28.9169], Tmin=(746.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=C(F)OC(=O)CF(11784)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {4,D} {5,S}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-699.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,326,540,652,719,1357,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35102,0.0645619,-7.44396e-05,5.06037e-08,-1.48493e-11,-83997.4,23.2434], Tmin=(100,'K'), Tmax=(806.522,'K')), NASAPolynomial(coeffs=[7.52541,0.033938,-1.74813e-05,3.51992e-09,-2.53822e-13,-84993.3,-5.21645], Tmin=(806.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(CdCFO) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)OC(=C)F(6001)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {4,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-372.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,326,540,652,719,1357,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,303.863,304.18,304.195],'cm^-1')),
        HinderedRotor(inertia=(0.12664,'amu*angstrom^2'), symmetry=1, barrier=(8.30843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126876,'amu*angstrom^2'), symmetry=1, barrier=(8.30907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30761,'amu*angstrom^2'), symmetry=1, barrier=(20.2121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3673.09,'J/mol'), sigma=(5.85032,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.73 K, Pc=41.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15728,0.0661994,-8.01239e-05,5.31159e-08,-1.43565e-11,-44705.6,22.452], Tmin=(100,'K'), Tmax=(895.094,'K')), NASAPolynomial(coeffs=[10.2684,0.0254826,-1.18891e-05,2.29347e-09,-1.61493e-13,-46336.7,-20.494], Tmin=(895.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-372.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(F)OC(O)=CF(7319)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {5,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {1,S} {3,S} {8,D}
7  C u0 p0 c0 {2,S} {5,D} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-537.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,326,540,652,719,1357,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,180,180,2174.14],'cm^-1')),
        HinderedRotor(inertia=(0.361808,'amu*angstrom^2'), symmetry=1, barrier=(8.31869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3628,'amu*angstrom^2'), symmetry=1, barrier=(8.34149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65311,'amu*angstrom^2'), symmetry=1, barrier=(38.0083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684257,0.0807532,-0.000128968,1.1458e-07,-4.04298e-11,-64537.6,25.5969], Tmin=(100,'K'), Tmax=(805.519,'K')), NASAPolynomial(coeffs=[8.1229,0.0320341,-1.63085e-05,3.18375e-09,-2.222e-13,-65353.8,-6.30908], Tmin=(805.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'C#COC(O)(F)CF(9988)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-584.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25682,'amu*angstrom^2'), symmetry=1, barrier=(28.8967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25487,'amu*angstrom^2'), symmetry=1, barrier=(28.8519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25608,'amu*angstrom^2'), symmetry=1, barrier=(28.8797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25638,'amu*angstrom^2'), symmetry=1, barrier=(28.8867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3718.53,'J/mol'), sigma=(6.08719,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.83 K, Pc=37.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.264929,0.087075,-0.000105067,5.94158e-08,-1.28119e-11,-70140.5,23.6722], Tmin=(100,'K'), Tmax=(1149.89,'K')), NASAPolynomial(coeffs=[22.0449,0.00946913,-3.83372e-06,7.24803e-10,-5.19385e-14,-75271.3,-87.076], Tmin=(1149.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'OC(F)(CF)[O+]=[C-]F(11116)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p1 c+1 {6,S} {8,D}
5  O u0 p2 c0 {6,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p1 c-1 {3,S} {4,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-846.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943917,0.0723286,-9.14671e-05,6.01624e-08,-1.59924e-11,-101735,22.4604], Tmin=(100,'K'), Tmax=(910.563,'K')), NASAPolynomial(coeffs=[11.8098,0.0245969,-1.28388e-05,2.59612e-09,-1.87625e-13,-103714,-28.9436], Tmin=(910.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-846.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CJ2_singlet-FO)"""),
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
    E0 = (-443.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-145.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-189.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-1.56564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-271.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-271.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-227.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (41.3023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-178.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (115.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-451.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-202.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-259.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-220.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-278.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-72.8404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-108.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-72.3964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-319.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-123.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-138.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-99.7596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-188.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-157.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-95.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-455.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-394.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-125.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-233.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-401.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (27.4674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-219.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-158.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-7.01775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['HF(38)', 'CH2CO(28)', 'O=C(F)CF(476)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.29418e+10,'s^-1'), n=0.745914, Ea=(126.868,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H2(8)', 'C=C(F)OC(O)(F)[C]F(11771)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'C=C(F)OC(O)F(4123)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.79957e-05,'m^3/(mol*s)'), n=3.10993, Ea=(23.7398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(S)(25)', 'C=C(F)OC(O)(F)F(3586)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(163.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['C=C(F)OCC(O)(F)F(11772)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['C=C(F)OC(F)(F)CO(11773)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC(F)=CF(3655)', 'C=C(O)F(213)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000170623,'m^3/(mol*s)'), n=2.84144, Ea=(220.872,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OR] for rate rule [Cd/monosub_Cd/disub;H_OCd]
Euclidian distance = 1.4142135623730951
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C(F)OF(479)', 'C=C(O)F(213)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.000342053,'m^3/(mol*s)'), n=2.805, Ea=(189.117,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OR]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', 'C=C(F)OC(F)=CF(10044)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OF(140)', 'C=C(F)OC(=C)F(1563)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2500,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['O=C(F)CC(O)(F)CF(11596)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(118.275,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(F)OC([O])(F)CF(11774)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(F)OC(O)(F)[CH]F(11775)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[C](F)OC([O])(F)CF(11776)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[C](F)OC(O)(F)[CH]F(11777)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'C=C(F)O[C](O)CF(11778)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', '[CH2]C(O)(F)OC(=C)F(11779)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'C=[C]OC(O)(F)CF(11780)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2CFO(79)', 'O[C](F)CF(9946)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(4.65139,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2CF(70)', '[O]C(O)(F)CF(9948)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['OH(7)', 'C=C(F)O[C](F)CF(10005)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'C=C(F)OC([O])(F)CF(10014)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2F(46)', 'C=C(F)O[C](O)F(3813)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', 'C=C(F)OC(O)(F)[CH]F(11781)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -1.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(5)', '[CH]=C(F)OC(O)(F)CF(11107)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -12.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['O=C(F)CF(476)', 'CC(=O)F(214)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(115.214,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(F)OC(O)(F)CF(11594)'],
    products = ['OC(F)=CF(3655)', 'CC(=O)F(214)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.66666e+07,'s^-1'), n=1.2, Ea=(175.241,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OC(F)(CF)[O+]=[C-]CF(11782)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(80.2203,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(F)OC(O)(F)CF(11783)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(13.9191,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

reaction(
    label = 'reaction30',
    reactants = ['HF(38)', 'C=C(F)OC(=O)CF(11784)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(169.838,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F2(78)', 'C=C(O)OC(=C)F(6001)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HF(38)', 'C=C(F)OC(O)=CF(7319)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(190.141,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction33',
    reactants = ['HF(38)', 'C#COC(O)(F)CF(9988)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.56721e+32,'m^3/(mol*s)'), n=-7.40875, Ea=(297.855,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_3COCdCddCtO2d->Ct',), comment="""Estimated from node HF_3COCdCddCtO2d->Ct"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(S)(25)', 'OC(F)(CF)[O+]=[C-]F(11116)'],
    products = ['C=C(F)OC(O)(F)CF(11594)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.1e+24,'cm^3/(mol*s)'), n=-3.8, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for CF
Exact match found for rate rule [CF]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #3799',
    isomers = [
        'C=C(F)OC(O)(F)CF(11594)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', 'O=C(F)CF(476)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3799',
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

