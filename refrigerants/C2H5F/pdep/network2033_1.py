species(
    label = 'COC(C)(F)OF(5723)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-535.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,375,361,584,565,722,1474,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,264.374,264.388,264.413],'cm^-1')),
        HinderedRotor(inertia=(0.168907,'amu*angstrom^2'), symmetry=1, barrier=(8.38139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168936,'amu*angstrom^2'), symmetry=1, barrier=(8.38132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543819,'amu*angstrom^2'), symmetry=1, barrier=(26.98,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168933,'amu*angstrom^2'), symmetry=1, barrier=(8.38161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953409,0.0730048,-9.39965e-05,7.04885e-08,-2.21274e-11,-64299.2,21.1915], Tmin=(100,'K'), Tmax=(767.371,'K')), NASAPolynomial(coeffs=[8.54527,0.0334315,-1.6642e-05,3.28567e-09,-2.33627e-13,-65464.3,-13.4247], Tmin=(767.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-535.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCFOO) + group(Cs-CsHHH) + group(Cs-OsHHH)"""),
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
    label = 'CH2O(20)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'COC(F)OF(2851)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-490.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,509,613,660,1171,1360,1414,3084,2750,2800,2850,1350,1500,750,1050,1375,1000,180,188.444,189.798],'cm^-1')),
        HinderedRotor(inertia=(0.274162,'amu*angstrom^2'), symmetry=1, barrier=(6.53147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259593,'amu*angstrom^2'), symmetry=1, barrier=(6.46414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0175,'amu*angstrom^2'), symmetry=1, barrier=(24.4839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0488,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2920.12,'J/mol'), sigma=(5.00118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=456.12 K, Pc=52.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57502,0.0436326,-9.01374e-05,1.80057e-07,-1.52567e-10,-58994.4,11.3996], Tmin=(10,'K'), Tmax=(368.022,'K')), NASAPolynomial(coeffs=[3.85145,0.0338158,-2.23602e-05,6.98187e-09,-8.28973e-13,-58968.7,10.9691], Tmin=(368.022,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-490.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""COC(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC(O)(F)OF(5746)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-575.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,375,361,584,565,722,1474,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.74446,'amu*angstrom^2'), symmetry=1, barrier=(17.1166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.74329,'amu*angstrom^2'), symmetry=1, barrier=(17.0897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742741,'amu*angstrom^2'), symmetry=1, barrier=(17.0771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0488,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77615,0.0149641,0.00015787,-4.21306e-07,3.15019e-10,-69155.3,10.4668], Tmin=(10,'K'), Tmax=(480.216,'K')), NASAPolynomial(coeffs=[6.67586,0.0318421,-2.30161e-05,7.74e-09,-9.71547e-13,-69906.9,-6.32167], Tmin=(480.216,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-575.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""CC(O)(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3OH(26)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-211.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1964.49,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0117392,'amu*angstrom^2'), symmetry=1, barrier=(6.09573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.0419,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.93,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.65851,-0.0162983,6.91938e-05,-7.58373e-08,2.80428e-11,-25612,-0.897331], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.52727,0.0103179,-3.62893e-06,5.77448e-10,-3.42183e-14,-26002.9,5.16759], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-211.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3OH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C(F)OF(1352)',
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
    label = 'C=C(F)OC(979)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-334.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.213121,'amu*angstrom^2'), symmetry=1, barrier=(4.90007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213355,'amu*angstrom^2'), symmetry=1, barrier=(4.90544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0696,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87192,0.00971488,0.000143084,-4.18647e-07,3.70296e-10,-40178.9,9.38863], Tmin=(10,'K'), Tmax=(382.222,'K')), NASAPolynomial(coeffs=[3.94797,0.0298295,-1.79163e-05,5.29829e-09,-6.12521e-13,-40337.5,7.09678], Tmin=(382.222,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-334.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(F)OC""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CO[C](C)OF(2981)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {1,S} {6,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-96.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,244.949,244.958,4000],'cm^-1')),
        HinderedRotor(inertia=(0.265178,'amu*angstrom^2'), symmetry=1, barrier=(11.2914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26519,'amu*angstrom^2'), symmetry=1, barrier=(11.2914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00918011,'amu*angstrom^2'), symmetry=1, barrier=(28.4332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667787,'amu*angstrom^2'), symmetry=1, barrier=(28.4332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0769,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72197,0.0548489,-5.91385e-05,3.88787e-08,-1.1165e-11,-11504.6,21.3678], Tmin=(100,'K'), Tmax=(820.067,'K')), NASAPolynomial(coeffs=[6.67712,0.0306799,-1.49314e-05,2.94167e-09,-2.09751e-13,-12317.3,-1.55506], Tmin=(820.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(Cs_P)"""),
)

species(
    label = 'COC(C)([O])F(5747)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-411.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,208.756,208.786],'cm^-1')),
        HinderedRotor(inertia=(0.194181,'amu*angstrom^2'), symmetry=1, barrier=(6.00719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194435,'amu*angstrom^2'), symmetry=1, barrier=(6.00761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194171,'amu*angstrom^2'), symmetry=1, barrier=(6.00745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0769,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23867,0.0467582,-1.84316e-05,-5.41098e-08,6.34479e-11,-49414.6,17.3512], Tmin=(100,'K'), Tmax=(463.383,'K')), NASAPolynomial(coeffs=[5.35142,0.0329407,-1.59549e-05,3.11383e-09,-2.19898e-13,-49843.2,3.21624], Tmin=(463.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = 'CH3O(27)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (10.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.15,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71181,-0.00280463,3.76551e-05,-4.73072e-08,1.86588e-11,1295.7,6.57241], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.75779,0.00744142,-2.69705e-06,4.38091e-10,-2.63537e-14,378.112,-1.9668], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(10.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH3O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'CC([O])(F)OF(3850)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6 C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-316.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,373,336,502,477,708,1220,1477,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.529573,'amu*angstrom^2'), symmetry=1, barrier=(12.1759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529293,'amu*angstrom^2'), symmetry=1, barrier=(12.1695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.72704,0.0223816,9.90217e-05,-3.51305e-07,3.23898e-10,-38024.3,11.0719], Tmin=(10,'K'), Tmax=(396.547,'K')), NASAPolynomial(coeffs=[5.49064,0.0292883,-2.05211e-05,6.71941e-09,-8.27658e-13,-38358.4,1.7465], Tmin=(396.547,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-316.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CC([O])(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CO[C](C)F(5111)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.207719,'amu*angstrom^2'), symmetry=1, barrier=(4.77587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207629,'amu*angstrom^2'), symmetry=1, barrier=(4.7738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207635,'amu*angstrom^2'), symmetry=1, barrier=(4.77393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0775,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66954,0.0320039,-1.60696e-05,2.38872e-09,3.29384e-13,-34298.8,10.4555], Tmin=(10,'K'), Tmax=(1317.55,'K')), NASAPolynomial(coeffs=[8.63575,0.0219387,-1.03165e-05,2.36478e-09,-2.13886e-13,-36042.5,-16.524], Tmin=(1317.55,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-285.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""CO[C](C)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CO[C](F)OF(2930)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {6,S}
5 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-283.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2800,2850,1350,1500,750,1050,1375,1000,482,586,761,1411,180,180,2231.05],'cm^-1')),
        HinderedRotor(inertia=(0.0762271,'amu*angstrom^2'), symmetry=1, barrier=(1.75261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750144,'amu*angstrom^2'), symmetry=1, barrier=(1.72473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02126,'amu*angstrom^2'), symmetry=1, barrier=(46.4729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4686,0.0377789,-3.60845e-05,2.08158e-08,-5.42975e-12,-33987.1,19.3424], Tmin=(100,'K'), Tmax=(876.652,'K')), NASAPolynomial(coeffs=[5.57455,0.0236074,-1.18371e-05,2.37703e-09,-1.71625e-13,-34531.7,4.76682], Tmin=(876.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsF1sO2sO2s)"""),
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
    label = '[CH2]C(F)(OC)OF(5748)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-320.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,364,474,579,619,655,1345,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.191713,'amu*angstrom^2'), symmetry=1, barrier=(4.40787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542116,'amu*angstrom^2'), symmetry=1, barrier=(12.4643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189608,'amu*angstrom^2'), symmetry=1, barrier=(4.35946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65587,'amu*angstrom^2'), symmetry=1, barrier=(38.0717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.593783,0.0822003,-0.000129965,1.1311e-07,-3.9304e-11,-38480.9,22.914], Tmin=(100,'K'), Tmax=(792.282,'K')), NASAPolynomial(coeffs=[9.11704,0.0304865,-1.56193e-05,3.06154e-09,-2.1439e-13,-39559,-14.5015], Tmin=(792.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCFOO) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(Csj(Cs-F1sO2sO2s)(H)(H))"""),
)

species(
    label = '[CH2]OC(C)(F)OF(5749)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {11,S} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-347.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,375,361,584,565,722,1474,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,272.127,274.783,275.193],'cm^-1')),
        HinderedRotor(inertia=(0.103081,'amu*angstrom^2'), symmetry=1, barrier=(5.51127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308385,'amu*angstrom^2'), symmetry=1, barrier=(16.5739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103988,'amu*angstrom^2'), symmetry=1, barrier=(5.50052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578657,'amu*angstrom^2'), symmetry=1, barrier=(30.5881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881989,0.074059,-9.73992e-05,6.93205e-08,-2.00795e-11,-41700.3,22.2395], Tmin=(100,'K'), Tmax=(837.111,'K')), NASAPolynomial(coeffs=[10.6485,0.0273939,-1.37857e-05,2.73525e-09,-1.9515e-13,-43335.5,-23.1423], Tmin=(837.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCFOO) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CsJOCH3)"""),
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
    label = 'COC(C)=O(297)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-425.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,871.044,871.091,871.16,871.435,874.258],'cm^-1')),
        HinderedRotor(inertia=(0.103236,'amu*angstrom^2'), symmetry=1, barrier=(2.37359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104136,'amu*angstrom^2'), symmetry=1, barrier=(2.3943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105986,'amu*angstrom^2'), symmetry=1, barrier=(2.43683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75061,0.0220814,1.21738e-05,-2.1921e-08,7.33604e-12,-51179.4,16.5003], Tmin=(100,'K'), Tmax=(1154.02,'K')), NASAPolynomial(coeffs=[5.62534,0.0261099,-1.125e-05,2.11764e-09,-1.47682e-13,-52774.6,-1.8174], Tmin=(1154.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)HHH) + group(Cs-OsHHH) + group(Cds-OdCsOs)"""),
)

species(
    label = 'C=C(OC)OF(645)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {5,S}
4  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-147.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,724.02,724.021],'cm^-1')),
        HinderedRotor(inertia=(0.543379,'amu*angstrom^2'), symmetry=1, barrier=(12.4934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543387,'amu*angstrom^2'), symmetry=1, barrier=(12.4935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0759146,'amu*angstrom^2'), symmetry=1, barrier=(12.4936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3075.29,'J/mol'), sigma=(5.21583,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.35 K, Pc=49.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56001,0.0555296,-6.08095e-05,3.59409e-08,-8.65624e-12,-17628.5,17.8726], Tmin=(100,'K'), Tmax=(999.394,'K')), NASAPolynomial(coeffs=[10.0651,0.0214887,-9.71767e-06,1.85937e-09,-1.30759e-13,-19328.5,-23.1546], Tmin=(999.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
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
    E0 = (-370.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (65.4823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-19.0353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-107.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-89.5233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (113.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-201.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-20.4385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-43.0573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-45.5948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-9.92739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (27.7939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1.09474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-199.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-104.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['COC(C)(F)OF(5723)'],
    products = ['HF(38)', 'CH2O(20)', 'CC(=O)F(499)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(27.897,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', 'COC(F)OF(2851)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.67081e+07,'m^3/(mol*s)'), n=-0.25099, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05783420702249951, var=0.11014630093021238, Tref=1000.0, N=12, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', 'CC(O)(F)OF(5746)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.67448e-06,'m^3/(mol*s)'), n=3.15288, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node OH_2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3OH(26)', 'C=C(F)OF(1352)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(93.6,'cm^3/(mol*s)'), n=2.85, Ea=(175.31,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OCmethyl]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(169)', 'C=C(F)OC(979)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1250,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'CO[C](C)OF(2981)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'COC(C)([O])F(5747)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3O(27)', 'C[C](F)OF(5141)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(19)', 'CC([O])(F)OF(3850)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21e+07,'m^3/(mol*s)'), n=-8.5357e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C_2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]F(160)', 'CO[C](C)F(5111)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH3(19)', 'CO[C](F)OF(2930)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.29446e+18,'m^3/(mol*s)'), n=-4.19701, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R
Ea raised from -10.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[CH2]C(F)(OC)OF(5748)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.29291e+24,'m^3/(mol*s)'), n=-6.27825, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_Sp-3C-2C_Ext-3C-R_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_Sp-3C-2C_Ext-3C-R_Ext-3C-R
Ea raised from -4.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[CH2]OC(C)(F)OF(5749)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.31504e+12,'m^3/(mol*s)'), n=-2.20453, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F2(77)', 'COC(C)=O(297)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(98.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'C=C(OC)OF(645)'],
    products = ['COC(C)(F)OF(5723)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(186.536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

network(
    label = 'PDepNetwork #2033',
    isomers = [
        'COC(C)(F)OF(5723)',
    ],
    reactants = [
        ('HF(38)', 'CH2O(20)', 'CC(=O)F(499)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2033',
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

