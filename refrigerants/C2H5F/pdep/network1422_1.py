species(
    label = '[CH2]C(F)OOC(=C)F(3748)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
6  C u1 p0 c0 {5,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-276.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.79733,'amu*angstrom^2'), symmetry=1, barrier=(41.3242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240024,'amu*angstrom^2'), symmetry=1, barrier=(5.51862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79884,'amu*angstrom^2'), symmetry=1, barrier=(41.3588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239316,'amu*angstrom^2'), symmetry=1, barrier=(5.50234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986101,0.0734595,-8.71765e-05,5.9752e-08,-1.74801e-11,-33180,27.7661], Tmin=(100,'K'), Tmax=(812.366,'K')), NASAPolynomial(coeffs=[8.47902,0.0365668,-1.90585e-05,3.85353e-09,-2.78423e-13,-34397.5,-6.82616], Tmin=(812.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CJCOOH)"""),
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
    label = 'FC1CO1(1534)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-280.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,2750,3150,900,1100,1311.49,1311.83,1311.95,1311.99],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2921.61,'J/mol'), sigma=(5.00203,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=456.35 K, Pc=52.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.11059,-0.00875873,9.68406e-05,-1.46154e-07,7.12277e-11,-33718.7,8.8538], Tmin=(10,'K'), Tmax=(645.87,'K')), NASAPolynomial(coeffs=[0.921679,0.0253893,-1.59063e-05,4.73976e-09,-5.4007e-13,-33607.1,20.5195], Tmin=(645.87,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-280.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1CO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)OCC(=O)F(4389)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-606.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02468,0.070604,-7.6504e-05,4.49692e-08,-1.10113e-11,-72898,26.6844], Tmin=(100,'K'), Tmax=(969.967,'K')), NASAPolynomial(coeffs=[10.5723,0.0312303,-1.56134e-05,3.11771e-09,-2.24239e-13,-74750.1,-19.0859], Tmin=(969.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(COCsFO) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(F)OO[CH]F(2054)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {8,S} {9,S}
7  C u1 p0 c0 {2,S} {4,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-242.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,2950,3100,1380,975,1025,1650,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.281406,'amu*angstrom^2'), symmetry=1, barrier=(6.47007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279759,'amu*angstrom^2'), symmetry=1, barrier=(6.43222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20474,'amu*angstrom^2'), symmetry=1, barrier=(50.6912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40291,0.0648407,-0.000106904,1.01036e-07,-3.78572e-11,-29105.1,24.0868], Tmin=(100,'K'), Tmax=(794.299,'K')), NASAPolynomial(coeffs=[5.47364,0.0306094,-1.63276e-05,3.24954e-09,-2.29337e-13,-29318.6,8.11207], Tmin=(794.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + group(Cds-CdsHH) + radical(CsF1sHO2s)"""),
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
    label = '[CH]C(F)OOC(=C)F(4390)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {10,S} {11,S}
8  C u2 p0 c0 {5,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-42.4717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180,1249.23,1249.31],'cm^-1')),
        HinderedRotor(inertia=(0.122291,'amu*angstrom^2'), symmetry=1, barrier=(2.8117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12237,'amu*angstrom^2'), symmetry=1, barrier=(2.81352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01545,'amu*angstrom^2'), symmetry=1, barrier=(46.3392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01549,'amu*angstrom^2'), symmetry=1, barrier=(46.3401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935378,0.0746026,-9.66457e-05,7.13805e-08,-2.21704e-11,-5004.28,26.9533], Tmin=(100,'K'), Tmax=(772.545,'K')), NASAPolynomial(coeffs=[8.73579,0.0342125,-1.82192e-05,3.69917e-09,-2.67313e-13,-6209.46,-8.66577], Tmin=(772.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.4717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'F[C]1CCC(F)OO1(4391)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {13,S}
8  C u1 p0 c0 {2,S} {4,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-413.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23889,0.0460367,1.0889e-05,-5.60008e-08,2.7987e-11,-49597.8,19.6625], Tmin=(100,'K'), Tmax=(914.665,'K')), NASAPolynomial(coeffs=[16.8759,0.0161883,-3.357e-06,4.43642e-10,-3.02278e-14,-54070.3,-63.1941], Tmin=(914.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCFHO) + ring(12dioxane) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH2]C1(F)CC(F)OO1(4392)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {5,S} {11,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-416.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898297,0.0528003,-6.50403e-07,-4.93653e-08,2.68098e-11,-49959.9,21.1602], Tmin=(100,'K'), Tmax=(914.594,'K')), NASAPolynomial(coeffs=[19.2793,0.01303,-2.04275e-06,2.09272e-10,-1.47471e-14,-55021,-75.164], Tmin=(914.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CJCOOH)"""),
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
    label = 'C=COOC(=C)F(1690)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {2,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-122.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,326,540,652,719,1357,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.699439,'amu*angstrom^2'), symmetry=1, barrier=(16.0815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242363,'amu*angstrom^2'), symmetry=1, barrier=(5.5724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40117,'amu*angstrom^2'), symmetry=1, barrier=(32.2157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6266,0.0543692,-4.66295e-05,2.09174e-08,-3.92893e-12,-14597,23.9587], Tmin=(100,'K'), Tmax=(1228.56,'K')), NASAPolynomial(coeffs=[10.1773,0.0265294,-1.26389e-05,2.47271e-09,-1.75625e-13,-16698,-19.0537], Tmin=(1228.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(F)O[O](1327)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-96.7758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.725193,'amu*angstrom^2'), symmetry=1, barrier=(16.6736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3364.7,'J/mol'), sigma=(5.40718,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.56 K, Pc=48.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91162,0.00529573,8.44313e-05,-1.79248e-07,1.11122e-10,-11633.8,10.2654], Tmin=(10,'K'), Tmax=(555.797,'K')), NASAPolynomial(coeffs=[4.21889,0.0228301,-1.61816e-05,5.35569e-09,-6.6581e-13,-11972.9,6.21972], Tmin=(555.797,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-96.7758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CDC(F)O[O]""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CHF(55)',
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
    label = 'C=C(F)OOC(=C)F(1464)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-307.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,262,390,483,597,572,732,631,807,1275,1439,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.121995,'amu*angstrom^2'), symmetry=1, barrier=(2.80491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125049,'amu*angstrom^2'), symmetry=1, barrier=(2.87513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385701,'amu*angstrom^2'), symmetry=1, barrier=(8.86802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11558,0.0719182,-0.000112429,1.07154e-07,-4.12128e-11,-36859.9,27.1961], Tmin=(100,'K'), Tmax=(778.721,'K')), NASAPolynomial(coeffs=[4.63588,0.0390318,-2.05667e-05,4.09767e-09,-2.90266e-13,-36959.3,13.9751], Tmin=(778.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-307.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(CdCFO) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C](F)OOC=C(1691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u1 p0 c0 {1,S} {3,S} {6,S}
5  C u0 p0 c0 {2,S} {7,D} {8,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (103.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,395,473,707,1436,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,207.132],'cm^-1')),
        HinderedRotor(inertia=(0.401638,'amu*angstrom^2'), symmetry=1, barrier=(12.228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949057,'amu*angstrom^2'), symmetry=1, barrier=(28.8942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401636,'amu*angstrom^2'), symmetry=1, barrier=(12.228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04156,'amu*angstrom^2'), symmetry=1, barrier=(31.7105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02716,0.0687382,-7.87916e-05,4.74893e-08,-1.161e-11,12496.6,26.9984], Tmin=(100,'K'), Tmax=(985.612,'K')), NASAPolynomial(coeffs=[11.822,0.0249286,-1.2118e-05,2.3915e-09,-1.7099e-13,10368.7,-24.9242], Tmin=(985.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CsCsF1sO2s) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C(F)OO[C]=C(1692)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u1 p0 c0 {4,S} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (148.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.513121,'amu*angstrom^2'), symmetry=1, barrier=(11.7977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50087,'amu*angstrom^2'), symmetry=1, barrier=(34.508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50111,'amu*angstrom^2'), symmetry=1, barrier=(34.5134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114782,'amu*angstrom^2'), symmetry=1, barrier=(2.63906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40267,0.0615739,-6.37862e-05,3.59501e-08,-8.496e-12,17917.3,27.9786], Tmin=(100,'K'), Tmax=(998.729,'K')), NASAPolynomial(coeffs=[9.65439,0.0285252,-1.41504e-05,2.81757e-09,-2.02376e-13,16269,-11.8207], Tmin=(998.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(C=CJO)"""),
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
    label = '[CH2][CH]F(839)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (116.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0021411,'amu*angstrom^2'), symmetry=1, barrier=(6.24533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46308e-06,3.53001e-10,6.36801e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95459,0.00674552,-2.7461e-06,4.76057e-10,-2.99989e-14,11484.3,-14.7486], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(Cs_P)"""),
)

species(
    label = 'CH2CF(69)',
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
    label = '[CH2]C(F)O[O](846)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-48.0821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.486055,'amu*angstrom^2'), symmetry=1, barrier=(11.1754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485983,'amu*angstrom^2'), symmetry=1, barrier=(11.1737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33352,0.040021,-5.12824e-05,3.75377e-08,-1.13901e-11,-5725.96,17.2582], Tmin=(100,'K'), Tmax=(795.131,'K')), NASAPolynomial(coeffs=[6.83717,0.0173642,-8.53954e-06,1.69956e-09,-1.21736e-13,-6442.14,-3.43674], Tmin=(795.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][C](F)OOC(=C)F(1461)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u1 p0 c0 {1,S} {3,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-82.1702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,395,473,707,1436,326,540,652,719,1357,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.162383,'amu*angstrom^2'), symmetry=1, barrier=(3.73351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514274,'amu*angstrom^2'), symmetry=1, barrier=(11.8242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161665,'amu*angstrom^2'), symmetry=1, barrier=(3.717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18868,'amu*angstrom^2'), symmetry=1, barrier=(50.322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51958,0.0863625,-0.000145194,1.34543e-07,-4.90431e-11,-9766.98,30.2135], Tmin=(100,'K'), Tmax=(804.691,'K')), NASAPolynomial(coeffs=[7.30524,0.0356481,-1.89995e-05,3.76598e-09,-2.64642e-13,-10309.2,2.36779], Tmin=(804.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.1702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CsCsF1sO2s) + radical(CJCOOH)"""),
)

species(
    label = '[CH]=C(F)OOC([CH2])F(4393)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
6  C u1 p0 c0 {5,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u1 p0 c0 {7,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-14.2732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000,293,496,537,1218,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.103573,'amu*angstrom^2'), symmetry=1, barrier=(2.38135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428387,'amu*angstrom^2'), symmetry=1, barrier=(9.84945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85208,'amu*angstrom^2'), symmetry=1, barrier=(42.583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85184,'amu*angstrom^2'), symmetry=1, barrier=(42.5775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659118,0.0807498,-0.000114307,8.98253e-08,-2.90488e-11,-1603.05,29.0414], Tmin=(100,'K'), Tmax=(749.486,'K')), NASAPolynomial(coeffs=[9.74741,0.0322474,-1.72392e-05,3.48683e-09,-2.50676e-13,-2965.41,-12.184], Tmin=(749.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.2732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CJCOOH) + radical(Cdj(Cd-F1sO2s)(H))"""),
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
    label = 'C#COOC([CH2])F(1513)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u1 p0 c0 {4,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (139.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.49473,'amu*angstrom^2'), symmetry=1, barrier=(34.3668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49562,'amu*angstrom^2'), symmetry=1, barrier=(34.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49704,'amu*angstrom^2'), symmetry=1, barrier=(34.4198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49677,'amu*angstrom^2'), symmetry=1, barrier=(34.4137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3238.92,'J/mol'), sigma=(5.49609,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=505.91 K, Pc=44.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1087,0.0670647,-7.60749e-05,4.4181e-08,-1.03782e-11,16939.1,23.5113], Tmin=(100,'K'), Tmax=(1023.38,'K')), NASAPolynomial(coeffs=[12.3003,0.0233208,-1.19578e-05,2.41264e-09,-1.74592e-13,14648.5,-30.7403], Tmin=(1023.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(CsCFHO) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCOOH)"""),
)

species(
    label = 'C=C(F)OO[C](C)F(4394)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u1 p0 c0 {1,S} {3,S} {5,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-296.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.26696,'amu*angstrom^2'), symmetry=1, barrier=(6.13795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567272,'amu*angstrom^2'), symmetry=1, barrier=(13.0427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266829,'amu*angstrom^2'), symmetry=1, barrier=(6.13494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1691,'amu*angstrom^2'), symmetry=1, barrier=(49.8718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622606,0.0832072,-0.000130607,1.20401e-07,-4.47082e-11,-35503.6,27.8485], Tmin=(100,'K'), Tmax=(782.517,'K')), NASAPolynomial(coeffs=[6.37099,0.0398866,-2.08512e-05,4.13442e-09,-2.91881e-13,-35976.5,4.25196], Tmin=(782.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH]=C(F)OOC(C)F(4395)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u1 p0 c0 {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-228.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.118526,'amu*angstrom^2'), symmetry=1, barrier=(2.72515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523792,'amu*angstrom^2'), symmetry=1, barrier=(12.043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74956,'amu*angstrom^2'), symmetry=1, barrier=(40.2257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75309,'amu*angstrom^2'), symmetry=1, barrier=(40.3071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870567,0.0761668,-9.39226e-05,6.67472e-08,-2.01183e-11,-27344.3,26.2965], Tmin=(100,'K'), Tmax=(792.104,'K')), NASAPolynomial(coeffs=[8.68185,0.0367235,-1.92339e-05,3.8901e-09,-2.80869e-13,-28581.8,-9.56843], Tmin=(792.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'C=C(F)OO[CH]CF(4396)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {5,S} {11,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-267.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.169874,'amu*angstrom^2'), symmetry=1, barrier=(3.90573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10647,'amu*angstrom^2'), symmetry=1, barrier=(48.4319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10784,'amu*angstrom^2'), symmetry=1, barrier=(48.4634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171233,'amu*angstrom^2'), symmetry=1, barrier=(3.93697,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651005,0.0821154,-0.000125945,1.14865e-07,-4.25227e-11,-32045.1,28.0079], Tmin=(100,'K'), Tmax=(776.83,'K')), NASAPolynomial(coeffs=[6.43452,0.0399819,-2.0735e-05,4.10404e-09,-2.89737e-13,-32570.9,3.96541], Tmin=(776.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFHH) + group(CdCFO) + group(Cds-CdsHH) + radical(CCsJOOC)"""),
)

species(
    label = 'C=[C]OOC(F)CF(4397)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
7  C u0 p0 c0 {8,D} {12,S} {13,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-239.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,2950,3100,1380,975,1025,1650,1685,370,237.385,239.319],'cm^-1')),
        HinderedRotor(inertia=(0.908811,'amu*angstrom^2'), symmetry=1, barrier=(36.889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523955,'amu*angstrom^2'), symmetry=1, barrier=(21.2189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344278,'amu*angstrom^2'), symmetry=1, barrier=(14.1621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902647,'amu*angstrom^2'), symmetry=1, barrier=(36.885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12929,0.0680836,-6.7755e-05,3.61187e-08,-8.08279e-12,-28754,28.2722], Tmin=(100,'K'), Tmax=(1047.98,'K')), NASAPolynomial(coeffs=[10.6003,0.0319338,-1.60128e-05,3.20311e-09,-2.30611e-13,-30739.1,-17.8639], Tmin=(1047.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
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
    E0 = (-131.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-96.3025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (219.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (250.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-168.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-168.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (83.1725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-140.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (14.2895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (256.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (302.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-154.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (100.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (129.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (210.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (278.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-184.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (227.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-30.5555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-98.9575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-37.6382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (21.6224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)OOC(=C)F(3748)'],
    products = ['CH2CFO(78)', 'FC1CO1(1534)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(63.953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for R2OO_S;C_pri_rad_intra;OOR
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)OOC(=C)F(3748)'],
    products = ['[CH2]C(F)OCC(=O)F(4389)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(99.4845,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(18)', 'C=C(F)OO[CH]F(2054)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(5)', '[CH]C(F)OOC(=C)F(4390)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)OOC(=C)F(3748)'],
    products = ['F[C]1CCC(F)OO1(4391)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.48374e+06,'s^-1'), n=1.00311, Ea=(26.8393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SSS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)OOC(=C)F(3748)'],
    products = ['[CH2]C1(F)CC(F)OO1(4392)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.79588e+07,'s^-1'), n=0.7, Ea=(27.3006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SSS_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'C=COOC(=C)F(1690)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.4037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C(F)O[O](1327)', 'CH2CHF(55)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.000486681,'m^3/(mol*s)'), n=2.49948, Ea=(27.9764,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Sp-2CNS=1R!H_Ext-1R!H-R_N-5R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Sp-2CNS=1R!H_Ext-1R!H-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C=C(F)OOC(=C)F(1464)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.01632e+22,'m^3/(mol*s)'), n=-4.30708, Ea=(28.8132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2][C](F)OOC=C(1691)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2]C(F)OO[C]=C(1692)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2CFO(78)', '[CH2]C([O])F(527)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(75.0814,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(F)O[O](1327)', '[CH2][CH]F(839)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2CF(69)', '[CH2]C(F)O[O](846)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH2][C](F)OOC(=C)F(1461)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.90705e+22,'m^3/(mol*s)'), n=-6.09098, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Sp-4C-2C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Sp-4C-2C_Ext-2C-R
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH]=C(F)OOC([CH2])F(4393)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.84097e+21,'m^3/(mol*s)'), n=-5.92191, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3347815616900648, var=90.07385025150133, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_Sp-3C-2C_Ext-3C-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(F)OOC(=C)F(3748)'],
    products = ['CH2CFO(78)', 'CC(=O)F(499)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(11.3907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'C#COOC([CH2])F(1513)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.56721e+32,'m^3/(mol*s)'), n=-7.40875, Ea=(287.854,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_3COCdCddCtO2d->Ct',), comment="""Estimated from node HF_3COCdCddCtO2d->Ct"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(F)OO[C](C)F(4394)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.465e+11,'s^-1'), n=0, Ea=(184.64,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 440 used for R2H_S;C_rad_out_noH;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(F)OOC(C)F(4395)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21892e+08,'s^-1'), n=0.873546, Ea=(48.341,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;Cd_rad_out_singleH;XH_out] + [R6H_RSSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_DSSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.605551275463989
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(F)OO[CH]CF(4396)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(148.797,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OOC(F)CF(4397)'],
    products = ['[CH2]C(F)OOC(=C)F(3748)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(180.585,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1422',
    isomers = [
        '[CH2]C(F)OOC(=C)F(3748)',
    ],
    reactants = [
        ('CH2CFO(78)', 'FC1CO1(1534)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1422',
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

