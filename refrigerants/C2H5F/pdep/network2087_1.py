species(
    label = 'CC(OF)C(C)(O)F(5983)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {13,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {4,S}
"""),
    E0 = (-598.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,333,384,448,608,1254,1480,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.697504,'amu*angstrom^2'), symmetry=1, barrier=(16.037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697515,'amu*angstrom^2'), symmetry=1, barrier=(16.0372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847071,'amu*angstrom^2'), symmetry=1, barrier=(16.0369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489851,'amu*angstrom^2'), symmetry=1, barrier=(11.2626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697504,'amu*angstrom^2'), symmetry=1, barrier=(16.037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.102,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230767,0.0857618,-8.99017e-05,5.00205e-08,-1.13667e-11,-71859.7,26.4272], Tmin=(100,'K'), Tmax=(1052.99,'K')), NASAPolynomial(coeffs=[13.9383,0.0336912,-1.57268e-05,3.05925e-09,-2.17259e-13,-74746.5,-40.4114], Tmin=(1052.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cs-CsHHH)"""),
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
    label = 'CH3CHO(36)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.65,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'CC([C]F)OF(2299)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p1 c0 {1,S} {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (19.1662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,617,898,1187,229.256],'cm^-1')),
        HinderedRotor(inertia=(0.247486,'amu*angstrom^2'), symmetry=1, barrier=(8.75617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234453,'amu*angstrom^2'), symmetry=1, barrier=(8.73007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68503,'amu*angstrom^2'), symmetry=1, barrier=(25.3952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80829,0.04975,-5.24911e-05,2.96425e-08,-6.84026e-12,2382.83,20.0099], Tmin=(100,'K'), Tmax=(1038.44,'K')), NASAPolynomial(coeffs=[9.6195,0.0196617,-9.02913e-06,1.74031e-09,-1.22924e-13,760.541,-17.9693], Tmin=(1038.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.1662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'CC(O)OF(1543)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {1,S} {4,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-336.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.346091,'amu*angstrom^2'), symmetry=1, barrier=(7.95732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346792,'amu*angstrom^2'), symmetry=1, barrier=(7.97344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896703,'amu*angstrom^2'), symmetry=1, barrier=(20.617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0583,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78744,0.0177066,0.000134802,-4.71341e-07,4.67563e-10,-40514.8,9.53186], Tmin=(10,'K'), Tmax=(354.107,'K')), NASAPolynomial(coeffs=[4.66427,0.0310735,-2.03978e-05,6.43887e-09,-7.76868e-13,-40722.8,4.15181], Tmin=(354.107,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-336.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""CC(O)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CC(O)(F)COF(5726)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {13,S}
4  O u0 p2 c0 {2,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-562.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,185.756,187.807],'cm^-1')),
        HinderedRotor(inertia=(0.515471,'amu*angstrom^2'), symmetry=1, barrier=(13.0301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504348,'amu*angstrom^2'), symmetry=1, barrier=(13.0207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523445,'amu*angstrom^2'), symmetry=1, barrier=(13.0389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863522,'amu*angstrom^2'), symmetry=1, barrier=(21.8143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.075,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3565.28,'J/mol'), sigma=(6.03172,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=556.89 K, Pc=36.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882636,0.0730326,-8.61579e-05,5.50449e-08,-1.43895e-11,-67504,22.1932], Tmin=(100,'K'), Tmax=(921.55,'K')), NASAPolynomial(coeffs=[11.1217,0.0285901,-1.38198e-05,2.71461e-09,-1.93319e-13,-69391.2,-26.3681], Tmin=(921.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-562.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cs-CsHHH)"""),
)

species(
    label = 'CC(OF)C(O)F(2878)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-544.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,190.539,1199.82],'cm^-1')),
        HinderedRotor(inertia=(0.384609,'amu*angstrom^2'), symmetry=1, barrier=(9.89508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690233,'amu*angstrom^2'), symmetry=1, barrier=(17.8487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691007,'amu*angstrom^2'), symmetry=1, barrier=(17.8504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691507,'amu*angstrom^2'), symmetry=1, barrier=(17.8484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.075,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3531.24,'J/mol'), sigma=(5.86554,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.57 K, Pc=39.71 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01043,0.0666709,-6.73404e-05,3.55312e-08,-7.61305e-12,-65323.9,23.4977], Tmin=(100,'K'), Tmax=(1116.31,'K')), NASAPolynomial(coeffs=[12.5853,0.0251955,-1.16098e-05,2.2487e-09,-1.59403e-13,-67908.2,-33.6181], Tmin=(1116.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-544.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(Cs-CsHHH)"""),
)

species(
    label = 'CC(F)C(C)(O)OF(5993)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {5,S} {16,S}
5  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {4,S}
"""),
    E0 = (-606.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.102,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5399,0.0976501,-0.000100682,4.88678e-08,-8.86559e-12,-72751.1,31.7618], Tmin=(100,'K'), Tmax=(1528.76,'K')), NASAPolynomial(coeffs=[27.654,0.00868258,-1.04288e-06,3.35612e-11,7.05516e-16,-80207,-116.665], Tmin=(1528.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(Cs-CsHHH) + group(Cs-CsHHH)"""),
)

species(
    label = 'FOF(1637)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {1,S} {2,S}
"""),
    E0 = (16.0257,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(53.9917,'amu')),
        NonlinearRotor(inertia=([8.34135,45.8552,54.1965],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([486.983,915.713,1034.58],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (53.9962,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03449,-0.00336126,4.20188e-05,-7.82672e-08,4.57369e-11,1928.28,6.37844], Tmin=(10,'K'), Tmax=(574.663,'K')), NASAPolynomial(coeffs=[3.91346,0.00598694,-4.58407e-06,1.5533e-09,-1.93094e-13,1801.75,5.67331], Tmin=(574.663,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(16.0257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""FOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC=C(C)O(5061)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {13,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-244.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.616792,'amu*angstrom^2'), symmetry=1, barrier=(14.1813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614455,'amu*angstrom^2'), symmetry=1, barrier=(14.1275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616224,'amu*angstrom^2'), symmetry=1, barrier=(14.1682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63688,0.0428964,-1.05625e-05,-1.89796e-08,1.10864e-11,-29308.6,17.5105], Tmin=(100,'K'), Tmax=(965.183,'K')), NASAPolynomial(coeffs=[12.6347,0.0187514,-6.34798e-06,1.11669e-09,-7.82684e-14,-32429.9,-40.3296], Tmin=(965.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH)"""),
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
    label = 'C=CC(C)(O)F(5191)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {13,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-433.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.02147,'amu*angstrom^2'), symmetry=1, barrier=(23.4855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02294,'amu*angstrom^2'), symmetry=1, barrier=(23.5195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01742,'amu*angstrom^2'), symmetry=1, barrier=(23.3925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0962,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4593,0.0544848,-4.35367e-05,1.80149e-08,-3.07084e-12,-52102.7,19.2831], Tmin=(100,'K'), Tmax=(1362.09,'K')), NASAPolynomial(coeffs=[11.5345,0.0248972,-1.09531e-05,2.06694e-09,-1.43716e-13,-54847.4,-32.4374], Tmin=(1362.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-433.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = 'CC(F)=C(C)OF(5074)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {1,S} {5,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (-268.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,323,467,575,827,1418,180],'cm^-1')),
        HinderedRotor(inertia=(0.174871,'amu*angstrom^2'), symmetry=1, barrier=(4.02062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174659,'amu*angstrom^2'), symmetry=1, barrier=(4.01575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174768,'amu*angstrom^2'), symmetry=1, barrier=(4.01825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20197,0.0690122,-0.000101023,9.31769e-08,-3.49969e-11,-32182.7,22.2404], Tmin=(100,'K'), Tmax=(791.851,'K')), NASAPolynomial(coeffs=[4.63102,0.0387753,-1.92806e-05,3.75933e-09,-2.63344e-13,-32320.8,9.05421], Tmin=(791.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'C=C(F)C(C)OF(5075)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-263.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,618.629],'cm^-1')),
        HinderedRotor(inertia=(0.453039,'amu*angstrom^2'), symmetry=1, barrier=(13.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0483961,'amu*angstrom^2'), symmetry=1, barrier=(13.1418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.571596,'amu*angstrom^2'), symmetry=1, barrier=(13.1421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948546,0.0616976,-5.39545e-05,2.37271e-08,-4.16423e-12,-31555.7,24.5707], Tmin=(100,'K'), Tmax=(1362.97,'K')), NASAPolynomial(coeffs=[14.9243,0.0206816,-8.81436e-06,1.64759e-09,-1.143e-13,-35365.4,-47.1821], Tmin=(1362.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(F)=COF(2323)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-250.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.0691364,'amu*angstrom^2'), symmetry=1, barrier=(1.58958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717506,'amu*angstrom^2'), symmetry=1, barrier=(16.4969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96351,0.00554679,0.000501419,-3.99572e-06,1.14263e-08,-30102,6.51516], Tmin=(10,'K'), Tmax=(87.5713,'K')), NASAPolynomial(coeffs=[3.29096,0.036261,-2.45805e-05,7.85978e-09,-9.53058e-13,-30090.3,8.12209], Tmin=(87.5713,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-250.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CC(F)DCOF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C[C](O)C(C)OF(5065)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
7  C u1 p0 c0 {3,S} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-178.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,298.199,298.202],'cm^-1')),
        HinderedRotor(inertia=(0.00189574,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169702,'amu*angstrom^2'), symmetry=1, barrier=(10.7087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169702,'amu*angstrom^2'), symmetry=1, barrier=(10.7087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169705,'amu*angstrom^2'), symmetry=1, barrier=(10.7087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169701,'amu*angstrom^2'), symmetry=1, barrier=(10.7087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.104,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518223,0.0810265,-9.91736e-05,6.84266e-08,-1.93435e-11,-21300.2,26.4258], Tmin=(100,'K'), Tmax=(856.599,'K')), NASAPolynomial(coeffs=[10.8275,0.0328865,-1.48755e-05,2.82009e-09,-1.96263e-13,-23066.3,-21.7146], Tmin=(856.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH)"""),
)

species(
    label = 'CC([O])C(C)(O)F(5994)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {15,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {2,S}
"""),
    E0 = (-502.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,333,384,448,608,1254,1480,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,308.506,308.848],'cm^-1')),
        HinderedRotor(inertia=(0.209371,'amu*angstrom^2'), symmetry=1, barrier=(14.1224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209987,'amu*angstrom^2'), symmetry=1, barrier=(14.1222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209197,'amu*angstrom^2'), symmetry=1, barrier=(14.1224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208856,'amu*angstrom^2'), symmetry=1, barrier=(14.1223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.104,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696242,0.0715793,-6.54194e-05,3.12949e-08,-6.10382e-12,-60344.9,24.384], Tmin=(100,'K'), Tmax=(1216.64,'K')), NASAPolynomial(coeffs=[13.6996,0.0288275,-1.27106e-05,2.41282e-09,-1.69016e-13,-63509,-40.8997], Tmin=(1216.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-502.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ)"""),
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
    label = 'C[CH]C(C)(O)F(5995)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {14,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {5,S} {13,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-360.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,788.596],'cm^-1')),
        HinderedRotor(inertia=(0.640322,'amu*angstrom^2'), symmetry=1, barrier=(14.7223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164198,'amu*angstrom^2'), symmetry=1, barrier=(3.77524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163538,'amu*angstrom^2'), symmetry=1, barrier=(3.76005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.640304,'amu*angstrom^2'), symmetry=1, barrier=(14.7219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05557,0.057342,-4.33545e-05,1.65885e-08,-2.55376e-12,-43272.6,23.435], Tmin=(100,'K'), Tmax=(1531.44,'K')), NASAPolynomial(coeffs=[14.5052,0.0222122,-8.94552e-06,1.60949e-09,-1.08479e-13,-47392,-47.1841], Tmin=(1531.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO)"""),
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
    label = 'C[C](F)C(C)OF(2305)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
7  C u1 p0 c0 {1,S} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-204.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,212,367,445,1450,283.759,283.767],'cm^-1')),
        HinderedRotor(inertia=(0.00209774,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226359,'amu*angstrom^2'), symmetry=1, barrier=(12.937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1543,'amu*angstrom^2'), symmetry=1, barrier=(8.79463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154365,'amu*angstrom^2'), symmetry=1, barrier=(8.79896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900817,0.0738087,-9.33011e-05,7.12431e-08,-2.29649e-11,-24438.4,25.1765], Tmin=(100,'K'), Tmax=(746.744,'K')), NASAPolynomial(coeffs=[7.96893,0.0359485,-1.72523e-05,3.35095e-09,-2.35969e-13,-25494,-6.85905], Tmin=(746.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CsCsCsF1s)"""),
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
    label = 'CC(OF)C(C)([O])F(5996)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {13,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-364.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,286,334,437,422,648,895,1187,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,224.893,225.141],'cm^-1')),
        HinderedRotor(inertia=(0.00333794,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302384,'amu*angstrom^2'), symmetry=1, barrier=(10.866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303269,'amu*angstrom^2'), symmetry=1, barrier=(10.8647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301687,'amu*angstrom^2'), symmetry=1, barrier=(10.864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.466282,0.0830519,-9.64964e-05,6.22132e-08,-1.65701e-11,-43732.6,26.732], Tmin=(100,'K'), Tmax=(901.941,'K')), NASAPolynomial(coeffs=[11.2771,0.0351058,-1.67557e-05,3.27134e-09,-2.32102e-13,-45682.7,-24.3079], Tmin=(901.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sCsCs))"""),
)

species(
    label = 'C[CH]OF(540)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (46.0665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,252.797],'cm^-1')),
        HinderedRotor(inertia=(0.00264785,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296597,'amu*angstrom^2'), symmetry=1, barrier=(13.4557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49569,0.0286678,-1.90921e-05,3.84508e-09,6.27455e-13,5598.43,14.8818], Tmin=(100,'K'), Tmax=(1090.57,'K')), NASAPolynomial(coeffs=[9.35696,0.01068,-4.22389e-06,7.913e-10,-5.5993e-14,3675.04,-20.7717], Tmin=(1090.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.0665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
)

species(
    label = 'C[C](O)F(1910)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {8,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-300.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436],'cm^-1')),
        HinderedRotor(inertia=(0.496378,'amu*angstrom^2'), symmetry=1, barrier=(11.4127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.495871,'amu*angstrom^2'), symmetry=1, barrier=(11.401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89273,0.00700432,9.57436e-05,-2.34846e-07,1.68061e-10,-36154.7,9.17454], Tmin=(10,'K'), Tmax=(488.816,'K')), NASAPolynomial(coeffs=[4.86387,0.0207301,-1.28815e-05,4.00364e-09,-4.85912e-13,-36508.5,2.53597], Tmin=(488.816,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-300.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""C[C](O)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CC(O)(F)[CH]OF(5818)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {12,S}
4  O u0 p2 c0 {2,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-378.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.763664,'amu*angstrom^2'), symmetry=1, barrier=(17.5581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472801,'amu*angstrom^2'), symmetry=1, barrier=(10.8706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0908322,'amu*angstrom^2'), symmetry=1, barrier=(17.5604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763743,'amu*angstrom^2'), symmetry=1, barrier=(17.56,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632726,0.077232,-9.78978e-05,6.30222e-08,-1.60904e-11,-45397.1,23.2888], Tmin=(100,'K'), Tmax=(955.518,'K')), NASAPolynomial(coeffs=[13.9328,0.0215665,-1.05305e-05,2.07852e-09,-1.48512e-13,-47939.3,-40.274], Tmin=(955.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJO)"""),
)

species(
    label = 'C[C](OF)C(C)(O)F(5862)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {15,S}
4  O u0 p2 c0 {2,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {8,S} {12,S} {13,S} {14,S}
8  C u1 p0 c0 {4,S} {5,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-418.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,316,385,515,654,689,1295,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,180.006,180.008],'cm^-1')),
        HinderedRotor(inertia=(0.0676435,'amu*angstrom^2'), symmetry=1, barrier=(1.56167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376121,'amu*angstrom^2'), symmetry=1, barrier=(8.64787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380032,'amu*angstrom^2'), symmetry=1, barrier=(8.7377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.328015,'amu*angstrom^2'), symmetry=1, barrier=(8.71922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.849845,'amu*angstrom^2'), symmetry=1, barrier=(20.4066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162974,0.0923006,-0.000124373,8.98027e-08,-2.47382e-11,-50200.7,27.1097], Tmin=(100,'K'), Tmax=(655.63,'K')), NASAPolynomial(coeffs=[10.9217,0.0361679,-1.76974e-05,3.44687e-09,-2.42408e-13,-51815.8,-21.8111], Tmin=(655.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-418.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJO)"""),
)

species(
    label = 'CC(OF)[C](O)F(3199)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {4,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-349.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,259.957,260.869],'cm^-1')),
        HinderedRotor(inertia=(0.173373,'amu*angstrom^2'), symmetry=1, barrier=(8.39454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17371,'amu*angstrom^2'), symmetry=1, barrier=(8.38754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383818,'amu*angstrom^2'), symmetry=1, barrier=(18.5627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173455,'amu*angstrom^2'), symmetry=1, barrier=(8.39302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961483,0.0736847,-9.9508e-05,6.73616e-08,-1.53365e-11,-41928.4,24.5053], Tmin=(100,'K'), Tmax=(630.059,'K')), NASAPolynomial(coeffs=[10.0093,0.0268124,-1.3079e-05,2.53398e-09,-1.77306e-13,-43278.3,-16.6303], Tmin=(630.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-349.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(Cs-CsHHH) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH2]C(OF)C(C)(O)F(5997)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {5,S} {15,S}
4  O u0 p2 c0 {2,S} {6,S}
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
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-386.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,333,384,448,608,1254,1480,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,277.865,282.067],'cm^-1')),
        HinderedRotor(inertia=(0.722371,'amu*angstrom^2'), symmetry=1, barrier=(18.3808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346296,'amu*angstrom^2'), symmetry=1, barrier=(18.4105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138167,'amu*angstrom^2'), symmetry=1, barrier=(7.6913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322392,'amu*angstrom^2'), symmetry=1, barrier=(18.3764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555971,'amu*angstrom^2'), symmetry=1, barrier=(13.0063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.033609,0.0914152,-0.00011101,7.06591e-08,-1.80713e-11,-46405.5,27.7733], Tmin=(100,'K'), Tmax=(949.013,'K')), NASAPolynomial(coeffs=[14.5293,0.0303173,-1.44393e-05,2.81983e-09,-2.00282e-13,-49156.8,-41.4015], Tmin=(949.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-386.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)(F)C(C)OF(5998)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (-387.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,256.49,257.686],'cm^-1')),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(8.53471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185207,'amu*angstrom^2'), symmetry=1, barrier=(8.54379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182711,'amu*angstrom^2'), symmetry=1, barrier=(8.53683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352766,'amu*angstrom^2'), symmetry=1, barrier=(16.4203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662762,'amu*angstrom^2'), symmetry=1, barrier=(31.085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0919393,0.0967127,-0.000133105,9.9689e-08,-3.02068e-11,-46492.4,27.8404], Tmin=(100,'K'), Tmax=(804.066,'K')), NASAPolynomial(coeffs=[12.4343,0.0343944,-1.68418e-05,3.28708e-09,-2.31713e-13,-48506.6,-29.8591], Tmin=(804.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'CC(=O)C(C)OF(3599)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-293.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,797.78,801.557],'cm^-1')),
        HinderedRotor(inertia=(0.0296771,'amu*angstrom^2'), symmetry=1, barrier=(13.0759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570481,'amu*angstrom^2'), symmetry=1, barrier=(13.1165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569998,'amu*angstrom^2'), symmetry=1, barrier=(13.1054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146832,'amu*angstrom^2'), symmetry=1, barrier=(3.37595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.096,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.78,'J/mol'), sigma=(5.79777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.16 K, Pc=40.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00172,0.0570395,-3.69733e-05,8.05708e-09,5.6478e-13,-35192.7,26.7375], Tmin=(100,'K'), Tmax=(1156.31,'K')), NASAPolynomial(coeffs=[14.31,0.0235866,-9.90177e-06,1.86088e-09,-1.30443e-13,-39111.6,-43.0377], Tmin=(1156.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs)"""),
)

species(
    label = 'CC(=O)C(C)(O)F(5328)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {14,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-661.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.749109,'amu*angstrom^2'), symmetry=1, barrier=(17.2235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.747847,'amu*angstrom^2'), symmetry=1, barrier=(17.1945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.748856,'amu*angstrom^2'), symmetry=1, barrier=(17.2177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.747769,'amu*angstrom^2'), symmetry=1, barrier=(17.1927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.096,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16759,0.0648001,-6.05312e-05,3.07769e-08,-6.54894e-12,-79429.1,22.0407], Tmin=(100,'K'), Tmax=(1103.3,'K')), NASAPolynomial(coeffs=[10.5854,0.0306555,-1.4109e-05,2.72596e-09,-1.92718e-13,-81507.3,-24.3205], Tmin=(1103.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-661.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs)"""),
)

species(
    label = 'CC(O)=C(C)OF(3618)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,S} {14,S}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-301.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,295.977],'cm^-1')),
        HinderedRotor(inertia=(0.184456,'amu*angstrom^2'), symmetry=1, barrier=(11.4634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184997,'amu*angstrom^2'), symmetry=1, barrier=(11.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183594,'amu*angstrom^2'), symmetry=1, barrier=(11.4615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184628,'amu*angstrom^2'), symmetry=1, barrier=(11.4627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.096,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568427,0.0726499,-7.58874e-05,4.10153e-08,-8.79376e-12,-36101.5,23.8714], Tmin=(100,'K'), Tmax=(1134.58,'K')), NASAPolynomial(coeffs=[14.9571,0.0219222,-8.8216e-06,1.60817e-09,-1.10571e-13,-39366.6,-47.3624], Tmin=(1134.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-301.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs)"""),
)

species(
    label = 'C=C(O)C(C)OF(3594)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {6,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-269.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,267.013,267.032],'cm^-1')),
        HinderedRotor(inertia=(0.322461,'amu*angstrom^2'), symmetry=1, barrier=(16.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322563,'amu*angstrom^2'), symmetry=1, barrier=(16.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322486,'amu*angstrom^2'), symmetry=1, barrier=(16.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322529,'amu*angstrom^2'), symmetry=1, barrier=(16.3182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.096,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3674.47,'J/mol'), sigma=(6.07294,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.94 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434212,0.0664334,-4.29923e-05,-2.20015e-09,8.37269e-12,-32259.6,26.372], Tmin=(100,'K'), Tmax=(966.25,'K')), NASAPolynomial(coeffs=[19.7723,0.0139995,-4.47215e-06,8.0641e-10,-5.94416e-14,-37286.1,-72.9313], Tmin=(966.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)F(498)',
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
    E0 = (-376.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (167.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (126.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (50.8058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (68.9446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-105.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (197.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-123.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-94.4436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-118.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (61.0874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (88.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-235.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-64.1577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (18.1944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (41.0473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-60.7114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-48.3752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-12.7989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-19.4047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (18.6854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (17.9473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-222.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-280.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-133.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-143.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-349.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(OF)C(C)(O)F(5983)'],
    products = ['HF(38)', 'CC(=O)F(499)', 'CH3CHO(36)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(28.494,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH3OH(26)', 'CC([C]F)OF(2299)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(166.188,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C[C]F(124)', 'CC(O)OF(1543)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(189.05,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(S)(25)', 'CC(O)(F)COF(5726)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(980214,'m^3/(mol*s)'), n=0.201831, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.043461663341785264, var=0.18543765059589312, Tref=1000.0, N=2, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(S)(25)', 'CC(OF)C(O)F(2878)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(490107,'m^3/(mol*s)'), n=0.201831, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.043461663341785264, var=0.18543765059589312, Tref=1000.0, N=2, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC(OF)C(C)(O)F(5983)'],
    products = ['CC(F)C(C)(O)OF(5993)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(1637)', 'CC=C(C)O(5061)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/Nd2_Cd/H/Nd;H_OH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OF(169)', 'C=CC(C)(O)F(5191)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(65.2,'cm^3/(mol*s)'), n=2.92, Ea=(212.129,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 9 used for Cd/H2_Cd/H/Nd;H_OH
Exact match found for rate rule [Cd/H2_Cd/H/Nd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', 'CC(F)=C(C)OF(5074)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H2O(3)', 'C=C(F)C(C)OF(5075)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2500,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH3OH(26)', 'CC(F)=COF(2323)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.79e-05,'cm^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_Cd;CH3OH] for rate rule [Cd/monosub_Cd/disub;CH3OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C[C](O)C(C)OF(5065)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'CC([O])C(C)(O)F(5994)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(160)', 'C[CH]C(C)(O)F(5995)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OH(7)', 'C[C](F)C(C)OF(2305)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'CC(OF)C(C)([O])F(5996)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]OF(540)', 'C[C](O)F(1910)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -0.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(19)', 'CC(O)(F)[CH]OF(5818)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'C[C](OF)C(C)(O)F(5862)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH3(19)', 'CC(OF)[C](O)F(3199)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.29446e+18,'m^3/(mol*s)'), n=-4.19701, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R
Ea raised from -10.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[CH2]C(OF)C(C)(O)F(5997)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[CH2]C(O)(F)C(C)OF(5998)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'CC(=O)C(C)OF(3599)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(158.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'CC(=O)C(C)(O)F(5328)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(467.884,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'CC(O)=C(C)OF(3618)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(254.749,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'C=C(O)C(C)OF(3594)'],
    products = ['CC(OF)C(C)(O)F(5983)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(213.36,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC(OF)C(C)(O)F(5983)'],
    products = ['HF(38)', 'C=C(O)F(498)', 'CH3CHO(36)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(55.2542,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #2087',
    isomers = [
        'CC(OF)C(C)(O)F(5983)',
    ],
    reactants = [
        ('HF(38)', 'CC(=O)F(499)', 'CH3CHO(36)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2087',
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

