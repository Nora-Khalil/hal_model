species(
    label = '[CH2]C=CCC[C]=O(1862)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u1 p0 c0 {5,S} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (152.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,261.019,4000],'cm^-1')),
        HinderedRotor(inertia=(0.790854,'amu*angstrom^2'), symmetry=1, barrier=(38.2159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15799,'amu*angstrom^2'), symmetry=1, barrier=(7.63498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30808,'amu*angstrom^2'), symmetry=1, barrier=(14.8898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398608,'amu*angstrom^2'), symmetry=1, barrier=(19.2647,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2004,0.0606472,-4.66902e-05,1.938e-08,-3.38164e-12,18422.7,27.3071], Tmin=(100,'K'), Tmax=(1317.42,'K')), NASAPolynomial(coeffs=[10.8801,0.0312573,-1.32271e-05,2.44636e-09,-1.68236e-13,15872.3,-22.0603], Tmin=(1317.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'CH2CO(27)',
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
    label = 'C=CC=C(61)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,D} {5,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {7,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80601,0.0102582,6.17268e-05,-9.01653e-08,3.59121e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.045,'K')), NASAPolynomial(coeffs=[12.4693,0.0100555,-2.41212e-06,4.5709e-10,-3.93171e-14,8010.8,-43.6372], Tmin=(946.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][C]=O(167)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.3074e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'CH2(S)(24)',
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
    label = '[CH2]C=CC[C]=O(1300)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u1 p0 c0 {4,S} {11,S} {12,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (182.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,315.93],'cm^-1')),
        HinderedRotor(inertia=(0.0017048,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312051,'amu*angstrom^2'), symmetry=1, barrier=(21.9843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313377,'amu*angstrom^2'), symmetry=1, barrier=(22.0269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84068,0.0373102,-1.18218e-06,-2.21776e-08,9.72893e-12,22011.4,22.3283], Tmin=(100,'K'), Tmax=(1088.78,'K')), NASAPolynomial(coeffs=[12.2928,0.0212541,-9.84398e-06,1.97424e-09,-1.44531e-13,18411.1,-35.0678], Tmin=(1088.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=CCC(=C)[O](1858)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {3,D} {7,S} {11,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (127.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924688,0.0562639,-1.95361e-05,-1.80404e-08,1.18946e-11,15451.9,27.5064], Tmin=(100,'K'), Tmax=(976.244,'K')), NASAPolynomial(coeffs=[15.5,0.0226347,-7.95274e-06,1.42522e-09,-1.00663e-13,11362.9,-48.8281], Tmin=(976.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C[C]=O(1863)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (210.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,375.832,3213.21],'cm^-1')),
        HinderedRotor(inertia=(0.0998903,'amu*angstrom^2'), symmetry=1, barrier=(10.0124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786918,'amu*angstrom^2'), symmetry=1, barrier=(78.8754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0998901,'amu*angstrom^2'), symmetry=1, barrier=(10.0124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0998911,'amu*angstrom^2'), symmetry=1, barrier=(10.0124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.66,'J/mol'), sigma=(6.12618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.10 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10296,0.0627311,-5.47256e-05,2.73666e-08,-5.73559e-12,25385.1,29.3539], Tmin=(100,'K'), Tmax=(1127.43,'K')), NASAPolynomial(coeffs=[10.1343,0.0306887,-1.20943e-05,2.15793e-09,-1.45717e-13,23348.6,-15.3004], Tmin=(1127.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=O(152)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08918,0.00200392,-1.61651e-05,2.55044e-08,-1.16417e-11,52802.7,4.52499], Tmin=(100,'K'), Tmax=(856.118,'K')), NASAPolynomial(coeffs=[0.961586,0.00569052,-3.48048e-06,7.19212e-10,-5.0805e-14,53738.7,21.4665], Tmin=(856.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C=CC[CH2](1193)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {5,S} {9,S}
4  C u1 p0 c0 {1,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (304.198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00401063,'amu*angstrom^2'), symmetry=1, barrier=(7.11751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.986217,'amu*angstrom^2'), symmetry=1, barrier=(22.6751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395285,'amu*angstrom^2'), symmetry=1, barrier=(22.6808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01079,0.0344277,7.32934e-06,-3.03227e-08,1.30186e-11,36666.4,21.0919], Tmin=(100,'K'), Tmax=(1010.22,'K')), NASAPolynomial(coeffs=[10.0278,0.0241677,-9.33579e-06,1.72607e-09,-1.22041e-13,33950.3,-23.0929], Tmin=(1010.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = 'CH2(T)(17)',
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
    label = '[CH]=CCC[C]=O(1503)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u1 p0 c0 {1,D} {3,S}
6  C u1 p0 c0 {4,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (283.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,216.955],'cm^-1')),
        HinderedRotor(inertia=(0.231983,'amu*angstrom^2'), symmetry=1, barrier=(7.73237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231444,'amu*angstrom^2'), symmetry=1, barrier=(7.73926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232248,'amu*angstrom^2'), symmetry=1, barrier=(7.73377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83123,0.0502801,-4.80732e-05,2.71327e-08,-6.57722e-12,34227.9,23.37], Tmin=(100,'K'), Tmax=(968.708,'K')), NASAPolynomial(coeffs=[7.48368,0.0269395,-1.19307e-05,2.25887e-09,-1.57784e-13,33132.8,-3.72001], Tmin=(968.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'O=C1CC=CCC1(1867)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {1,D} {2,S} {4,S}
6  C u0 p0 c0 {3,S} {7,D} {15,S}
7  C u0 p0 c0 {4,S} {6,D} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-134.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36165,0.0163311,8.69628e-05,-1.09507e-07,3.75179e-11,-16105.6,16.2397], Tmin=(100,'K'), Tmax=(1060.54,'K')), NASAPolynomial(coeffs=[10.6507,0.0379764,-1.84842e-05,3.81863e-09,-2.84935e-13,-20839.2,-38.2657], Tmin=(1060.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane)"""),
)

species(
    label = 'CC=CCC=C=O(2074)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  C u0 p0 c0 {2,S} {7,D} {15,S}
7  C u0 p0 c0 {1,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-33.7596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44166,0.0596875,-4.82585e-05,2.27303e-08,-4.73786e-12,-3971.11,23.0818], Tmin=(100,'K'), Tmax=(1085.85,'K')), NASAPolynomial(coeffs=[7.72778,0.0365311,-1.62705e-05,3.09114e-09,-2.16292e-13,-5336.28,-7.76287], Tmin=(1085.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.7596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=C=CCCC=O(1871)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {7,D} {12,S}
5  C u0 p0 c0 {1,D} {3,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (17.4773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5502,0.0570761,-3.98604e-05,1.49461e-08,-2.44499e-12,2187.17,24.0021], Tmin=(100,'K'), Tmax=(1332.53,'K')), NASAPolynomial(coeffs=[8.59202,0.0359378,-1.60653e-05,3.04136e-09,-2.11497e-13,310.493,-11.9921], Tmin=(1332.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.4773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=[C]CCC1[CH]C1(2075)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {4,S} {15,S}
7  C u1 p0 c0 {1,D} {5,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (264.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32932,0.0511673,-2.57647e-05,2.88401e-09,1.06921e-12,31946.1,28.1619], Tmin=(100,'K'), Tmax=(1236.86,'K')), NASAPolynomial(coeffs=[11.2168,0.0294895,-1.19641e-05,2.17712e-09,-1.48535e-13,28712.5,-24.8255], Tmin=(1236.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CCC1=O(2076)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {3,S} {13,S}
6  C u0 p0 c0 {1,D} {2,S} {4,S}
7  C u1 p0 c0 {2,S} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (166.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62976,0.0377604,2.73272e-05,-5.93546e-08,2.44918e-11,20081.6,24.4922], Tmin=(100,'K'), Tmax=(985.077,'K')), NASAPolynomial(coeffs=[12.6278,0.0278107,-1.03746e-05,1.92936e-09,-1.38897e-13,16230.8,-36.9496], Tmin=(985.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][CH]C1CCC1=O(2077)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {1,D} {2,S} {4,S}
6  C u1 p0 c0 {2,S} {7,S} {13,S}
7  C u1 p0 c0 {6,S} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (235.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82985,0.0340223,3.10139e-05,-5.91234e-08,2.33937e-11,28462.2,27.189], Tmin=(100,'K'), Tmax=(1001.14,'K')), NASAPolynomial(coeffs=[11.3151,0.0292285,-1.14026e-05,2.15007e-09,-1.54794e-13,24904,-26.868], Tmin=(1001.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = 'CO(14)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-118.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2193.04],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C=C[CH2](74)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,D} {3,S} {5,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u1 p0 c0 {1,S} {7,S} {10,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.209873,'amu*angstrom^2'), symmetry=1, barrier=(25.2262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778586,'amu*angstrom^2'), symmetry=1, barrier=(93.3137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56317,0.0223429,1.87066e-05,-3.93097e-08,1.63982e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.265,'K')), NASAPolynomial(coeffs=[9.82996,0.0151966,-5.22271e-06,9.67654e-10,-7.07858e-14,30607.7,-26.985], Tmin=(974.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'H(6)',
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
    label = '[CH2]C=CCC=C=O(2078)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  C u0 p0 c0 {1,D} {5,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (117.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.601726,'amu*angstrom^2'), symmetry=1, barrier=(13.8349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0466388,'amu*angstrom^2'), symmetry=1, barrier=(9.36357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02553,'amu*angstrom^2'), symmetry=1, barrier=(23.5789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23481,0.0597159,-4.88265e-05,2.14386e-08,-3.91701e-12,14261.1,24.0971], Tmin=(100,'K'), Tmax=(1275.07,'K')), NASAPolynomial(coeffs=[11.1768,0.028527,-1.21357e-05,2.25487e-09,-1.55699e-13,11725.8,-26.2831], Tmin=(1275.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CCC[C]=O(2079)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {6,D} {13,S} {14,S}
6  C u0 p0 c0 {4,D} {5,D}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (177.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,193.747,193.812],'cm^-1')),
        HinderedRotor(inertia=(0.396092,'amu*angstrom^2'), symmetry=1, barrier=(10.5063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396602,'amu*angstrom^2'), symmetry=1, barrier=(10.5067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396334,'amu*angstrom^2'), symmetry=1, barrier=(10.5058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3804,0.0609701,-5.7359e-05,3.1645e-08,-7.52254e-12,21432.4,25.447], Tmin=(100,'K'), Tmax=(983.61,'K')), NASAPolynomial(coeffs=[8.23605,0.0330905,-1.48426e-05,2.82834e-09,-1.98312e-13,20083.7,-7.51422], Tmin=(983.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C=C(1759)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 C u2 p0 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,227.33,229.304,232.839],'cm^-1')),
        HinderedRotor(inertia=(1.31513,'amu*angstrom^2'), symmetry=1, barrier=(50.4662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817973,3.34731e-05,-4.36187e-08,1.5821e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.759,'K')), NASAPolynomial(coeffs=[5.3676,0.0170742,-6.35103e-06,1.16618e-09,-8.27612e-14,44095,-3.44632], Tmin=(983.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C[C]=O(169)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {7,S} {8,S}
4 C u1 p0 c0 {1,D} {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (167.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.182442,'amu*angstrom^2'), symmetry=1, barrier=(4.19471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180479,'amu*angstrom^2'), symmetry=1, barrier=(4.14956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0632,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63264,0.0331121,-4.63435e-05,4.02609e-08,-1.40729e-11,20157.3,15.987], Tmin=(100,'K'), Tmax=(833.292,'K')), NASAPolynomial(coeffs=[4.91249,0.0168159,-7.3739e-06,1.37542e-09,-9.40241e-14,19963.2,6.51888], Tmin=(833.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=CCC=C[O](2080)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,D} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {4,D} {15,S}
7  C u1 p0 c0 {5,S} {13,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (136.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949904,0.0501068,1.01319e-05,-5.57758e-08,2.67744e-11,16589.7,27.1271], Tmin=(100,'K'), Tmax=(955.03,'K')), NASAPolynomial(coeffs=[18.4852,0.0184132,-5.66341e-06,1.02494e-09,-7.69619e-14,11336.4,-66.6317], Tmin=(955.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=CCC[C]=O(2081)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u1 p0 c0 {4,S} {5,D}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (238.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78564,0.0576868,-2.65034e-05,-5.08729e-08,6.16772e-11,28776.6,25.5507], Tmin=(100,'K'), Tmax=(485.475,'K')), NASAPolynomial(coeffs=[5.89027,0.0394065,-1.80343e-05,3.42915e-09,-2.38528e-13,28194.9,6.82838], Tmin=(485.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=C[CH]CC=O(2082)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u1 p0 c0 {2,S} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {7,S} {11,S}
6  C u0 p0 c0 {1,D} {2,S} {13,S}
7  C u1 p0 c0 {5,S} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (192.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,346.838,1871.44],'cm^-1')),
        HinderedRotor(inertia=(0.108182,'amu*angstrom^2'), symmetry=1, barrier=(9.2361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108196,'amu*angstrom^2'), symmetry=1, barrier=(9.2361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.320181,'amu*angstrom^2'), symmetry=1, barrier=(27.331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010997,'amu*angstrom^2'), symmetry=1, barrier=(27.3311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05788,0.0569622,-3.72418e-05,1.21239e-08,-1.59597e-12,23237.5,29.5197], Tmin=(100,'K'), Tmax=(1748.51,'K')), NASAPolynomial(coeffs=[14.8761,0.0253504,-1.01224e-05,1.78372e-09,-1.17528e-13,18405.4,-44.8661], Tmin=(1748.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC[C]=O(2083)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u1 p0 c0 {2,S} {5,D}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (238.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.262153,'amu*angstrom^2'), symmetry=1, barrier=(6.0274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262177,'amu*angstrom^2'), symmetry=1, barrier=(6.02796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262572,'amu*angstrom^2'), symmetry=1, barrier=(6.03706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262163,'amu*angstrom^2'), symmetry=1, barrier=(6.02764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78564,0.0576868,-2.65034e-05,-5.08729e-08,6.16772e-11,28776.6,25.5507], Tmin=(100,'K'), Tmax=(485.475,'K')), NASAPolynomial(coeffs=[5.89027,0.0394065,-1.80343e-05,3.42915e-09,-2.38528e-13,28194.9,6.82838], Tmin=(485.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=[C]CCC=O(2084)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {7,D} {12,S}
5  C u0 p0 c0 {1,D} {2,S} {13,S}
6  C u1 p0 c0 {4,S} {14,S} {15,S}
7  C u1 p0 c0 {3,S} {4,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (230.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1685,370,244.072,244.073],'cm^-1')),
        HinderedRotor(inertia=(0.113995,'amu*angstrom^2'), symmetry=1, barrier=(4.81909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113999,'amu*angstrom^2'), symmetry=1, barrier=(4.81909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00282978,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639021,'amu*angstrom^2'), symmetry=1, barrier=(27.0129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49274,0.0592249,-4.5467e-05,1.99225e-08,-3.89795e-12,27775.1,25.9199], Tmin=(100,'K'), Tmax=(1134.35,'K')), NASAPolynomial(coeffs=[7.56517,0.0378118,-1.71511e-05,3.28072e-09,-2.3022e-13,26397.5,-4.14143], Tmin=(1134.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]C=CC[C]=O(2085)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u1 p0 c0 {3,S} {6,S} {13,S}
6  C u0 p0 c0 {4,D} {5,S} {15,S}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (149.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,829.985,830.307],'cm^-1')),
        HinderedRotor(inertia=(0.72995,'amu*angstrom^2'), symmetry=1, barrier=(16.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147487,'amu*angstrom^2'), symmetry=1, barrier=(3.39102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343842,'amu*angstrom^2'), symmetry=1, barrier=(16.7889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150487,'amu*angstrom^2'), symmetry=1, barrier=(16.7852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36151,0.0467316,-3.9003e-06,-2.24189e-08,9.72248e-12,18056.7,25.5773], Tmin=(100,'K'), Tmax=(1120.66,'K')), NASAPolynomial(coeffs=[12.9783,0.0295868,-1.35031e-05,2.65778e-09,-1.91472e-13,13925.8,-38.6038], Tmin=(1120.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CCCC=O(2086)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,D} {12,S}
5  C u0 p0 c0 {1,D} {3,S} {13,S}
6  C u1 p0 c0 {7,S} {14,S} {15,S}
7  C u1 p0 c0 {4,D} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (230.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1685,370,244.072,244.073],'cm^-1')),
        HinderedRotor(inertia=(0.113995,'amu*angstrom^2'), symmetry=1, barrier=(4.81909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113999,'amu*angstrom^2'), symmetry=1, barrier=(4.81909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00282978,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639021,'amu*angstrom^2'), symmetry=1, barrier=(27.0129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49274,0.0592249,-4.5467e-05,1.99225e-08,-3.89795e-12,27775.1,25.9199], Tmin=(100,'K'), Tmax=(1134.35,'K')), NASAPolynomial(coeffs=[7.56517,0.0378118,-1.71511e-05,3.28072e-09,-2.3022e-13,26397.5,-4.14143], Tmin=(1134.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=CC[CH][C]=O(2087)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {14,S}
5  C u0 p0 c0 {3,S} {4,D} {13,S}
6  C u1 p0 c0 {2,S} {7,S} {15,S}
7  C u1 p0 c0 {1,D} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (168.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,255.59,256.392],'cm^-1')),
        HinderedRotor(inertia=(0.00245468,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173277,'amu*angstrom^2'), symmetry=1, barrier=(7.91676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248567,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325739,'amu*angstrom^2'), symmetry=1, barrier=(14.4668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24409,0.058217,-4.29307e-05,1.72312e-08,-2.90965e-12,20350.4,26.6305], Tmin=(100,'K'), Tmax=(1360.1,'K')), NASAPolynomial(coeffs=[10.6868,0.0304464,-1.23036e-05,2.21897e-09,-1.50258e-13,17781.8,-21.8292], Tmin=(1360.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC1CCC1=O(2088)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {1,D} {2,S} {4,S}
6  C u0 p0 c0 {2,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-42.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84948,0.0343186,2.85892e-05,-5.46711e-08,2.12631e-11,-5035.21,22.6428], Tmin=(100,'K'), Tmax=(1016.73,'K')), NASAPolynomial(coeffs=[10.6265,0.0309466,-1.24038e-05,2.34865e-09,-1.68421e-13,-8430.47,-27.7667], Tmin=(1016.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone)"""),
)

species(
    label = 'C=CC=CCC=O(2089)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {7,D} {12,S}
6  C u0 p0 c0 {1,D} {2,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-41.3668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25993,0.0458328,7.57118e-06,-3.98008e-08,1.69909e-11,-4864.2,24.141], Tmin=(100,'K'), Tmax=(1049.02,'K')), NASAPolynomial(coeffs=[14.9619,0.0264618,-1.17389e-05,2.34545e-09,-1.72885e-13,-9547.83,-51.2405], Tmin=(1049.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.3668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCCC=C=O(2090)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u0 p0 c0 {1,D} {5,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-21.4616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06199,0.0573791,-3.83733e-05,1.29098e-08,-1.76247e-12,-2469.58,26.3395], Tmin=(100,'K'), Tmax=(1688.18,'K')), NASAPolynomial(coeffs=[14.3948,0.0257886,-1.03043e-05,1.82539e-09,-1.21016e-13,-6971.25,-44.9654], Tmin=(1688.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.4616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C1C[CH][CH]CC1(2091)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {1,D} {2,S} {4,S}
6  C u1 p0 c0 {3,S} {7,S} {14,S}
7  C u1 p0 c0 {4,S} {6,S} {15,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (129.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04364,0.0342051,1.86383e-05,-3.43503e-08,1.12686e-11,15665.6,22.7044], Tmin=(100,'K'), Tmax=(1203.9,'K')), NASAPolynomial(coeffs=[8.54985,0.0380917,-1.79809e-05,3.52442e-09,-2.50525e-13,12250.8,-17.5675], Tmin=(1203.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclohexanone) + radical(cyclohexane) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CC=CC[C]=O(2092)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (118.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.78347,'amu*angstrom^2'), symmetry=1, barrier=(18.0135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782872,'amu*angstrom^2'), symmetry=1, barrier=(17.9998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782314,'amu*angstrom^2'), symmetry=1, barrier=(17.9869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19677,0.0487089,-7.28985e-06,-2.55703e-08,1.26837e-11,14375.7,25.1849], Tmin=(100,'K'), Tmax=(1042.88,'K')), NASAPolynomial(coeffs=[15.5436,0.021984,-9.5591e-06,1.90322e-09,-1.40485e-13,9844.15,-52.0117], Tmin=(1042.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC[CH]C[C]=O(2056)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {3,S} {12,S}
5  C u0 p0 c0 {2,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (212.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,226.256,1084.44,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00329293,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237761,'amu*angstrom^2'), symmetry=1, barrier=(8.63669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499253,'amu*angstrom^2'), symmetry=1, barrier=(18.1366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499293,'amu*angstrom^2'), symmetry=1, barrier=(18.1366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32543,0.0579883,-4.30361e-05,1.72033e-08,-2.90315e-12,25712.4,29.5118], Tmin=(100,'K'), Tmax=(1351.01,'K')), NASAPolynomial(coeffs=[10.4471,0.0309813,-1.30507e-05,2.40681e-09,-1.65103e-13,23247.7,-17.2395], Tmin=(1351.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]CCC[C]=O(2093)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {6,D} {14,S} {15,S}
6  C u1 p0 c0 {3,S} {5,D}
7  C u1 p0 c0 {1,D} {4,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (250.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,325.917,325.918,3137.29],'cm^-1')),
        HinderedRotor(inertia=(0.141606,'amu*angstrom^2'), symmetry=1, barrier=(10.6739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17403,'amu*angstrom^2'), symmetry=1, barrier=(13.1179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141605,'amu*angstrom^2'), symmetry=1, barrier=(10.6739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141605,'amu*angstrom^2'), symmetry=1, barrier=(10.6739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17595,0.0666815,-7.18026e-05,4.93287e-08,-1.47947e-11,30276.4,28.0413], Tmin=(100,'K'), Tmax=(790.439,'K')), NASAPolynomial(coeffs=[6.86186,0.0379079,-1.71991e-05,3.27508e-09,-2.28772e-13,29377.5,1.94725], Tmin=(790.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CCC[CH][C]=O(2094)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u1 p0 c0 {3,S} {7,S} {13,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (180.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,298.472,298.617,1488.11],'cm^-1')),
        HinderedRotor(inertia=(0.000978494,'amu*angstrom^2'), symmetry=1, barrier=(10.9569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255513,'amu*angstrom^2'), symmetry=1, barrier=(16.1517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336858,'amu*angstrom^2'), symmetry=1, barrier=(21.355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43577,'amu*angstrom^2'), symmetry=1, barrier=(90.7087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01996,0.0600114,-4.50256e-05,1.80224e-08,-2.96958e-12,21833.9,27.8781], Tmin=(100,'K'), Tmax=(1418.62,'K')), NASAPolynomial(coeffs=[12.5362,0.0275399,-1.06914e-05,1.88736e-09,-1.26149e-13,18566.5,-31.708], Tmin=(1418.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CCCC[C]=O(2095)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u1 p0 c0 {1,D} {4,S}
7  C u1 p0 c0 {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (260.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21383,0.0647017,-6.05903e-05,3.34196e-08,-7.94899e-12,31389.3,27.8331], Tmin=(100,'K'), Tmax=(983.083,'K')), NASAPolynomial(coeffs=[8.44359,0.0352851,-1.57061e-05,2.98191e-09,-2.08643e-13,29967.8,-6.923], Tmin=(983.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=CCCC=O(2096)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {7,S} {14,S}
6  C u0 p0 c0 {1,D} {3,S} {13,S}
7  C u2 p0 c0 {5,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (211.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47633,0.0579113,-3.26975e-05,8.88337e-09,-9.95681e-13,25531.9,26.307], Tmin=(100,'K'), Tmax=(1878.3,'K')), NASAPolynomial(coeffs=[11.5542,0.0364497,-1.55584e-05,2.80023e-09,-1.86024e-13,21746.1,-28.6657], Tmin=(1878.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-28.4299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (569.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (217.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (162.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (562.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (484.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-20.8987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-3.45664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-3.45664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (197.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (28.9227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (55.1352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (35.2167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (67.9345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (149.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (225.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (254.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (363.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (128.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (179.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (204.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (219.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (150.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (386.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (93.7491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (76.8241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-20.5221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (49.8172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (60.5387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (31.2835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (94.7023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (152.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (157.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (270.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (124.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (239.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (55.9322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['CH2CO(27)', 'C=CC=C(61)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH2]C=CC[C]=O(1300)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(148.957,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['[CH2]C=CCC(=C)[O](1858)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C[C]=O(1863)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(133.469,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=O(152)', '[CH2]C=CC[CH2](1193)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(17)', '[CH]=CCC[C]=O(1503)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['O=C1CC=CCC1(1867)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['CC=CCC=C=O(2074)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['C=C=CCCC=O(1871)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['O=[C]CCC1[CH]C1(2075)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['[CH2]C1[CH]CCC1=O(2076)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['[CH2][CH]C1CCC1=O(2077)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(83.5651,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 79.7 to 83.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CO(14)', '[CH2]C=CC[CH2](1193)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(30.522,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2CO(27)', '[CH2]C=C[CH2](74)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.002829,'m^3/(mol*s)'), n=2.61489, Ea=(34.8017,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4622630802209839, var=1.4932365737855162, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_N-Sp-5R!H-4R!H_Ext-2R!H-R_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_N-Sp-5R!H-4R!H_Ext-2R!H-R_Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(6)', '[CH2]C=CCC=C=O(2078)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(0.885011,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(6)', 'C=C=CCC[C]=O(2079)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.08158e+21,'m^3/(mol*s)'), n=-4.30708, Ea=(17.4721,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]=O(167)', '[CH2]C=C[CH2](74)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.0455e+08,'m^3/(mol*s)'), n=-0.533333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.09789675142796564, var=0.2171391628300253, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_3BrCO->O_Ext-2CF-R_N-5R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_3BrCO->O_Ext-2CF-R_N-5R!H->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C=C(1759)', '[CH2]C[C]=O(169)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['[CH2]C=CCC=C[O](2080)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['C[C]=CCC[C]=O(2081)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C=C[CH]CC=O(2082)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.32587e+10,'s^-1'), n=0.723333, Ea=(193.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/Cd;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC=[C]CC[C]=O(2083)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C=[C]CCC=O(2084)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]C=CC[C]=O(2085)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=CCCC=O(2086)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;CO_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC=CC[CH][C]=O(2087)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['C=CC1CCC1=O(2088)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['C=CC=CCC=O(2089)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['C=CCCC=C=O(2090)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['O=C1C[CH][CH]CC1(2091)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.25216e+09,'s^-1'), n=0.38532, Ea=(59.7134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=O(167)', 'C=CC=C(61)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(18.8247,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(6)', 'C=CC=CC[C]=O(2092)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.667852,'m^3/(mol*s)'), n=2.28533, Ea=(2.61635,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.25112379235625565, var=1.3735217093431273, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC[CH]C[C]=O(2056)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]CCC[C]=O(2093)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=CCC[CH][C]=O(2094)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C=CCC[C]=O(1862)'],
    products = ['[CH]=CCCC[C]=O(2095)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C=CCCC=O(2096)'],
    products = ['[CH2]C=CCC[C]=O(1862)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.546e+09,'s^-1'), n=0.732, Ea=(25.1375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1098',
    isomers = [
        '[CH2]C=CCC[C]=O(1862)',
    ],
    reactants = [
        ('CH2CO(27)', 'C=CC=C(61)'),
        ('[CH2][C]=O(167)', 'C=CC=C(61)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1098',
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

