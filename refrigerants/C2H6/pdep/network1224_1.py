species(
    label = 'C=[C]C(C)C(=C)[O](2142)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (216.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,187.843,187.845,187.845],'cm^-1')),
        HinderedRotor(inertia=(0.00477797,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406476,'amu*angstrom^2'), symmetry=1, barrier=(10.1774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406478,'amu*angstrom^2'), symmetry=1, barrier=(10.1774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86918,0.0689848,-6.69141e-05,3.62265e-08,-8.08248e-12,26123.7,25.6832], Tmin=(100,'K'), Tmax=(1071.36,'K')), NASAPolynomial(coeffs=[11.3137,0.0299888,-1.23156e-05,2.25161e-09,-1.54397e-13,23885.8,-25.4256], Tmin=(1071.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
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
    label = 'C=C=CC(1692)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.759585,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74634,0.0218191,8.22302e-06,-2.14762e-08,8.55596e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[6.82084,0.0192337,-7.45616e-06,1.36535e-09,-9.53184e-14,16028,-10.4336], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C]CC(=C)[O](1431)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (249.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,400.022,400.338,400.737],'cm^-1')),
        HinderedRotor(inertia=(0.0760986,'amu*angstrom^2'), symmetry=1, barrier=(8.67354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761993,'amu*angstrom^2'), symmetry=1, barrier=(8.67861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.78,'J/mol'), sigma=(6.06207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.87 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54215,0.0477416,-3.21904e-05,6.04368e-09,1.69299e-12,30142.6,23.6772], Tmin=(100,'K'), Tmax=(1008.06,'K')), NASAPolynomial(coeffs=[12.024,0.0183644,-6.65275e-06,1.17509e-09,-8.06673e-14,27408.7,-30.0539], Tmin=(1008.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)=CC(=C)[O](2671)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (60.8083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753118,0.06288,-3.80581e-05,6.82713e-10,5.72964e-12,7438.14,24.513], Tmin=(100,'K'), Tmax=(956.44,'K')), NASAPolynomial(coeffs=[14.7682,0.0240138,-8.07404e-06,1.37016e-09,-9.264e-14,3853.99,-47.1993], Tmin=(956.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.8083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(=CC)C(=C)[O](2421)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,D} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u1 p0 c0 {3,S} {14,S} {15,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (95.7161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462914,0.0694054,-5.31693e-05,1.38006e-08,1.49844e-12,11646.8,24.1211], Tmin=(100,'K'), Tmax=(975.079,'K')), NASAPolynomial(coeffs=[16.3564,0.0219777,-7.54728e-06,1.29949e-09,-8.85418e-14,7702.54,-56.4864], Tmin=(975.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.7161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=C([O])CC(2672)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {7,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (157.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634257,0.0661591,-4.65674e-05,7.41842e-09,4.04234e-12,19070.2,25.9488], Tmin=(100,'K'), Tmax=(939.612,'K')), NASAPolynomial(coeffs=[15.4472,0.0222825,-7.14721e-06,1.17772e-09,-7.84531e-14,15439.7,-49.0988], Tmin=(939.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[C]=C(57)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u2 p0 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.36689e-08,8.45455e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.85899e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([O])=CC(574)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {5,S}
5  C u1 p0 c0 {4,S} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (52.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,270.949],'cm^-1')),
        HinderedRotor(inertia=(0.126895,'amu*angstrom^2'), symmetry=1, barrier=(6.6168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126614,'amu*angstrom^2'), symmetry=1, barrier=(6.61684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96215,0.0395246,-2.15042e-05,-2.39262e-09,4.75075e-12,6363.74,18.5338], Tmin=(100,'K'), Tmax=(938.6,'K')), NASAPolynomial(coeffs=[10.4572,0.016349,-5.28591e-06,8.75308e-10,-5.83597e-14,4195.22,-24.9683], Tmin=(938.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'O(7)',
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
    label = 'C=[C]C(C)[C]=C(2528)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {6,D} {12,S} {13,S}
4  C u0 p0 c0 {5,D} {11,S} {14,S}
5  C u1 p0 c0 {1,S} {4,D}
6  C u1 p0 c0 {1,S} {3,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (530.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,246.5,248.077],'cm^-1')),
        HinderedRotor(inertia=(0.228446,'amu*angstrom^2'), symmetry=1, barrier=(10.0527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233159,'amu*angstrom^2'), symmetry=1, barrier=(10.0708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237074,'amu*angstrom^2'), symmetry=1, barrier=(10.0758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55319,0.0577252,-5.7179e-05,3.71685e-08,-1.0816e-11,63932.8,21.4596], Tmin=(100,'K'), Tmax=(804.847,'K')), NASAPolynomial(coeffs=[6.00243,0.0356122,-1.59653e-05,3.02938e-09,-2.11427e-13,63216.6,0.960658], Tmin=(804.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC(=C)C1C(2424)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (38.1574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60154,0.0379193,2.65896e-05,-5.88323e-08,2.42081e-11,4688.75,19.851], Tmin=(100,'K'), Tmax=(994.35,'K')), NASAPolynomial(coeffs=[13.2545,0.026702,-1.02814e-05,1.95357e-09,-1.42303e-13,608.473,-45.166], Tmin=(994.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.1574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C=C(C)C(=C)O(2673)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {15,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {7,D}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
"""),
    E0 = (-16.9836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299673,0.0805484,-8.03292e-05,3.98939e-08,-7.60249e-12,-1875.92,24.1366], Tmin=(100,'K'), Tmax=(1385.52,'K')), NASAPolynomial(coeffs=[20.5919,0.0151854,-4.09927e-06,5.84416e-10,-3.50009e-14,-7180.42,-81.7162], Tmin=(1385.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.9836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=[C]C(C)[C]1CO1(2674)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (362.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73001,0.0647025,-5.79156e-05,3.03657e-08,-6.38446e-12,43781,27.0386], Tmin=(100,'K'), Tmax=(1278.51,'K')), NASAPolynomial(coeffs=[11.3549,0.0265457,-7.38171e-06,1.00822e-09,-5.58915e-14,41465.9,-25.2599], Tmin=(1278.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]1CC(=O)C1C(2575)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {3,S} {7,S}
6  C u0 p0 c0 {1,D} {2,S} {3,S}
7  C u1 p0 c0 {5,S} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (212.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65537,0.0455706,-1.70821e-05,-5.91173e-10,1.31111e-12,25612.6,26.1822], Tmin=(100,'K'), Tmax=(1329.87,'K')), NASAPolynomial(coeffs=[8.55271,0.0334586,-1.31591e-05,2.32412e-09,-1.54674e-13,23014.6,-11.9302], Tmin=(1329.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([O])C(=C)C1C(2675)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {3,S} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (376.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0412,0.0529442,-9.9235e-06,-2.65901e-08,1.44377e-11,45353.1,24.6002], Tmin=(100,'K'), Tmax=(982.168,'K')), NASAPolynomial(coeffs=[14.965,0.0242002,-8.72979e-06,1.58664e-09,-1.12734e-13,41269.3,-49.1897], Tmin=(982.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'CH3(18)',
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
    label = 'C=C=CC(=C)[O](1143)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (158.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1309,'amu*angstrom^2'), symmetry=1, barrier=(26.0016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976079,0.0558196,-5.47944e-05,2.68999e-08,-5.05266e-12,19170.5,20.4568], Tmin=(100,'K'), Tmax=(1432.3,'K')), NASAPolynomial(coeffs=[15.4408,0.0107351,-2.66835e-06,3.52119e-10,-1.99511e-14,15507.9,-52.8449], Tmin=(1432.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
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
    label = 'C=C=C(C)C(=C)[O](2676)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {7,D}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (120.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.873377,'amu*angstrom^2'), symmetry=1, barrier=(20.0807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87168,'amu*angstrom^2'), symmetry=1, barrier=(20.0416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443575,0.0722261,-7.30793e-05,3.84332e-08,-7.96348e-12,14664.6,22.964], Tmin=(100,'K'), Tmax=(1180.19,'K')), NASAPolynomial(coeffs=[15.7366,0.0203932,-7.20011e-06,1.21908e-09,-8.0334e-14,11054.9,-53.3497], Tmin=(1180.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=CC(1693)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u1 p0 c0 {4,S} {9,S} {10,S}
4  C u1 p0 c0 {2,D} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12014e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.76,'K')), NASAPolynomial(coeffs=[10.7464,0.0143241,-5.20138e-06,8.69083e-10,-5.48388e-14,40045.6,-31.3798], Tmin=(2065.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C(=C)[O](2677)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (119.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.781168,'amu*angstrom^2'), symmetry=1, barrier=(17.9606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.780627,'amu*angstrom^2'), symmetry=1, barrier=(17.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781673,'amu*angstrom^2'), symmetry=1, barrier=(17.9722,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478971,0.0716402,-7.27903e-05,3.85345e-08,-8.03983e-12,14561.6,23.9386], Tmin=(100,'K'), Tmax=(1172.3,'K')), NASAPolynomial(coeffs=[15.5148,0.0203372,-7.1471e-06,1.20484e-09,-7.91444e-14,11036.3,-50.9909], Tmin=(1172.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C=C(C)C(=C)[O](2141)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (62.2727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,364.236,364.679],'cm^-1')),
        HinderedRotor(inertia=(0.159215,'amu*angstrom^2'), symmetry=1, barrier=(15.0636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160502,'amu*angstrom^2'), symmetry=1, barrier=(15.0553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856122,'amu*angstrom^2'), symmetry=1, barrier=(80.7861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753129,0.0628798,-3.80576e-05,6.8207e-10,5.72993e-12,7614.26,24.357], Tmin=(100,'K'), Tmax=(956.437,'K')), NASAPolynomial(coeffs=[14.7682,0.0240139,-8.07407e-06,1.37016e-09,-9.26407e-14,4030.13,-47.3551], Tmin=(956.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.2727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC(C)C(=C)[O](2146)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  C u1 p0 c0 {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (225.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,250.313,251.622],'cm^-1')),
        HinderedRotor(inertia=(0.241225,'amu*angstrom^2'), symmetry=1, barrier=(10.6174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237104,'amu*angstrom^2'), symmetry=1, barrier=(10.6083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235389,'amu*angstrom^2'), symmetry=1, barrier=(10.5944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678248,0.0697169,-6.52479e-05,3.28424e-08,-6.68525e-12,27246.5,26.2946], Tmin=(100,'K'), Tmax=(1181.37,'K')), NASAPolynomial(coeffs=[13.4569,0.0264496,-1.0311e-05,1.84062e-09,-1.24692e-13,24227.2,-37.4849], Tmin=(1181.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)[O](1859)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {6,D} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u1 p0 c0 {2,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {12,S} {13,S}
7  C u0 p0 c0 {4,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (183.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,323.791,331.341,2023.52],'cm^-1')),
        HinderedRotor(inertia=(0.166079,'amu*angstrom^2'), symmetry=1, barrier=(13.0535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172042,'amu*angstrom^2'), symmetry=1, barrier=(13.0533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17049,'amu*angstrom^2'), symmetry=1, barrier=(12.9812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.52,'J/mol'), sigma=(6.40826,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.29 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648497,0.0693305,-6.62579e-05,3.47966e-08,-7.35886e-12,22195.4,27.7785], Tmin=(100,'K'), Tmax=(1145.9,'K')), NASAPolynomial(coeffs=[13.1664,0.0256346,-9.06e-06,1.52019e-09,-9.90665e-14,19326.5,-34.3182], Tmin=(1145.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]=C(C)C(=C)O(2678)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {15,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {7,D}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  C u1 p0 c0 {7,S} {13,S} {14,S}
7  C u1 p0 c0 {3,D} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
"""),
    E0 = (162.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.111007,0.0768382,-7.37003e-05,3.60251e-08,-6.803e-12,19680.9,26.4368], Tmin=(100,'K'), Tmax=(1403.4,'K')), NASAPolynomial(coeffs=[18.8766,0.0175565,-4.82001e-06,6.83001e-10,-4.0235e-14,14859.9,-69.7908], Tmin=(1403.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(C)[C]=C(2679)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {14,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u0 p0 c0 {6,D} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {5,D}
7  C u1 p0 c0 {4,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (325.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.64381,'amu*angstrom^2'), symmetry=1, barrier=(14.8024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643394,'amu*angstrom^2'), symmetry=1, barrier=(14.7929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64351,'amu*angstrom^2'), symmetry=1, barrier=(14.7956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64364,'amu*angstrom^2'), symmetry=1, barrier=(14.7986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266278,0.0784406,-8.30582e-05,4.57659e-08,-9.96865e-12,39293.4,26.746], Tmin=(100,'K'), Tmax=(1120.93,'K')), NASAPolynomial(coeffs=[15.8404,0.0228646,-8.68712e-06,1.53386e-09,-1.03542e-13,35801.9,-50.1679], Tmin=(1120.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([C]=C)C(=C)O(2145)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {15,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
"""),
    E0 = (283.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,199.996,200.025],'cm^-1')),
        HinderedRotor(inertia=(0.543499,'amu*angstrom^2'), symmetry=1, barrier=(15.4305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00421386,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543459,'amu*angstrom^2'), symmetry=1, barrier=(15.4303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.543356,'amu*angstrom^2'), symmetry=1, barrier=(15.4299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23495,0.0780817,-8.41943e-05,4.79074e-08,-1.07256e-11,34242.3,28.2349], Tmin=(100,'K'), Tmax=(1096.62,'K')), NASAPolynomial(coeffs=[15.6371,0.0219023,-7.35156e-06,1.19351e-09,-7.6266e-14,30864.1,-47.4925], Tmin=(1096.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C(C)C=C(2147)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  C u1 p0 c0 {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (225.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,250.772,250.804],'cm^-1')),
        HinderedRotor(inertia=(0.237638,'amu*angstrom^2'), symmetry=1, barrier=(10.6065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237627,'amu*angstrom^2'), symmetry=1, barrier=(10.6068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237703,'amu*angstrom^2'), symmetry=1, barrier=(10.6068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678248,0.0697169,-6.52479e-05,3.28424e-08,-6.68525e-12,27246.5,26.2946], Tmin=(100,'K'), Tmax=(1181.37,'K')), NASAPolynomial(coeffs=[13.4569,0.0264496,-1.0311e-05,1.84062e-09,-1.24692e-13,24227.2,-37.4849], Tmin=(1181.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(C)C(=C)O(2680)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {14,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {7,D}
7  C u1 p0 c0 {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (325.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.64381,'amu*angstrom^2'), symmetry=1, barrier=(14.8024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643394,'amu*angstrom^2'), symmetry=1, barrier=(14.7929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64351,'amu*angstrom^2'), symmetry=1, barrier=(14.7956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64364,'amu*angstrom^2'), symmetry=1, barrier=(14.7986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266278,0.0784406,-8.30582e-05,4.57659e-08,-9.96865e-12,39293.4,26.746], Tmin=(100,'K'), Tmax=(1120.93,'K')), NASAPolynomial(coeffs=[15.8404,0.0228646,-8.68712e-06,1.53386e-09,-1.03542e-13,35801.9,-50.1679], Tmin=(1120.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(C)C[C]=O(2069)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {6,D} {14,S} {15,S}
6  C u1 p0 c0 {2,S} {5,D}
7  C u1 p0 c0 {1,D} {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (242.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,261.734,261.735],'cm^-1')),
        HinderedRotor(inertia=(0.00246084,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172416,'amu*angstrom^2'), symmetry=1, barrier=(8.38182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172415,'amu*angstrom^2'), symmetry=1, barrier=(8.38185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172422,'amu*angstrom^2'), symmetry=1, barrier=(8.38185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.66,'J/mol'), sigma=(6.12618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.10 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27872,0.0628602,-5.68104e-05,3.03764e-08,-7.03817e-12,29315.5,27.4237], Tmin=(100,'K'), Tmax=(1004.69,'K')), NASAPolynomial(coeffs=[8.24352,0.0351311,-1.54108e-05,2.90554e-09,-2.02504e-13,27916,-6.21003], Tmin=(1004.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CC(=C)OC(2681)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {7,D} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u1 p0 c0 {7,S} {14,S} {15,S}
7  C u1 p0 c0 {4,D} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (218.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,361.375,361.428,361.501],'cm^-1')),
        HinderedRotor(inertia=(0.991971,'amu*angstrom^2'), symmetry=1, barrier=(91.9176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1971,'amu*angstrom^2'), symmetry=1, barrier=(18.2632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196954,'amu*angstrom^2'), symmetry=1, barrier=(18.2624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197072,'amu*angstrom^2'), symmetry=1, barrier=(18.2624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648507,0.0635621,-3.57576e-05,-4.05424e-09,7.72077e-12,26419.4,24.7393], Tmin=(100,'K'), Tmax=(966.179,'K')), NASAPolynomial(coeffs=[16.227,0.0222364,-7.57007e-06,1.31638e-09,-9.11192e-14,22327.7,-55.4784], Tmin=(966.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = 'C=[C]C(C)[C]=O(2654)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (266.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,302.177],'cm^-1')),
        HinderedRotor(inertia=(0.129157,'amu*angstrom^2'), symmetry=1, barrier=(8.27106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12744,'amu*angstrom^2'), symmetry=1, barrier=(8.28055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128227,'amu*angstrom^2'), symmetry=1, barrier=(8.30737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94246,0.0487255,-4.70948e-05,2.88356e-08,-7.86968e-12,32127.7,22.577], Tmin=(100,'K'), Tmax=(853.763,'K')), NASAPolynomial(coeffs=[6.02769,0.0295858,-1.3468e-05,2.57823e-09,-1.81031e-13,31430.1,3.51394], Tmin=(853.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=C1CC(=O)C1C(2582)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {1,D} {2,S} {3,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-36.2847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77167,0.0326478,4.07371e-05,-6.90751e-08,2.59325e-11,-4269.69,19.8406], Tmin=(100,'K'), Tmax=(1041.38,'K')), NASAPolynomial(coeffs=[12.8946,0.0301197,-1.35191e-05,2.7231e-09,-2.01773e-13,-8765.89,-44.7367], Tmin=(1041.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.2847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C=C(C)C(C)=O(2682)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {7,D}
5  C u0 p0 c0 {1,D} {3,S} {4,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-31.3577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36527,0.0581975,-4.00173e-05,1.38698e-08,-2.00157e-12,-3677.49,21.8454], Tmin=(100,'K'), Tmax=(1543.77,'K')), NASAPolynomial(coeffs=[11.7054,0.0314059,-1.39855e-05,2.62822e-09,-1.81104e-13,-6870.05,-32.5296], Tmin=(1543.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.3577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2][C]1OC(=C)C1C(2605)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {2,S} {6,S}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u1 p0 c0 {4,S} {14,S} {15,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (287.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631489,0.0606009,-1.79311e-05,-3.38022e-08,2.2136e-11,34689.7,22.1569], Tmin=(100,'K'), Tmax=(900.422,'K')), NASAPolynomial(coeffs=[19.5348,0.0137832,-1.83815e-06,1.13105e-10,-5.17946e-15,29779.3,-75.422], Tmin=(900.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C2CsJOC(O)) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=[C]C(C)=C(C)[O](2683)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,D} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {4,D}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (142.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410109,0.0746505,-7.60167e-05,4.18578e-08,-9.18563e-12,17277.4,24.4404], Tmin=(100,'K'), Tmax=(1111.46,'K')), NASAPolynomial(coeffs=[14.2861,0.0247123,-8.621e-06,1.43287e-09,-9.28317e-14,14192.8,-43.9695], Tmin=(1111.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)=O(2151)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (263.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,284.408,2341.18],'cm^-1')),
        HinderedRotor(inertia=(0.156252,'amu*angstrom^2'), symmetry=1, barrier=(8.96009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155628,'amu*angstrom^2'), symmetry=1, barrier=(8.96524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156348,'amu*angstrom^2'), symmetry=1, barrier=(8.971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223558,'amu*angstrom^2'), symmetry=1, barrier=(12.8763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.760429,0.0777028,-0.000106618,9.04779e-08,-3.15747e-11,31800.4,27.5938], Tmin=(100,'K'), Tmax=(803.239,'K')), NASAPolynomial(coeffs=[6.72392,0.0391904,-1.82368e-05,3.46057e-09,-2.38982e-13,31126.7,1.90027], Tmin=(803.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(C)C(C)=O(2684)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,D} {2,S} {4,S}
6  C u1 p0 c0 {2,S} {7,D}
7  C u1 p0 c0 {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (300.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.248962,'amu*angstrom^2'), symmetry=1, barrier=(5.72413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248932,'amu*angstrom^2'), symmetry=1, barrier=(5.72344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248942,'amu*angstrom^2'), symmetry=1, barrier=(5.72367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248968,'amu*angstrom^2'), symmetry=1, barrier=(5.72427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50462,0.0631101,-4.74808e-05,-7.96387e-09,2.81422e-11,36171.9,25.0752], Tmin=(100,'K'), Tmax=(527.815,'K')), NASAPolynomial(coeffs=[6.43215,0.0392767,-1.81412e-05,3.47106e-09,-2.42659e-13,35463.6,2.66881], Tmin=(527.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (-2.69897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (449.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (141.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (91.7757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (358.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (433.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (554.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (5.58535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (75.5481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (228.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (123.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (157.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (106.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (115.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (127.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (91.0886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (120.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (302.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (106.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (111.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (149.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (154.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (298.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (147.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (50.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (139.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (269.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (389.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (428.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (5.58535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (75.5481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (239.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (99.5107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (127.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (114.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['CH2CO(27)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', 'C=[C]CC(=C)[O](1431)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.2329e+07,'m^3/(mol*s)'), n=-0.144618, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_Sp-7R!H=6R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2]C(C)=CC(=C)[O](2671)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2]C(=CC)C(=C)[O](2421)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=[C]C=C([O])CC(2672)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3431.13,'s^-1'), n=2.77015, Ea=(361.127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_2R!H->C"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=C(57)', '[CH2]C([O])=CC(574)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(7)', 'C=[C]C(C)[C]=C(2528)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=C1OC(=C)C1C(2424)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=C=C(C)C(=C)O(2673)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=[C]C(C)[C]1CO1(2674)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2][C]1CC(=O)C1C(2575)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2]C1([O])C(=C)C1C(2675)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(159.843,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 157.6 to 159.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(18)', 'C=C=CC(=C)[O](1143)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0143,'m^3/(mol*s)'), n=2.41, Ea=(31.2265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Ext-2R!H-R_N-4R!H-inRing_N-Sp-5R!H-2R!H_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Ext-2R!H-R_N-4R!H-inRing_N-Sp-5R!H-2R!H_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(167)', 'C=C=CC(1692)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00742,'m^3/(mol*s)'), n=2.41, Ea=(28.5908,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-2C-R_Sp-5R!H-1R!H_N-Sp-6R!H-2C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-2C-R_Sp-5R!H-1R!H_N-Sp-6R!H-2C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(6)', 'C=C=C(C)C(=C)[O](2676)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(58.9142,'m^3/(mol*s)'), n=1.72001, Ea=(13.9497,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2675934934080712, var=0.0038470988293168424, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2CO(27)', '[CH2][C]=CC(1693)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.23758e-13,'m^3/(mol*s)'), n=5.21655, Ea=(9.81748,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5310967300932249, var=1.2868917210605535, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-1R!H-R_N-7R!H-inRing_Sp-2R!H=1R!H_Sp-5R!H=4R!H_Sp-6R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-1R!H-R_N-7R!H-inRing_Sp-2R!H=1R!H_Sp-5R!H=4R!H_Sp-6R!H-3R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(6)', 'C#CC(C)C(=C)[O](2677)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(172591,'m^3/(mol*s)'), n=0.983154, Ea=(7.4342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=O(167)', '[CH2][C]=CC(1693)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=C(C)C(=C)[O](2141)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.03e+09,'s^-1'), n=1.31, Ea=(263.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 156 used for R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC(C)C(=C)[O](2146)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2]C(C=C)C(=C)[O](1859)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2][C]=C(C)C(=C)O(2678)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(O)C(C)[C]=C(2679)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([C]=C)C(=C)O(2145)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C([O])C(C)C=C(2147)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]C(C)C(=C)O(2680)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=CC(=C)OC(2681)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(389.776,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_N-2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_N-2R!H->C"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(T)(17)', 'C=[C]C(C)[C]=O(2654)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=C1CC(=O)C1C(2582)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=C=C(C)C(C)=O(2682)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['[CH2][C]1OC(=C)C1C(2605)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(C)C(=C)[O](2142)'],
    products = ['C=[C]C(C)=C(C)[O](2683)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(50.7042,'s^-1'), n=3.11103, Ea=(102.21,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([C]=C)C(C)=O(2151)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]C(C)C(C)=O(2684)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1224',
    isomers = [
        'C=[C]C(C)C(=C)[O](2142)',
    ],
    reactants = [
        ('CH2CO(27)', 'C=C=CC(1692)'),
        ('[CH2][C]=O(167)', 'C=C=CC(1692)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1224',
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

