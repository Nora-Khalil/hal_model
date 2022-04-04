species(
    label = '[CH2]C(C=C)O[C]=C(1861)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (331.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,331.347,331.347,331.347,331.347],'cm^-1')),
        HinderedRotor(inertia=(0.215428,'amu*angstrom^2'), symmetry=1, barrier=(16.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215428,'amu*angstrom^2'), symmetry=1, barrier=(16.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215428,'amu*angstrom^2'), symmetry=1, barrier=(16.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215428,'amu*angstrom^2'), symmetry=1, barrier=(16.784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286392,0.0707372,-5.02297e-05,5.80443e-09,5.2661e-12,40072.7,30.4205], Tmin=(100,'K'), Tmax=(973.687,'K')), NASAPolynomial(coeffs=[19.0646,0.0178115,-6.00259e-06,1.06628e-09,-7.57746e-14,35267.9,-65.5676], Tmin=(973.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CJO)"""),
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
    label = 'C=[C]O[CH]CC=C(2097)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (328.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239843,0.0725759,-5.62508e-05,1.24664e-08,2.83788e-12,39640.2,29.745], Tmin=(100,'K'), Tmax=(977.425,'K')), NASAPolynomial(coeffs=[18.8413,0.0182045,-6.1928e-06,1.09279e-09,-7.68202e-14,34964.8,-64.8868], Tmin=(977.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC[CH]C=C(1860)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u1 p0 c0 {2,S} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (251.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.389331,0.0671197,-3.70248e-05,-6.52675e-09,8.98921e-12,30374.5,26.9344], Tmin=(100,'K'), Tmax=(985.202,'K')), NASAPolynomial(coeffs=[18.534,0.0205569,-7.40139e-06,1.35419e-09,-9.71111e-14,25483.7,-67.0091], Tmin=(985.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=CJO)"""),
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
    label = '[CH2]C([O])C=C(1930)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (252.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,529.352,529.574],'cm^-1')),
        HinderedRotor(inertia=(0.0693079,'amu*angstrom^2'), symmetry=1, barrier=(13.7795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0692816,'amu*angstrom^2'), symmetry=1, barrier=(13.7798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95358,0.036167,-6.30155e-06,-1.92919e-08,1.03178e-11,30464.9,21.1322], Tmin=(100,'K'), Tmax=(988.725,'K')), NASAPolynomial(coeffs=[12.0826,0.0154802,-5.70142e-06,1.06015e-09,-7.65758e-14,27470.1,-32.6352], Tmin=(988.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = 'C=[C]O[CH]C=C(1324)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {3,S} {4,D} {7,S}
3  C u1 p0 c0 {1,S} {2,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (274.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,803.538,805.882,808.916],'cm^-1')),
        HinderedRotor(inertia=(0.0828863,'amu*angstrom^2'), symmetry=1, barrier=(17.715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768929,'amu*angstrom^2'), symmetry=1, barrier=(17.6792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52124,'amu*angstrom^2'), symmetry=1, barrier=(34.9762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62293,0.0439765,-1.78828e-05,-8.61649e-09,6.47691e-12,33135.9,24.3831], Tmin=(100,'K'), Tmax=(1018.23,'K')), NASAPolynomial(coeffs=[12.2576,0.0195857,-7.56396e-06,1.39647e-09,-9.86729e-14,30068.9,-31.5409], Tmin=(1018.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=CJO)"""),
)

species(
    label = 'C=CC1CC(=C)O1(1865)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {7,D} {11,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (12.0051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37464,0.0317946,7.43374e-05,-1.34033e-07,5.84597e-11,1562.56,19.8668], Tmin=(100,'K'), Tmax=(915.873,'K')), NASAPolynomial(coeffs=[22.0823,0.00885878,1.34562e-06,-4.28247e-10,2.41062e-14,-5061.73,-93.6722], Tmin=(915.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.0051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=COC(=C)C=C(2098)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,D}
3  C u0 p0 c0 {2,S} {6,D} {8,S}
4  C u0 p0 c0 {1,S} {7,D} {9,S}
5  C u0 p0 c0 {2,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {12,S} {13,S}
7  C u0 p0 c0 {4,D} {14,S} {15,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (16.4855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.680567,0.0614951,-2.92391e-05,-1.1282e-08,1.02204e-11,2112.56,24.9728], Tmin=(100,'K'), Tmax=(974.364,'K')), NASAPolynomial(coeffs=[16.9248,0.02102,-7.28041e-06,1.30249e-09,-9.2325e-14,-2297.22,-59.3591], Tmin=(974.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.4855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]OC1[CH]CC1(2099)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {4,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (352.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43,0.03307,6.22489e-05,-1.11967e-07,4.71935e-11,42463.4,27.2449], Tmin=(100,'K'), Tmax=(946.021,'K')), NASAPolynomial(coeffs=[19.5848,0.0155035,-3.75902e-06,6.93948e-10,-5.8395e-14,36379.5,-73.3346], Tmin=(946.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1[CH]CC(=C)O1(2100)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {3,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u1 p0 c0 {2,S} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (224.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37512,0.0321879,7.03458e-05,-1.2611e-07,5.40389e-11,27087.4,23.9926], Tmin=(100,'K'), Tmax=(931.784,'K')), NASAPolynomial(coeffs=[21.5239,0.0110456,-8.25616e-07,8.43043e-11,-1.51066e-14,20495.5,-87.0138], Tmin=(931.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1CC1O[C]=C(2101)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (352.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02362,0.0459093,2.70348e-05,-7.96055e-08,3.74881e-11,42541.8,27.0026], Tmin=(100,'K'), Tmax=(924.924,'K')), NASAPolynomial(coeffs=[20.1989,0.013393,-1.98541e-06,2.38011e-10,-2.03496e-14,36838.4,-75.6669], Tmin=(924.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1OC(=C)C1[CH2](2102)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u0 p0 c0 {4,D} {14,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (292.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752884,0.0536927,1.10775e-05,-7.17857e-08,3.83773e-11,35274.9,21.3153], Tmin=(100,'K'), Tmax=(882.563,'K')), NASAPolynomial(coeffs=[21.6549,0.00898703,2.03228e-06,-7.25941e-10,5.51645e-14,29637,-87.9519], Tmin=(882.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Isobutyl) + radical(CJC(C)OC)"""),
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
    label = 'C2H3(28)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=[C]OC=C(472)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {3,D} {6,S}
3  C u0 p0 c0 {2,D} {7,S} {8,S}
4  C u0 p0 c0 {5,D} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {4,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (209.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,343.919,343.919,343.919],'cm^-1')),
        HinderedRotor(inertia=(0.221761,'amu*angstrom^2'), symmetry=1, barrier=(18.6133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22176,'amu*angstrom^2'), symmetry=1, barrier=(18.6133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0818,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21803,0.0314889,-3.55293e-06,-1.99569e-08,1.05741e-11,25275.1,20.564], Tmin=(100,'K'), Tmax=(960.345,'K')), NASAPolynomial(coeffs=[10.9617,0.0135903,-4.52395e-06,7.98529e-10,-5.66263e-14,22741.7,-25.7119], Tmin=(960.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
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
    label = 'C=[C]OC(=C)C=C(2103)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (256.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,1685,370,337.472,337.472,337.472,337.472],'cm^-1')),
        HinderedRotor(inertia=(0.194117,'amu*angstrom^2'), symmetry=1, barrier=(15.6879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194117,'amu*angstrom^2'), symmetry=1, barrier=(15.6879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194117,'amu*angstrom^2'), symmetry=1, barrier=(15.6879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.917244,0.0628892,-5.4552e-05,2.49004e-08,-4.59091e-12,30932.3,27.2542], Tmin=(100,'K'), Tmax=(1296.2,'K')), NASAPolynomial(coeffs=[13.5522,0.0238984,-9.43071e-06,1.69345e-09,-1.14946e-13,27656.9,-36.9802], Tmin=(1296.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C#COC([CH2])C=C(2104)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,T}
7  C u0 p0 c0 {6,T} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (301.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0441,'amu*angstrom^2'), symmetry=1, barrier=(24.0059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04367,'amu*angstrom^2'), symmetry=1, barrier=(23.9961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04438,'amu*angstrom^2'), symmetry=1, barrier=(24.0124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04362,'amu*angstrom^2'), symmetry=1, barrier=(23.9948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0833032,0.0745418,-5.92021e-05,1.04339e-08,4.84286e-12,36422.3,25.4773], Tmin=(100,'K'), Tmax=(962.83,'K')), NASAPolynomial(coeffs=[21.4124,0.0124509,-3.78477e-06,6.69143e-10,-4.93762e-14,31085.8,-82.9989], Tmin=(962.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C=C(C)O[C]=C(2105)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (283.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14418,0.060938,-4.63582e-05,1.87425e-08,-3.16115e-12,34170,28.5085], Tmin=(100,'K'), Tmax=(1365.49,'K')), NASAPolynomial(coeffs=[11.5712,0.0303935,-1.28048e-05,2.36085e-09,-1.61921e-13,31322.4,-25.044], Tmin=(1365.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([CH2])C=C(2106)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u1 p0 c0 {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (339.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912613,'amu*angstrom^2'), symmetry=1, barrier=(20.9828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911512,'amu*angstrom^2'), symmetry=1, barrier=(20.9575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911324,'amu*angstrom^2'), symmetry=1, barrier=(20.9531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911377,'amu*angstrom^2'), symmetry=1, barrier=(20.9544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0211896,0.0722837,-3.92396e-05,-1.6584e-08,1.59871e-11,40970.6,28.6471], Tmin=(100,'K'), Tmax=(938.714,'K')), NASAPolynomial(coeffs=[23.3244,0.0110566,-2.23867e-06,3.43347e-10,-2.74619e-14,34918.2,-91.2384], Tmin=(938.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C([CH2])OC=C(2107)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {7,D} {9,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u0 p0 c0 {4,D} {14,S} {15,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (202.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666618,0.0624702,-3.25015e-05,-7.56589e-09,8.89692e-12,24474.3,27.7839], Tmin=(100,'K'), Tmax=(971.853,'K')), NASAPolynomial(coeffs=[16.5688,0.0216334,-7.46283e-06,1.31885e-09,-9.24448e-14,20221,-54.4612], Tmin=(971.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]OC(C)[C]=C(2108)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {6,D} {12,S} {13,S}
5  C u0 p0 c0 {7,D} {14,S} {15,S}
6  C u1 p0 c0 {2,S} {4,D}
7  C u1 p0 c0 {1,S} {5,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (359.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,372.179,372.179,372.179,372.179],'cm^-1')),
        HinderedRotor(inertia=(0.166519,'amu*angstrom^2'), symmetry=1, barrier=(16.368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166519,'amu*angstrom^2'), symmetry=1, barrier=(16.368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166519,'amu*angstrom^2'), symmetry=1, barrier=(16.368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166519,'amu*angstrom^2'), symmetry=1, barrier=(16.368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614926,0.0624543,-3.10904e-05,-1.0231e-08,9.92328e-12,43349,29.963], Tmin=(100,'K'), Tmax=(982.352,'K')), NASAPolynomial(coeffs=[17.6777,0.0198412,-7.04246e-06,1.28705e-09,-9.25034e-14,38700.4,-58.6486], Tmin=(982.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([C]=C)OC=C(2109)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (330.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.857121,'amu*angstrom^2'), symmetry=1, barrier=(19.7069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856772,'amu*angstrom^2'), symmetry=1, barrier=(19.6989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856521,'amu*angstrom^2'), symmetry=1, barrier=(19.6931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857006,'amu*angstrom^2'), symmetry=1, barrier=(19.7043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0610095,0.0733777,-4.74641e-05,-4.50952e-09,1.08213e-11,39854.3,28.5744], Tmin=(100,'K'), Tmax=(946.449,'K')), NASAPolynomial(coeffs=[21.89,0.0133966,-3.55406e-06,5.91931e-10,-4.37203e-14,34276.7,-83.1733], Tmin=(946.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)O[C]=C(2110)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,D} {12,S}
5  C u0 p0 c0 {6,D} {13,S} {14,S}
6  C u1 p0 c0 {1,S} {5,D}
7  C u1 p0 c0 {4,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (368.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,309.924,309.927,309.927],'cm^-1')),
        HinderedRotor(inertia=(0.256639,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25663,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256639,'amu*angstrom^2'), symmetry=1, barrier=(17.4931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256645,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580623,0.0613017,-2.26984e-05,-2.2454e-08,1.51156e-11,44465.1,30.0155], Tmin=(100,'K'), Tmax=(965.726,'K')), NASAPolynomial(coeffs=[19.0542,0.0175976,-5.78188e-06,1.05128e-09,-7.72986e-14,39366.9,-66.3867], Tmin=(965.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(C)C=C(2111)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  C u1 p0 c0 {1,S} {7,D}
7  C u1 p0 c0 {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (368.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,309.924,309.927,309.927],'cm^-1')),
        HinderedRotor(inertia=(0.256639,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25663,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256639,'amu*angstrom^2'), symmetry=1, barrier=(17.4931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256645,'amu*angstrom^2'), symmetry=1, barrier=(17.493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580623,0.0613017,-2.26984e-05,-2.2454e-08,1.51156e-11,44465.1,30.0155], Tmin=(100,'K'), Tmax=(965.726,'K')), NASAPolynomial(coeffs=[19.0542,0.0175976,-5.78188e-06,1.05128e-09,-7.72986e-14,39366.9,-66.3867], Tmin=(965.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])OC=C(2112)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {7,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  C u1 p0 c0 {3,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (339.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912613,'amu*angstrom^2'), symmetry=1, barrier=(20.9828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911512,'amu*angstrom^2'), symmetry=1, barrier=(20.9575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911324,'amu*angstrom^2'), symmetry=1, barrier=(20.9531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911377,'amu*angstrom^2'), symmetry=1, barrier=(20.9544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0211896,0.0722837,-3.92396e-05,-1.6584e-08,1.59871e-11,40970.6,28.6471], Tmin=(100,'K'), Tmax=(938.714,'K')), NASAPolynomial(coeffs=[23.3244,0.0110566,-2.23867e-06,3.43347e-10,-2.74619e-14,34918.2,-91.2384], Tmin=(938.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_P)"""),
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
    E0 = (62.4496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (221.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (254.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (583.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (386.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (70.7339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (151.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (186.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (105.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (99.6872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (150.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (81.9726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (62.4496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (243.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (241.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (252.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (184.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (189.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (175.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (204.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (241.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (161.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (143.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (132.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (102.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['CH2CO(27)', 'C=CC=C(61)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['C=[C]O[CH]CC=C(2097)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['C=[C]OC[CH]C=C(1860)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(191.946,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=C(57)', '[CH2]C([O])C=C(1930)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(17)', 'C=[C]O[CH]C=C(1324)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['C=CC1CC(=C)O1(1865)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['C=COC(=C)C=C(2098)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['C=[C]OC1[CH]CC1(2099)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['[CH2]C1[CH]CC(=C)O1(2100)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.47e+07,'s^-1'), n=0.85, Ea=(43.5136,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 7 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['[CH2]C1CC1O[C]=C(2101)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['[CH2]C1OC(=C)C1[CH2](2102)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.65009e+08,'s^-1'), n=1.00067, Ea=(87.7294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2CO(27)', '[CH2]C=C[CH2](74)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.0083,'m^3/(mol*s)'), n=2.41, Ea=(137.617,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-2R!H-R_N-7R!H-inRing_Sp-6R!H-3R_Sp-2R!H=1R!H_Sp-5R!H=4R!H_N-Sp-7R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-2R!H-R_N-7R!H-inRing_Sp-6R!H-3R_Sp-2R!H=1R!H_Sp-5R!H=4R!H_N-Sp-7R!H-2R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=O(167)', 'C=CC=C(61)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.1206e-20,'m^3/(mol*s)'), n=7.37188, Ea=(75.3495,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing
Multiplied by reaction path degeneracy 2.0
Ea raised from 73.0 to 75.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['C2H3(28)', 'C=[C]OC=C(472)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(161.758,'m^3/(mol*s)'), n=1.13767, Ea=(16.7537,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.011839139166081196, var=0.8671401158115553, Tref=1000.0, N=160, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(6)', 'C=[C]OC(=C)C=C(2103)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.68706e+22,'m^3/(mol*s)'), n=-4.56653, Ea=(42.6271,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.23436359975141804, var=0.013647226894766704, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(6)', 'C#COC([CH2])C=C(2104)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(125.384,'m^3/(mol*s)'), n=2.03877, Ea=(8.36527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]=O(167)', '[CH2]C=C[CH2](74)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(19.4451,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['[CH2]C=C(C)O[C]=C(2105)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.08414e+09,'s^-1'), n=1.12201, Ea=(127.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_Cd] + [R2H_S;C_rad_out_2H;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=COC([CH2])C=C(2106)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C)O[C]=C(1861)'],
    products = ['[CH2]C=C([CH2])OC=C(2107)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]OC(C)[C]=C(2108)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([C]=C)OC=C(2109)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.8284271247461903
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC(C)O[C]=C(2110)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]OC(C)C=C(2111)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=CC([CH2])OC=C(2112)'],
    products = ['[CH2]C(C=C)O[C]=C(1861)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 3.605551275463989
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1097',
    isomers = [
        '[CH2]C(C=C)O[C]=C(1861)',
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
    label = 'PDepNetwork #1097',
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

