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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27872,0.0628602,-5.68104e-05,3.03764e-08,-7.03817e-12,29315.5,27.4237], Tmin=(100,'K'), Tmax=(1004.69,'K')), NASAPolynomial(coeffs=[8.24352,0.0351311,-1.54108e-05,2.90554e-09,-2.02504e-13,27916,-6.21003], Tmin=(1004.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = 'C=[C]CC[C]=O(1433)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {1,D} {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (274.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.227495,'amu*angstrom^2'), symmetry=1, barrier=(5.23055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227471,'amu*angstrom^2'), symmetry=1, barrier=(5.23001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227544,'amu*angstrom^2'), symmetry=1, barrier=(5.23169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3431.29,'J/mol'), sigma=(5.7779,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.96 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05041,0.0483781,-4.02502e-05,6.28079e-09,1.08062e-11,33103.9,22.7018], Tmin=(100,'K'), Tmax=(558.731,'K')), NASAPolynomial(coeffs=[5.91183,0.0295609,-1.343e-05,2.55485e-09,-1.78231e-13,32534.6,5.08733], Tmin=(558.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=CC)C[C]=O(2423)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {13,S}
6  C u1 p0 c0 {4,S} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (143.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939432,0.0577971,-3.2494e-05,4.50347e-09,1.09047e-12,17345.8,26.3324], Tmin=(100,'K'), Tmax=(1266.22,'K')), NASAPolynomial(coeffs=[14.9879,0.0270627,-1.22491e-05,2.35475e-09,-1.65548e-13,12694.3,-49.0784], Tmin=(1266.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C)=CC[C]=O(2563)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {4,D} {13,S}
6  C u1 p0 c0 {4,S} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (143.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939432,0.0577971,-3.2494e-05,4.50347e-09,1.09047e-12,17345.8,26.3324], Tmin=(100,'K'), Tmax=(1266.22,'K')), NASAPolynomial(coeffs=[14.9879,0.0270627,-1.22491e-05,2.35475e-09,-1.65548e-13,12694.3,-49.0784], Tmin=(1266.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

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
    label = 'C[CH]C[C]=O(560)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {3,S} {11,S}
5  C u1 p0 c0 {1,D} {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (132.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,204.848],'cm^-1')),
        HinderedRotor(inertia=(1.26364e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63642,'amu*angstrom^2'), symmetry=1, barrier=(18.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229607,'amu*angstrom^2'), symmetry=1, barrier=(6.83811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35404,0.0356784,-2.29653e-05,7.91175e-09,-1.1666e-12,16026.6,20.8557], Tmin=(100,'K'), Tmax=(1491.53,'K')), NASAPolynomial(coeffs=[7.59581,0.0216208,-8.82768e-06,1.59263e-09,-1.07417e-13,14463,-6.52845], Tmin=(1491.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(C)[C]=C(378)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {10,S} {11,S}
4  C u0 p0 c0 {5,D} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (394.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.0389829,'amu*angstrom^2'), symmetry=1, barrier=(13.0455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574506,'amu*angstrom^2'), symmetry=1, barrier=(13.209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125946,'amu*angstrom^2'), symmetry=1, barrier=(79.8742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94624,0.0392083,-1.10823e-05,-1.04211e-08,6.34926e-12,47544.7,21.7037], Tmin=(100,'K'), Tmax=(985.328,'K')), NASAPolynomial(coeffs=[8.76023,0.0246347,-8.82085e-06,1.52962e-09,-1.03282e-13,45566.5,-14.2931], Tmin=(985.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C(=O)CC1C(2564)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,S} {7,D}
6  C u0 p0 c0 {1,D} {3,S} {5,S}
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
    E0 = (-43.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91898,0.0330705,2.9767e-05,-5.27137e-08,1.94795e-11,-5195.87,20.8668], Tmin=(100,'K'), Tmax=(1056.51,'K')), NASAPolynomial(coeffs=[10.0401,0.0332178,-1.43049e-05,2.77392e-09,-1.99815e-13,-8636.1,-26.9192], Tmin=(1056.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C=C(C)CC=O(2565)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u0 p0 c0 {1,D} {2,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (8.38493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09217,0.0562613,-3.15525e-05,6.16665e-09,-3.71262e-15,1119.46,23.7552], Tmin=(100,'K'), Tmax=(1488.98,'K')), NASAPolynomial(coeffs=[15.9522,0.0268509,-1.25117e-05,2.38172e-09,-1.64118e-13,-4470.79,-57.7633], Tmin=(1488.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.38493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(C)C=C=O(2059)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u0 p0 c0 {1,D} {5,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-31.8108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9074,0.0635453,-5.0945e-05,2.16043e-08,-3.75027e-12,-3711.03,24.6205], Tmin=(100,'K'), Tmax=(1355.76,'K')), NASAPolynomial(coeffs=[13.117,0.0275222,-1.10891e-05,2.00589e-09,-1.36338e-13,-7021.68,-37.9999], Tmin=(1355.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.8108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d)"""),
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
    label = 'C=C=CC[C]=O(1297)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u0 p0 c0 {5,D} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {4,D}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (207.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.835049,'amu*angstrom^2'), symmetry=1, barrier=(19.1994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836186,'amu*angstrom^2'), symmetry=1, barrier=(19.2256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87225,0.0394272,-1.81596e-05,-1.99219e-09,2.4259e-12,25027.1,20.9957], Tmin=(100,'K'), Tmax=(1224.98,'K')), NASAPolynomial(coeffs=[12.1376,0.0189863,-9.14502e-06,1.81781e-09,-1.30449e-13,21530.9,-34.6172], Tmin=(1224.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
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
    label = 'C=C=C(C)C[C]=O(2566)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {6,D} {13,S} {14,S}
6  C u0 p0 c0 {4,D} {5,D}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (168.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,217.885],'cm^-1')),
        HinderedRotor(inertia=(0.459175,'amu*angstrom^2'), symmetry=1, barrier=(15.4292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457736,'amu*angstrom^2'), symmetry=1, barrier=(15.4279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458248,'amu*angstrom^2'), symmetry=1, barrier=(15.4288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25896,0.0567586,-3.94324e-05,1.29234e-08,-1.69122e-12,20348.3,23.9487], Tmin=(100,'K'), Tmax=(1751.06,'K')), NASAPolynomial(coeffs=[15.8184,0.0235002,-1.09427e-05,2.07681e-09,-1.42666e-13,15249.4,-54.4489], Tmin=(1751.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
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
    label = 'C=[C]C(C)C=C=O(2567)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,D} {12,S}
5  C u0 p0 c0 {6,D} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {5,D}
7  C u0 p0 c0 {1,D} {4,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (206.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.570782,'amu*angstrom^2'), symmetry=1, barrier=(13.1234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570657,'amu*angstrom^2'), symmetry=1, barrier=(13.1205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570398,'amu*angstrom^2'), symmetry=1, barrier=(13.1146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14668,0.0648709,-6.41109e-05,3.57785e-08,-8.34073e-12,24880.8,24.2417], Tmin=(100,'K'), Tmax=(1018.93,'K')), NASAPolynomial(coeffs=[9.9501,0.0303115,-1.3235e-05,2.49149e-09,-1.7359e-13,23086.8,-18.3948], Tmin=(1018.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C[C]=O(2568)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {7,T}
6  C u1 p0 c0 {1,D} {3,S}
7  C u0 p0 c0 {5,T} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (170.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1855,455,950,750,770,3400,2100,314.112],'cm^-1')),
        HinderedRotor(inertia=(0.133999,'amu*angstrom^2'), symmetry=1, barrier=(9.36699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133704,'amu*angstrom^2'), symmetry=1, barrier=(9.36003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134136,'amu*angstrom^2'), symmetry=1, barrier=(9.36179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09486,'amu*angstrom^2'), symmetry=1, barrier=(76.9158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09392,0.0643991,-6.60314e-05,3.92747e-08,-9.64674e-12,20557.4,25.4152], Tmin=(100,'K'), Tmax=(979.881,'K')), NASAPolynomial(coeffs=[9.86495,0.0285948,-1.12226e-05,1.98546e-09,-1.33059e-13,18838.5,-16.7217], Tmin=(979.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C](C)C[C]=O(2066)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (137.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,2622.57],'cm^-1')),
        HinderedRotor(inertia=(0.526983,'amu*angstrom^2'), symmetry=1, barrier=(12.2948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0290077,'amu*angstrom^2'), symmetry=1, barrier=(12.3963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985294,'amu*angstrom^2'), symmetry=1, barrier=(2.26538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.70685,'amu*angstrom^2'), symmetry=1, barrier=(80.4307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31246,0.0562455,-3.875e-05,1.43709e-08,-2.24777e-12,16586.7,25.0716], Tmin=(100,'K'), Tmax=(1449.26,'K')), NASAPolynomial(coeffs=[10.5775,0.0306735,-1.22826e-05,2.19567e-09,-1.47507e-13,13901.2,-23.0646], Tmin=(1449.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C(C)C=C[O](2569)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {1,S} {4,D} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (225.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5951,0.0663934,-4.98611e-05,1.50179e-08,-3.14007e-13,27274.4,26.3748], Tmin=(100,'K'), Tmax=(1048.84,'K')), NASAPolynomial(coeffs=[15.3199,0.0240421,-9.03613e-06,1.61831e-09,-1.11426e-13,23426.2,-48.986], Tmin=(1048.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C[C]=O(2072)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {7,D} {14,S}
6  C u1 p0 c0 {1,D} {3,S}
7  C u1 p0 c0 {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (252.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,339.321],'cm^-1')),
        HinderedRotor(inertia=(0.112887,'amu*angstrom^2'), symmetry=1, barrier=(9.25966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111767,'amu*angstrom^2'), symmetry=1, barrier=(9.25869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112338,'amu*angstrom^2'), symmetry=1, barrier=(9.2612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113826,'amu*angstrom^2'), symmetry=1, barrier=(9.26362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14129,0.0630321,-5.34796e-05,2.51758e-08,-4.98379e-12,30435.8,27.838], Tmin=(100,'K'), Tmax=(1179.35,'K')), NASAPolynomial(coeffs=[10.4522,0.0314522,-1.33134e-05,2.47046e-09,-1.70664e-13,28239.6,-18.6177], Tmin=(1179.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)[CH][C]=O(2068)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u1 p0 c0 {2,S} {7,S} {13,S}
6  C u0 p0 c0 {4,D} {14,S} {15,S}
7  C u1 p0 c0 {1,D} {5,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (172.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0479,0.0570796,-3.3289e-05,3.75541e-09,2.48021e-12,20876.4,27.5295], Tmin=(100,'K'), Tmax=(1044.09,'K')), NASAPolynomial(coeffs=[12.5452,0.0270082,-1.01651e-05,1.81076e-09,-1.23879e-13,17713.8,-32.0825], Tmin=(1044.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = 'C=[C][C](C)CC=O(2570)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {3,S} {7,S}
5  C u0 p0 c0 {1,D} {2,S} {13,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (214.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68908,0.0539868,-3.52325e-05,1.26455e-08,-2.03158e-12,25935,23.3703], Tmin=(100,'K'), Tmax=(1320.7,'K')), NASAPolynomial(coeffs=[7.23922,0.0371771,-1.61406e-05,3.00822e-09,-2.07299e-13,24469,-4.94978], Tmin=(1320.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=C)CC=O(2071)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  C u0 p0 c0 {1,D} {3,S} {11,S}
6  C u0 p0 c0 {7,D} {14,S} {15,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (288.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,2082.25,2083.72],'cm^-1')),
        HinderedRotor(inertia=(0.28708,'amu*angstrom^2'), symmetry=1, barrier=(11.8562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285789,'amu*angstrom^2'), symmetry=1, barrier=(11.8655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289659,'amu*angstrom^2'), symmetry=1, barrier=(11.8648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68791,'amu*angstrom^2'), symmetry=1, barrier=(69.0832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28387,0.0625788,-5.77002e-05,3.29767e-08,-8.24906e-12,34742.3,28.3687], Tmin=(100,'K'), Tmax=(934.842,'K')), NASAPolynomial(coeffs=[7.43233,0.0362701,-1.54857e-05,2.87148e-09,-1.97992e-13,33592.7,-0.879639], Tmin=(934.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(C)CC=O(2571)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,D} {3,S} {14,S}
6  C u1 p0 c0 {2,S} {7,D}
7  C u1 p0 c0 {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (330.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,271.362],'cm^-1')),
        HinderedRotor(inertia=(0.00228846,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191864,'amu*angstrom^2'), symmetry=1, barrier=(10.0282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191962,'amu*angstrom^2'), symmetry=1, barrier=(10.0284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191964,'amu*angstrom^2'), symmetry=1, barrier=(10.029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33588,0.0626844,-5.56464e-05,2.96212e-08,-6.97104e-12,39792.5,26.8065], Tmin=(100,'K'), Tmax=(977.978,'K')), NASAPolynomial(coeffs=[7.50474,0.0374535,-1.69483e-05,3.24173e-09,-2.27741e-13,38585.9,-2.8174], Tmin=(977.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (22.4975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (473.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (116.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (167.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (268.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (512.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (613.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (30.7818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (100.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (111.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (85.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (116.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (152.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (160.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (119.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (198.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (165.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (300.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (205.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (181.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (137.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (164.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (174.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (127.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (142.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (142.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['CH2CO(27)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', 'C=[C]CC[C]=O(1433)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(973284,'m^3/(mol*s)'), n=0.225671, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['[CH2]C(=CC)C[C]=O(2423)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['[CH2]C(C)=CC[C]=O(2563)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=[C]C(C)C(=C)[O](2142)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=C(57)', 'C[CH]C[C]=O(560)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=O(152)', '[CH2]C(C)[C]=C(378)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=C1C(=O)CC1C(2564)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=C=C(C)CC=O(2565)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=CC(C)C=C=O(2059)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO(14)', '[CH2]C(C)[C]=C(378)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(30.522,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(167)', 'C=C=CC(1692)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(30.7822,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(18)', 'C=C=CC[C]=O(1297)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0290254,'m^3/(mol*s)'), n=2.199, Ea=(28.943,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.12990640417617208, var=2.9840165643953642, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Ext-2R!H-R_N-4R!H-inRing_N-Sp-5R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Ext-2R!H-R_N-4R!H-inRing_N-Sp-5R!H-2R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(6)', 'C=C=C(C)C[C]=O(2566)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0.854359,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2CO(27)', '[CH2][C]=CC(1693)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00415,'m^3/(mol*s)'), n=2.41, Ea=(39.3938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-2R!H-R_N-7R!H-inRing_Sp-6R!H-3R_Sp-2R!H=1R!H_Sp-5R!H=4R!H_N-Sp-7R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-2R!H-R_N-7R!H-inRing_Sp-6R!H-3R_Sp-2R!H=1R!H_Sp-5R!H=4R!H_N-Sp-7R!H-2R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(6)', 'C=[C]C(C)C=C=O(2567)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(0.985,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(6)', 'C#CC(C)C[C]=O(2568)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(172591,'m^3/(mol*s)'), n=0.983154, Ea=(4.09865,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=O(167)', '[CH2][C]=CC(1693)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=C[C](C)C[C]=O(2066)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=[C]C(C)C=C[O](2569)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC(C)C[C]=O(2072)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=CC(C)[CH][C]=O(2068)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['[CH2]C(C=C)C[C]=O(1863)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C[C]=O(2069)'],
    products = ['C=[C][C](C)CC=O(2570)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.98524e+06,'s^-1'), n=1.75122, Ea=(105.056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([C]=C)CC=O(2071)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(405.74,'s^-1'), n=2.51765, Ea=(74.9333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]C(C)CC=O(2571)'],
    products = ['C=[C]C(C)C[C]=O(2069)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1230',
    isomers = [
        'C=[C]C(C)C[C]=O(2069)',
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
    label = 'PDepNetwork #1230',
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

