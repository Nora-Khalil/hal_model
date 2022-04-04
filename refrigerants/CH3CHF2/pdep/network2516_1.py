species(
    label = '[CH]=C(C)C[C]=O(6347)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {1,D} {2,S}
6  C u1 p0 c0 {4,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (274.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.118816,'amu*angstrom^2'), symmetry=1, barrier=(12.9087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118729,'amu*angstrom^2'), symmetry=1, barrier=(12.9082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118729,'amu*angstrom^2'), symmetry=1, barrier=(12.9096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62794,0.0468405,-3.212e-05,1.02103e-08,-1.27794e-12,33148.1,22.1808], Tmin=(100,'K'), Tmax=(1855.64,'K')), NASAPolynomial(coeffs=[15.6016,0.0167189,-7.77124e-06,1.46262e-09,-9.94022e-14,27962.1,-53.8727], Tmin=(1855.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = 'C#CC(5603)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (172.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.64],'cm^-1')),
        HinderedRotor(inertia=(0.630825,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30522,0.0109266,1.31984e-05,-2.25157e-08,8.87389e-12,20738.2,7.19167], Tmin=(100,'K'), Tmax=(969.992,'K')), NASAPolynomial(coeffs=[5.80059,0.0113535,-4.03487e-06,7.19054e-10,-5.02177e-14,19749.9,-7.36988], Tmin=(969.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CC[C]=O(542)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3 C u0 p0 c0 {2,S} {5,D} {8,S}
4 C u1 p0 c0 {1,D} {2,S}
5 C u1 p0 c0 {3,D} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (313.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.182146,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182136,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3285.42,'J/mol'), sigma=(5.46087,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.18 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43871,0.0274412,-4.68074e-06,-1.14208e-08,5.23481e-12,37817.6,18.4996], Tmin=(100,'K'), Tmax=(1137.06,'K')), NASAPolynomial(coeffs=[10.2746,0.0146221,-7.22317e-06,1.47543e-09,-1.08342e-13,35082.3,-24.5022], Tmin=(1137.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C)C(=C)[O](6343)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {6,D}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  C u1 p0 c0 {3,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (227.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06833,0.0594635,-5.64811e-05,2.44128e-08,-2.97418e-12,27452.6,20.265], Tmin=(100,'K'), Tmax=(950.332,'K')), NASAPolynomial(coeffs=[14.1614,0.0156261,-5.08015e-06,8.35639e-10,-5.51698e-14,24455.1,-44.9125], Tmin=(950.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[C]=O(128)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.49],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.0892,0.0020037,-1.61643e-05,2.55031e-08,-1.16411e-11,52802.7,4.52493], Tmin=(100,'K'), Tmax=(856.126,'K')), NASAPolynomial(coeffs=[0.961549,0.00569059,-3.48052e-06,7.19222e-10,-5.08058e-14,53738.7,21.4667], Tmin=(856.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(6403)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  C u2 p0 c0 {2,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (336.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,376.905,377.215,377.226],'cm^-1')),
        HinderedRotor(inertia=(0.510427,'amu*angstrom^2'), symmetry=1, barrier=(51.472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.510854,'amu*angstrom^2'), symmetry=1, barrier=(51.4778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52506,0.0264751,7.08405e-06,-2.07737e-08,8.00454e-12,40556.1,14.6774], Tmin=(100,'K'), Tmax=(1061.32,'K')), NASAPolynomial(coeffs=[6.48667,0.0249945,-9.83289e-06,1.79345e-09,-1.24009e-13,38957.7,-8.23953], Tmin=(1061.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1=CC(=O)C1(6471)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {1,D} {2,S} {6,S}
6  C u0 p0 c0 {4,D} {5,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-0.0136073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45851,0.0226891,3.2671e-05,-4.90191e-08,1.68379e-11,63.1711,15.5769], Tmin=(100,'K'), Tmax=(1114.14,'K')), NASAPolynomial(coeffs=[9.60024,0.027029,-1.35353e-05,2.78137e-09,-2.04963e-13,-3388.95,-28.0003], Tmin=(1114.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.0136073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C(C)C=C=O(6472)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u0 p0 c0 {1,D} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-13.6411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68156,0.0546109,-6.10694e-05,4.20194e-08,-1.23728e-11,-1560.41,19.1792], Tmin=(100,'K'), Tmax=(810.322,'K')), NASAPolynomial(coeffs=[6.9412,0.0286473,-1.30069e-05,2.4769e-09,-1.72983e-13,-2412.8,-5.08921], Tmin=(810.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.6411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'CO(13)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=[C]C(5604)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (490.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859241,'amu*angstrom^2'), symmetry=1, barrier=(1.97556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04067e-06,1.20986e-09,-3.01512e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.69,'K')), NASAPolynomial(coeffs=[6.96312,0.0100051,-3.7071e-06,6.32809e-10,-4.05654e-14,57370.1,-10.4646], Tmin=(1807.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C(C)C=C=O(6473)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u1 p0 c0 {3,D} {11,S}
6  C u0 p0 c0 {1,D} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (233.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.646468,'amu*angstrom^2'), symmetry=1, barrier=(14.8636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649413,'amu*angstrom^2'), symmetry=1, barrier=(14.9313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47371,0.0598601,-8.44902e-05,6.93062e-08,-2.30691e-11,28165,20.3189], Tmin=(100,'K'), Tmax=(808.325,'K')), NASAPolynomial(coeffs=[7.5737,0.0251157,-1.15564e-05,2.17734e-09,-1.49603e-13,27327.8,-6.89081], Tmin=(808.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=O(136)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = 'C#CC[C]=O(553)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3 C u0 p0 c0 {2,S} {5,T}
4 C u1 p0 c0 {1,D} {2,S}
5 C u0 p0 c0 {3,T} {8,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (231.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,1855,455,950,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.41231,'amu*angstrom^2'), symmetry=1, barrier=(32.4717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41208,'amu*angstrom^2'), symmetry=1, barrier=(32.4666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0659,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86091,0.0115829,0.000101893,-3.28908e-07,3.12745e-10,27895.3,10.1603], Tmin=(10,'K'), Tmax=(356.529,'K')), NASAPolynomial(coeffs=[3.91283,0.025347,-1.63757e-05,5.10676e-09,-6.10596e-13,27800.4,8.68463], Tmin=(356.529,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(231.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""C#CC[C]DO""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = '[CH]C(C)=CC=O(6474)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {1,D} {4,S} {11,S}
6  C u2 p0 c0 {3,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (214.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23707,0.0466622,-2.57112e-05,6.15675e-09,-5.72063e-13,25824.6,18.3924], Tmin=(100,'K'), Tmax=(2346.45,'K')), NASAPolynomial(coeffs=[15.5117,0.024033,-1.12453e-05,2.04678e-09,-1.34173e-13,19594.9,-56.9721], Tmin=(2346.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)=C[C]=O(6475)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {6,S} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u1 p0 c0 {1,D} {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (155.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54945,0.0590457,-6.39397e-05,3.95624e-08,-1.04434e-11,18810.4,19.6195], Tmin=(100,'K'), Tmax=(894.505,'K')), NASAPolynomial(coeffs=[8.09987,0.0297541,-1.48207e-05,2.95469e-09,-2.12236e-13,17638.5,-11.2523], Tmin=(894.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C(=C)C[C]=O(6476)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {9,S} {12,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (179.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66358,0.0415021,-1.03515e-05,-1.46005e-08,7.52907e-12,21653.2,21.2213], Tmin=(100,'K'), Tmax=(1094.22,'K')), NASAPolynomial(coeffs=[13.0595,0.0203006,-9.33095e-06,1.86332e-09,-1.36064e-13,17934.6,-40.3798], Tmin=(1094.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C(=C)CC=O(6477)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {6,S}
4  C u0 p0 c0 {1,D} {2,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u2 p0 c0 {3,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (238.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72854,0.040921,-2.71141e-06,-1.80493e-08,7.28651e-12,28773,21.6991], Tmin=(100,'K'), Tmax=(1175.75,'K')), NASAPolynomial(coeffs=[10.8668,0.0297602,-1.38967e-05,2.70866e-09,-1.9247e-13,25246.6,-29.7251], Tmin=(1175.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (44.9826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (503.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (290.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (545.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (53.2669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (108.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (64.7426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (207.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (221.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (134.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (171.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (421.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (162.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (237.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (174.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (146.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['CH2CO(28)', 'C#CC(5603)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[CH]=CC[C]=O(542)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.16873e+07,'m^3/(mol*s)'), n=-0.471498, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07514139785501385, var=0.03091755307632206, Tref=1000.0, N=6, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['[CH]=C(C)C(=C)[O](6343)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=O(128)', '[CH]C(=C)C(6403)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['CC1=CC(=O)C1(6471)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['C=C(C)C=C=O(6472)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO(13)', '[CH]C(=C)C(6403)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(76.6447,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CO(28)', '[CH]=[C]C(5604)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0154,'m^3/(mol*s)'), n=2.41, Ea=(6.94651,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH]=C(C)C=C=O(6473)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1306.26,'m^3/(mol*s)'), n=1.34286, Ea=(6.52113,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2935929756704599, var=0.2435564102543282, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=O(136)', 'C#CC(5603)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(32.2309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH3(19)', 'C#CC[C]=O(553)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00276756,'m^3/(mol*s)'), n=2.82919, Ea=(33.4671,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.07604778385317046, var=1.4951309799213957, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(136)', '[CH]=[C]C(5604)'],
    products = ['[CH]=C(C)C[C]=O(6347)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.88995e+09,'m^3/(mol*s)'), n=-0.866642, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=4.73600691749923e-16, var=0.18719384195212546, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_3BrCO->O_Ext-2CF-R_N-5R!H->O_Ext-2CF-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_3BrCO->O_Ext-2CF-R_N-5R!H->O_Ext-2CF-R
Ea raised from -13.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['[CH]C(C)=CC=O(6474)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.27137e+08,'s^-1'), n=1.53496, Ea=(117.681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['[CH2]C(C)=C[C]=O(6475)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['[CH2]C(=C)C[C]=O(6476)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(C)C[C]=O(6347)'],
    products = ['[CH]C(=C)CC=O(6477)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52029e+07,'s^-1'), n=1.31708, Ea=(101.589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2516',
    isomers = [
        '[CH]=C(C)C[C]=O(6347)',
    ],
    reactants = [
        ('CH2CO(28)', 'C#CC(5603)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2516',
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

