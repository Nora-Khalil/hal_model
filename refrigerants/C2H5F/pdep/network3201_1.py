species(
    label = '[CH2]C(=CF)C(=C)[O](8478)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {3,D} {8,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-60.0492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,279.899,280.012],'cm^-1')),
        HinderedRotor(inertia=(0.400065,'amu*angstrom^2'), symmetry=1, barrier=(22.2709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400502,'amu*angstrom^2'), symmetry=1, barrier=(22.2717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.373245,0.0688993,-7.27295e-05,3.78398e-08,-7.51254e-12,-7082.4,23.4934], Tmin=(100,'K'), Tmax=(1336.43,'K')), NASAPolynomial(coeffs=[18.1436,0.0111242,-2.73409e-06,3.54577e-10,-1.98552e-14,-11422.5,-65.859], Tmin=(1336.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.0492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
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
    label = 'C=C=CF(5941)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C([O])C[C]=CF(8480)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {7,D} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (58.0469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1685,370,306.527,306.538,306.544,4000],'cm^-1')),
        HinderedRotor(inertia=(0.246596,'amu*angstrom^2'), symmetry=1, barrier=(16.445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246715,'amu*angstrom^2'), symmetry=1, barrier=(16.4504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3729.08,'J/mol'), sigma=(6.01568,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.47 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16902,0.0585594,-5.7575e-05,2.95577e-08,-6.03825e-12,7086.5,25.8717], Tmin=(100,'K'), Tmax=(1187.85,'K')), NASAPolynomial(coeffs=[13.0177,0.0186596,-7.1901e-06,1.27969e-09,-8.67376e-14,4271.6,-33.3312], Tmin=(1187.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.0469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
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
    label = 'C=C([O])[C]=CF(2679)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {6,D} {9,S}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (25.0377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47248,'amu*angstrom^2'), symmetry=1, barrier=(33.8553,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54819,0.0530411,-6.23292e-05,3.36366e-08,-5.52875e-12,3100.63,17.9969], Tmin=(100,'K'), Tmax=(817.484,'K')), NASAPolynomial(coeffs=[12.3345,0.0104489,-2.86653e-06,3.85885e-10,-2.12992e-14,996.757,-33.9487], Tmin=(817.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.0377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCFH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'O(6)',
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
    label = '[CH2]C([C]=C)=CF(8828)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u1 p0 c0 {2,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,D} {7,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (215.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.80824,'amu*angstrom^2'), symmetry=1, barrier=(41.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80277,'amu*angstrom^2'), symmetry=1, barrier=(41.4492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01568,0.056136,-5.40796e-05,2.66274e-08,-5.07585e-12,25978.2,20.0421], Tmin=(100,'K'), Tmax=(1383.51,'K')), NASAPolynomial(coeffs=[14.6826,0.0131418,-3.69184e-06,5.28901e-10,-3.12927e-14,22529.6,-49.1257], Tmin=(1383.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1OCC1=CF(8490)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {1,S} {4,D} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-148.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49798,0.0413216,5.79463e-06,-4.24187e-08,2.05216e-11,-17795.9,17.1383], Tmin=(100,'K'), Tmax=(962.635,'K')), NASAPolynomial(coeffs=[15.8993,0.0148073,-4.82087e-06,8.97281e-10,-6.77704e-14,-22112.7,-59.8119], Tmin=(962.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'C=C([O])[C]1CC1F(8975)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {6,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (21.4886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86207,0.0334944,2.03472e-05,-5.21223e-08,2.26726e-11,2673.61,23.0287], Tmin=(100,'K'), Tmax=(970.239,'K')), NASAPolynomial(coeffs=[13.6946,0.0176951,-6.21845e-06,1.16858e-09,-8.67922e-14,-1174.91,-41.6994], Tmin=(970.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(C=C(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C([CH]F)=C1CO1(8976)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {7,S}
6  C u1 p0 c0 {1,S} {5,S} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (70.9726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1928,0.0460984,1.93674e-06,-4.61909e-08,2.38488e-11,8651.41,21.1919], Tmin=(100,'K'), Tmax=(944.981,'K')), NASAPolynomial(coeffs=[19.0286,0.00962748,-2.11936e-06,3.7352e-10,-3.18721e-14,3538.01,-73.066], Tmin=(944.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.9726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(Allyl_P) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[O]C1=C([CH]F)CC1(8977)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u1 p0 c0 {1,S} {5,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1.07215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62575,0.0382142,1.23727e-05,-4.83028e-08,2.24633e-11,-30.789,22.9376], Tmin=(100,'K'), Tmax=(959.564,'K')), NASAPolynomial(coeffs=[15.4764,0.0148169,-4.733e-06,8.76631e-10,-6.63717e-14,-4269.85,-51.5499], Tmin=(959.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.07215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[CH2][C]1C(=C)OC1F(8978)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (44.9575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66242,0.0292105,5.89398e-05,-1.14587e-07,5.1853e-11,5512.18,18.8411], Tmin=(100,'K'), Tmax=(900.226,'K')), NASAPolynomial(coeffs=[20.8048,0.00320036,3.89433e-06,-9.63585e-10,6.5201e-14,-326.862,-84.7863], Tmin=(900.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.9575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([O])CC1=CF(8944)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (227.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40735,0.0492872,-2.87905e-05,1.43019e-09,3.13123e-12,27446.8,26.1114], Tmin=(100,'K'), Tmax=(1039.67,'K')), NASAPolynomial(coeffs=[12.7257,0.0202595,-7.85704e-06,1.43863e-09,-1.00579e-13,24308.7,-32.7071], Tmin=(1039.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cd(Cd-F)-Cs-Cs) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1([CH]F)OC1=C(8932)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {3,S} {8,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (150.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.418035,0.0784632,-8.85855e-05,4.62648e-08,-8.88877e-12,18224.3,24.2676], Tmin=(100,'K'), Tmax=(1479.8,'K')), NASAPolynomial(coeffs=[22.7251,0.00261533,1.76959e-06,-5.10501e-10,3.88913e-14,12830,-91.5377], Tmin=(1479.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CJC(C)OC) + radical(CsCsF1sH)"""),
)

species(
    label = '[CH2][C]=O(140)',
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
    label = 'C=[C][CH]F(5942)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (196.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.900551,'amu*angstrom^2'), symmetry=1, barrier=(20.7054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79302e-05,1.96926e-08,-6.19852e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.897,'K')), NASAPolynomial(coeffs=[4.81186,0.0161179,-7.47821e-06,1.46101e-09,-1.03887e-13,23340.4,3.53844], Tmin=(747.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C(C)=[C]F(8979)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  C u1 p0 c0 {1,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (45.0278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.567985,'amu*angstrom^2'), symmetry=1, barrier=(13.0591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567551,'amu*angstrom^2'), symmetry=1, barrier=(13.0491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.761715,0.0699029,-8.357e-05,5.15997e-08,-1.24844e-11,5533.38,22.472], Tmin=(100,'K'), Tmax=(1016.5,'K')), NASAPolynomial(coeffs=[14.0078,0.0177771,-6.64825e-06,1.14947e-09,-7.62387e-14,2840.54,-41.6489], Tmin=(1016.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.0278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CdCdF1s)"""),
)

species(
    label = '[CH]=C(O)C([CH2])=CF(8980)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {11,S}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {3,D} {8,S}
7  C u1 p0 c0 {4,D} {12,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (49.2422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.10412,'amu*angstrom^2'), symmetry=1, barrier=(25.386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10234,'amu*angstrom^2'), symmetry=1, barrier=(25.345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10488,'amu*angstrom^2'), symmetry=1, barrier=(25.4034,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.235267,0.078398,-8.89506e-05,4.74394e-08,-9.42612e-12,6087.61,24.5783], Tmin=(100,'K'), Tmax=(1392.49,'K')), NASAPolynomial(coeffs=[21.8904,0.00522829,2.22931e-07,-2.10532e-10,1.86944e-14,857.602,-86.1446], Tmin=(1392.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.2422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C(=C)[C]F(8981)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {12,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,D} {8,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u2 p0 c0 {1,S} {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (37.5112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,228.234,228.305,228.34,228.394],'cm^-1')),
        HinderedRotor(inertia=(0.517862,'amu*angstrom^2'), symmetry=1, barrier=(19.1196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.517591,'amu*angstrom^2'), symmetry=1, barrier=(19.1198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.517054,'amu*angstrom^2'), symmetry=1, barrier=(19.1197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353875,0.0735887,-8.41589e-05,4.67961e-08,-9.99259e-12,4648.51,22.1908], Tmin=(100,'K'), Tmax=(1159.25,'K')), NASAPolynomial(coeffs=[18.3189,0.0115997,-3.94829e-06,6.6786e-10,-4.46521e-14,483.35,-67.1348], Tmin=(1159.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]=C([O])C(C)=CF(8982)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (35.5478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,194,682,905,1196,1383,3221,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.66546,'amu*angstrom^2'), symmetry=1, barrier=(15.3002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66662,'amu*angstrom^2'), symmetry=1, barrier=(15.3269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658234,0.0707274,-8.34944e-05,5.01763e-08,-1.17399e-11,4398.16,22.5916], Tmin=(100,'K'), Tmax=(1053.69,'K')), NASAPolynomial(coeffs=[15.0518,0.0160864,-5.70859e-06,9.61086e-10,-6.2929e-14,1364.92,-47.6017], Tmin=(1053.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.5478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C(=C)OF(8983)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,D} {8,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u2 p0 c0 {3,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (370.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,275.929,275.93,275.931,275.931,275.931],'cm^-1')),
        HinderedRotor(inertia=(0.929991,'amu*angstrom^2'), symmetry=1, barrier=(50.2465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92999,'amu*angstrom^2'), symmetry=1, barrier=(50.2465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.929988,'amu*angstrom^2'), symmetry=1, barrier=(50.2465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827105,0.0694221,-6.94312e-05,3.7428e-08,-8.22675e-12,44727,22.8065], Tmin=(100,'K'), Tmax=(1091.53,'K')), NASAPolynomial(coeffs=[12.2945,0.0273991,-1.16828e-05,2.15766e-09,-1.48622e-13,42223.5,-33.5217], Tmin=(1091.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(CF)C(=C)[O](8984)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (51.8383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.780383,'amu*angstrom^2'), symmetry=1, barrier=(17.9425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779628,'amu*angstrom^2'), symmetry=1, barrier=(17.9252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751513,0.0694296,-8.00066e-05,4.69509e-08,-1.07943e-11,6353.47,23.0528], Tmin=(100,'K'), Tmax=(1067.5,'K')), NASAPolynomial(coeffs=[14.6618,0.0173065,-6.76571e-06,1.21102e-09,-8.23897e-14,3383.61,-44.9651], Tmin=(1067.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.8383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(F)C(=C)[O](8479)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (67.0442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.883489,'amu*angstrom^2'), symmetry=1, barrier=(20.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.887086,'amu*angstrom^2'), symmetry=1, barrier=(20.3959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3671.58,'J/mol'), sigma=(5.97962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.49 K, Pc=38.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23828,0.0653743,-8.35037e-05,6.15886e-08,-1.88195e-11,8158.78,23.2569], Tmin=(100,'K'), Tmax=(791.209,'K')), NASAPolynomial(coeffs=[8.53884,0.0284626,-1.35191e-05,2.61484e-09,-1.83751e-13,7003.63,-10.2537], Tmin=(791.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.0442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]F(837)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'C=[C]C(=C)[O](715)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {5,S}
3 C u0 p0 c0 {2,D} {6,S} {7,S}
4 C u0 p0 c0 {5,D} {8,S} {9,S}
5 C u1 p0 c0 {2,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (216.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44423,'amu*angstrom^2'), symmetry=1, barrier=(33.2056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65189,0.0455512,-4.93437e-05,2.74298e-08,-5.8039e-12,26168.2,16.7593], Tmin=(100,'K'), Tmax=(1302.66,'K')), NASAPolynomial(coeffs=[11.8538,0.00923595,-1.78225e-06,1.49146e-10,-4.09098e-15,23933.6,-33.5317], Tmin=(1302.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1OC(F)C1=C(8489)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-189.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33466,0.0431498,7.17503e-06,-4.67764e-08,2.24006e-11,-22672.8,14.6027], Tmin=(100,'K'), Tmax=(970.987,'K')), NASAPolynomial(coeffs=[17.6048,0.0134532,-4.61516e-06,9.11269e-10,-7.14538e-14,-27592.1,-72.4737], Tmin=(970.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = '[CH2]C1=C([O])CC1F(8985)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-28.2183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48993,0.0427811,-1.00112e-06,-3.37273e-08,1.7062e-11,-3292.33,21.2754], Tmin=(100,'K'), Tmax=(968.832,'K')), NASAPolynomial(coeffs=[15.1927,0.0159343,-5.46131e-06,1.01288e-09,-7.49031e-14,-7342.63,-51.5988], Tmin=(968.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1=C([CH]F)CO1(8986)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {4,D} {7,S}
6  C u1 p0 c0 {1,S} {4,S} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (53.1587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3197,0.0430518,9.02379e-06,-5.20107e-08,2.55286e-11,6504.54,21.9166], Tmin=(100,'K'), Tmax=(946.238,'K')), NASAPolynomial(coeffs=[18.3528,0.01071,-2.57999e-06,4.61337e-10,-3.81288e-14,1505.47,-68.6995], Tmin=(946.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.1587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1([O])C(=C)C1F(8961)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (221.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34634,0.0536593,-4.31144e-05,1.76626e-08,-2.92488e-12,26790.5,24.4579], Tmin=(100,'K'), Tmax=(1424.85,'K')), NASAPolynomial(coeffs=[12.9419,0.0211071,-8.84533e-06,1.62863e-09,-1.11617e-13,23486.1,-35.5896], Tmin=(1424.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cd(Cd)-Cs-Cs(F)) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH]C(=CF)C(=C)O(8987)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {11,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {1,S} {3,D} {8,S}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u2 p0 c0 {3,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (21.3313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08112,'amu*angstrom^2'), symmetry=1, barrier=(47.8491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08131,'amu*angstrom^2'), symmetry=1, barrier=(47.8535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08113,'amu*angstrom^2'), symmetry=1, barrier=(47.8494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.113317,0.0766696,-7.79953e-05,3.92298e-08,-7.53956e-12,2725.42,24.5696], Tmin=(100,'K'), Tmax=(1395.93,'K')), NASAPolynomial(coeffs=[19.7637,0.0135123,-3.46707e-06,4.54892e-10,-2.54297e-14,-2219.89,-75.7922], Tmin=(1395.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.3313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C(=C)CF(8988)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (51.8383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.780383,'amu*angstrom^2'), symmetry=1, barrier=(17.9425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779628,'amu*angstrom^2'), symmetry=1, barrier=(17.9252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751513,0.0694296,-8.00066e-05,4.69509e-08,-1.07943e-11,6353.47,23.0528], Tmin=(100,'K'), Tmax=(1067.5,'K')), NASAPolynomial(coeffs=[14.6618,0.0173065,-6.76571e-06,1.21102e-09,-8.23897e-14,3383.61,-44.9651], Tmin=(1067.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.8383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=CF)C[C]=O(8484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {4,D} {10,S}
7  C u1 p0 c0 {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-12.5249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,1855,455,950,572.999],'cm^-1')),
        HinderedRotor(inertia=(0.0966575,'amu*angstrom^2'), symmetry=1, barrier=(22.8676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0130383,'amu*angstrom^2'), symmetry=1, barrier=(3.14094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995359,'amu*angstrom^2'), symmetry=1, barrier=(22.8853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3528.4,'J/mol'), sigma=(5.75757,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.13 K, Pc=41.95 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07395,0.054536,-4.20173e-05,1.51593e-08,-2.13739e-12,-1392.76,24.9094], Tmin=(100,'K'), Tmax=(1692.01,'K')), NASAPolynomial(coeffs=[17.8618,0.0148494,-6.83495e-06,1.29742e-09,-8.92847e-14,-7073.91,-64.9116], Tmin=(1692.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.5249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=O)=CF(8989)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 C u0 p0 c0 {1,S} {3,D} {7,S}
6 C u1 p0 c0 {2,D} {3,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (3.48227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0728568,'amu*angstrom^2'), symmetry=1, barrier=(15.0125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.944075,'amu*angstrom^2'), symmetry=1, barrier=(21.7061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11906,0.0445391,-4.47775e-05,2.33036e-08,-5.02496e-12,483.921,17.4653], Tmin=(100,'K'), Tmax=(1090.57,'K')), NASAPolynomial(coeffs=[9.01903,0.0192312,-9.96819e-06,2.02453e-09,-1.46952e-13,-1021.05,-16.4212], Tmin=(1090.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.48227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(CdCFH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C=O)CJ) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'O=C1CCC1=CF(8783)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {1,S} {5,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-203.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17538,0.0318727,7.17643e-06,-2.39553e-08,8.85109e-12,-24457.4,20.1938], Tmin=(100,'K'), Tmax=(1145.43,'K')), NASAPolynomial(coeffs=[9.12645,0.0272297,-1.24511e-05,2.43071e-09,-1.7353e-13,-27337.6,-19.9064], Tmin=(1145.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C1([CH]F)CC1=O(8849)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  C u1 p0 c0 {1,S} {3,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (152.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712977,0.0730308,-9.06467e-05,5.85114e-08,-1.48972e-11,18496.3,22.4031], Tmin=(100,'K'), Tmax=(963.434,'K')), NASAPolynomial(coeffs=[13.4598,0.0201078,-8.2481e-06,1.4933e-09,-1.01489e-13,16040.2,-38.6179], Tmin=(963.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(CsCsFHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CJC(C)2C=O) + radical(CsCsF1sH)"""),
)

species(
    label = 'C=C([C]F)C(C)=O(8990)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u2 p0 c0 {1,S} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (25.2365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,180,180,180,1156.91],'cm^-1')),
        HinderedRotor(inertia=(0.00341688,'amu*angstrom^2'), symmetry=1, barrier=(3.2468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141302,'amu*angstrom^2'), symmetry=1, barrier=(3.2488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141386,'amu*angstrom^2'), symmetry=1, barrier=(3.25074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68378,0.0568063,-7.52803e-05,6.9965e-08,-2.81456e-11,3113.02,22.3393], Tmin=(100,'K'), Tmax=(723.307,'K')), NASAPolynomial(coeffs=[3.80978,0.0374427,-1.93492e-05,3.87457e-09,-2.77131e-13,3004.44,14.1466], Tmin=(723.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.2365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]C(=C)C(=O)CF(8991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u2 p0 c0 {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (36.9806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,350,440,435,1725,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,460.883,461.145,461.224,461.449,461.884],'cm^-1')),
        HinderedRotor(inertia=(0.341759,'amu*angstrom^2'), symmetry=1, barrier=(51.6737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343095,'amu*angstrom^2'), symmetry=1, barrier=(51.6841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34119,'amu*angstrom^2'), symmetry=1, barrier=(51.6784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19861,0.0478468,-3.27564e-06,-7.55579e-08,7.40588e-11,4504.13,21.3608], Tmin=(100,'K'), Tmax=(447.137,'K')), NASAPolynomial(coeffs=[3.98849,0.0424852,-2.10172e-05,4.16354e-09,-2.98089e-13,4237.6,12.9759], Tmin=(447.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.9806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1C(=O)CC1F(8796)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-212.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14757,0.0357598,-6.36014e-06,-8.70778e-09,3.34053e-12,-25462.3,17.7538], Tmin=(100,'K'), Tmax=(1365.92,'K')), NASAPolynomial(coeffs=[10.2159,0.0270063,-1.30815e-05,2.54496e-09,-1.7816e-13,-29054,-28.7665], Tmin=(1365.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH]C(=CF)C(C)=O(8992)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u2 p0 c0 {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (9.05652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,375,552.5,462.5,1710,194,682,905,1196,1383,3221,499.674,499.689,499.713,499.72],'cm^-1')),
        HinderedRotor(inertia=(0.296173,'amu*angstrom^2'), symmetry=1, barrier=(52.479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296162,'amu*angstrom^2'), symmetry=1, barrier=(52.4784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000675109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93958,0.0510165,-3.60943e-05,1.57412e-08,-3.52558e-12,1158.63,22.1443], Tmin=(100,'K'), Tmax=(910.639,'K')), NASAPolynomial(coeffs=[4.17224,0.0412095,-1.99403e-05,3.91511e-09,-2.78934e-13,751.998,11.582], Tmin=(910.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.05652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCFH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-119.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (160.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (338.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (390.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-119.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-34.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (58.6943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (2.60321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-1.60597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (159.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (83.1103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (119.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (67.5129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (288.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (133.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (173.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (14.0275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (12.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (351.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (159.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (93.7266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (363.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-119.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-1.20248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-1.60597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (154.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (287.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-2.1525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (28.3545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (165.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (317.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-119.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (86.1946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1.75273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (151.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-119.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-14.4272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['CH2CO(28)', 'C=C=CF(5941)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(8.69051,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 8.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C[C]=CF(8480)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(18)', 'C=C([O])[C]=CF(2679)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', '[CH2]C([C]=C)=CF(8828)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['C=C1OCC1=CF(8490)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['C=C([O])[C]1CC1F(8975)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C([CH]F)=C1CO1(8976)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[O]C1=C([CH]F)CC1(8977)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2][C]1C(=C)OC1F(8978)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1([O])CC1=CF(8944)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(287.423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 285.1 to 287.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1([CH]F)OC1=C(8932)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(210.952,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(140)', 'C=C=CF(5941)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.0142,'m^3/(mol*s)'), n=2.41, Ea=(17.6024,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_N-Sp-5R!H-1R!H_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_N-Sp-5R!H-1R!H_Ext-3R-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CO(28)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(140)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C(C)=[C]F(8979)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(O)C([CH2])=CF(8980)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(O)C(=C)[C]F(8981)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_single;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C([O])C(C)=CF(8982)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C(=C)C(=C)OF(8983)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(48.3926,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(CF)C(=C)[O](8984)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(175.382,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(F)C(=C)[O](8479)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]F(837)', 'C=[C]C(=C)[O](715)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['C=C1OC(F)C1=C(8489)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1=C([O])CC1F(8985)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1=C([CH]F)CO1(8986)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1([O])C(=C)C1F(8961)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(281.973,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 280.3 to 282.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CHF(40)', 'C=[C]C(=C)[O](715)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=CF)C(=C)O(8987)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C([O])C(=C)CF(8988)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=CF)C[C]=O(8484)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(T)(18)', '[CH2]C([C]=O)=CF(8989)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['O=C1CCC1=CF(8783)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['[CH2]C1([CH]F)CC1=O(8849)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(214.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(CO)_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([C]F)C(C)=O(8990)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C(=C)C(=O)CF(8991)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(182.724,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=CF)C(=C)[O](8478)'],
    products = ['C=C1C(=O)CC1F(8796)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_1H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C(=CF)C(C)=O(8992)'],
    products = ['[CH2]C(=CF)C(=C)[O](8478)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3201',
    isomers = [
        '[CH2]C(=CF)C(=C)[O](8478)',
    ],
    reactants = [
        ('CH2CO(28)', 'C=C=CF(5941)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3201',
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

