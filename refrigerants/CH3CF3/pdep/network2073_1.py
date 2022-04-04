species(
    label = '[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {13,S}
11 C u1 p0 c0 {8,S} {10,S} {14,S}
12 C u0 p0 c0 {7,D} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1037.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,1380,1390,370,380,2900,435,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.396716,'amu*angstrom^2'), symmetry=1, barrier=(9.12129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396465,'amu*angstrom^2'), symmetry=1, barrier=(9.1155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05051,'amu*angstrom^2'), symmetry=1, barrier=(24.1532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36482,'amu*angstrom^2'), symmetry=1, barrier=(31.3799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33582,0.128518,-0.000213316,1.85001e-07,-6.3268e-11,-124576,39.6398], Tmin=(100,'K'), Tmax=(794.543,'K')), NASAPolynomial(coeffs=[14.1069,0.0386198,-2.06516e-05,4.09151e-09,-2.87319e-13,-126646,-28.8962], Tmin=(794.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1037.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'O=CC(=O)F(4234)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=CC(F)F(344)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {3,S} {4,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-792.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.142102,'amu*angstrom^2'), symmetry=1, barrier=(3.2672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2858.21,'J/mol'), sigma=(4.57471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=446.45 K, Pc=67.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74532,0.0233589,0.000133412,-5.9431e-07,7.11464e-10,-95378,11.8568], Tmin=(10,'K'), Tmax=(302.26,'K')), NASAPolynomial(coeffs=[4.97982,0.0307765,-2.12832e-05,6.89182e-09,-8.42907e-13,-95561.1,5.58312], Tmin=(302.26,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-792.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DCC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C(F)(F)[CH]C(F)F(6773)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {12,S}
10 C u1 p0 c0 {7,S} {9,S} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-930.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.229665,'amu*angstrom^2'), symmetry=1, barrier=(5.28046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304717,'amu*angstrom^2'), symmetry=1, barrier=(7.00604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232156,'amu*angstrom^2'), symmetry=1, barrier=(5.33773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3542.99,'J/mol'), sigma=(5.98741,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.41 K, Pc=37.45 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0428076,0.0955271,-0.000140277,1.07017e-07,-3.14437e-11,-111761,32.5531], Tmin=(100,'K'), Tmax=(657.887,'K')), NASAPolynomial(coeffs=[12.0353,0.0319542,-1.66289e-05,3.30323e-09,-2.34179e-13,-113542,-21.8185], Tmin=(657.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-930.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C(F)(C=O)C([C](F)F)C(F)F(6758)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {14,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 C u0 p0 c0 {7,D} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1040.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,215.766,215.767,215.768,215.768],'cm^-1')),
        HinderedRotor(inertia=(0.364851,'amu*angstrom^2'), symmetry=1, barrier=(12.0536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99582,'amu*angstrom^2'), symmetry=1, barrier=(32.8986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681578,'amu*angstrom^2'), symmetry=1, barrier=(22.5174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364849,'amu*angstrom^2'), symmetry=1, barrier=(12.0535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4111.74,'J/mol'), sigma=(6.43753,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=642.24 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31926,0.127491,-0.000208995,1.7888e-07,-6.04712e-11,-125014,39.2676], Tmin=(100,'K'), Tmax=(791.651,'K')), NASAPolynomial(coeffs=[14.559,0.0374532,-1.98088e-05,3.9108e-09,-2.74252e-13,-127221,-31.685], Tmin=(791.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1040.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(=COF)C(F)(F)[CH]C(F)F(7470)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {12,S}
7  O u1 p2 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {7,S} {8,S} {12,D}
12 C u0 p0 c0 {6,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-790.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.256677,'amu*angstrom^2'), symmetry=1, barrier=(5.90151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256623,'amu*angstrom^2'), symmetry=1, barrier=(5.90027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256832,'amu*angstrom^2'), symmetry=1, barrier=(5.90508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834733,'amu*angstrom^2'), symmetry=1, barrier=(19.1921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62287,0.135775,-0.000231974,2.01697e-07,-6.82245e-11,-94886.6,39.4349], Tmin=(100,'K'), Tmax=(820.973,'K')), NASAPolynomial(coeffs=[15.2572,0.0364875,-1.94264e-05,3.81251e-09,-2.65115e-13,-97083.9,-35.174], Tmin=(820.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-790.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'O=C(F)[CH]OC(F)(F)[CH]C(F)F(6759)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u1 p0 c0 {6,S} {12,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1098.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3000,3050,390,425,1340,1360,335,370,611,648,830,1210,1753,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3665.26,'J/mol'), sigma=(6.00103,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.50 K, Pc=38.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31452,0.129622,-0.000200363,1.51022e-07,-4.13613e-11,-131915,38.5594], Tmin=(100,'K'), Tmax=(632.916,'K')), NASAPolynomial(coeffs=[16.1812,0.0362794,-1.99736e-05,4.02367e-09,-2.86276e-13,-134475,-40.5707], Tmin=(632.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1098.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]C(F)F-2(957)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-68.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,200.056,1350.78,1616.29],'cm^-1')),
        HinderedRotor(inertia=(0.00421206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96601,0.0195049,-1.18279e-05,2.09006e-09,3.19831e-13,-8190.14,12.8831], Tmin=(100,'K'), Tmax=(1211.04,'K')), NASAPolynomial(coeffs=[7.87163,0.00800485,-3.40879e-06,6.6204e-10,-4.73288e-14,-9723.19,-13.1468], Tmin=(1211.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C(F)(C=O)[C](F)F(4474)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-550.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16106,'amu*angstrom^2'), symmetry=1, barrier=(26.695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16613,'amu*angstrom^2'), symmetry=1, barrier=(26.8117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84721,0.0763401,-0.000131581,1.15087e-07,-3.91064e-11,-66131.8,25.4983], Tmin=(100,'K'), Tmax=(821.483,'K')), NASAPolynomial(coeffs=[10.2224,0.0205072,-1.10382e-05,2.17219e-09,-1.5119e-13,-67328.5,-15.7965], Tmin=(821.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-550.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCsF1sF1s)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C=C(F)C(F)(F)[CH]C(F)F(7471)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-916.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.871991,'amu*angstrom^2'), symmetry=1, barrier=(20.0488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293183,'amu*angstrom^2'), symmetry=1, barrier=(6.74085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294002,'amu*angstrom^2'), symmetry=1, barrier=(6.75968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.470561,0.107317,-0.000158625,1.24551e-07,-3.94437e-11,-110026,34.4443], Tmin=(100,'K'), Tmax=(770.353,'K')), NASAPolynomial(coeffs=[13.3166,0.0357201,-1.92002e-05,3.87854e-09,-2.78121e-13,-112150,-28.472], Tmin=(770.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-916.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S)"""),
)

species(
    label = 'O=CC1(F)OC(C(F)F)C1(F)F(6766)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {6,S} {9,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
12 C u0 p0 c0 {7,D} {10,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1333.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0173122,0.0936756,-0.000102524,5.82536e-08,-1.35569e-11,-160302,31.8681], Tmin=(100,'K'), Tmax=(1024.8,'K')), NASAPolynomial(coeffs=[14.6636,0.0365072,-1.88446e-05,3.81646e-09,-2.76694e-13,-163304,-39.1501], Tmin=(1024.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1333.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCCFO) + group(CsCsFFH) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(O2s-Cs-Cs-Cs(F))"""),
)

species(
    label = 'O=CC(O)(F)C(F)(F)C=C(F)F(6772)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {15,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
10 C u0 p0 c0 {9,S} {12,D} {13,S}
11 C u0 p0 c0 {7,D} {8,S} {14,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1326.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00276,0.119401,-0.000179067,1.40545e-07,-4.41378e-11,-159389,35.312], Tmin=(100,'K'), Tmax=(778.943,'K')), NASAPolynomial(coeffs=[15.1409,0.0364959,-1.94102e-05,3.89411e-09,-2.77833e-13,-161904,-38.5382], Tmin=(778.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1326.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF)"""),
)

species(
    label = 'FC(F)[CH]C(F)(F)C1(F)[CH]OO1(7472)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {13,S}
11 C u1 p0 c0 {9,S} {10,S} {14,S}
12 C u1 p0 c0 {7,S} {8,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-780.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.688467,0.112589,-0.000162725,1.27102e-07,-4.03099e-11,-93684.7,36.7176], Tmin=(100,'K'), Tmax=(767.542,'K')), NASAPolynomial(coeffs=[13.192,0.040248,-2.13427e-05,4.29387e-09,-3.07467e-13,-95815.3,-26.5745], Tmin=(767.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-780.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCsJOO)"""),
)

species(
    label = '[O]C1(F)[CH]OC(C(F)F)C1(F)F(7412)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {7,S} {9,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
12 C u1 p0 c0 {6,S} {10,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1067.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15161,0.11058,-0.000138022,8.47822e-08,-1.95975e-11,-128248,30.9524], Tmin=(100,'K'), Tmax=(898.98,'K')), NASAPolynomial(coeffs=[21.3637,0.0212577,-7.10077e-06,1.13011e-09,-7.10714e-14,-132735,-77.7136], Tmin=(898.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1067.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Tetrahydrofuran) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOCs) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C1OC1(F)C(F)(F)[CH]C(F)F(7473)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
10 C u0 p0 c0 {6,S} {7,S} {8,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {12,S} {14,S}
12 C u1 p0 c0 {9,S} {11,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-975.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0999088,0.103964,-0.00012955,6.51616e-08,3.5187e-12,-117170,34.6282], Tmin=(100,'K'), Tmax=(542.203,'K')), NASAPolynomial(coeffs=[10.4887,0.0459961,-2.49199e-05,5.04589e-09,-3.6211e-13,-118614,-12.7052], Tmin=(542.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-975.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(CsCsFFH) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C1C(C(F)F)C(F)(F)C1([O])F(7474)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {11,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {7,S} {8,S} {11,S} {14,S}
11 C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
12 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-970.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.344122,0.0959571,-0.000102044,5.42535e-08,-1.15059e-11,-116536,33.5655], Tmin=(100,'K'), Tmax=(1137,'K')), NASAPolynomial(coeffs=[18.5404,0.02952,-1.43957e-05,2.86137e-09,-2.05839e-13,-120830,-59.9661], Tmin=(1137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-970.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFF) + group(CsCsFFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(O2sj(Cs-F1sCsCs)) + radical(CC(C)OJ) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[O]C(F)(C=O)C(F)=CC(F)F(7475)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {13,S}
11 C u0 p0 c0 {6,D} {7,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-887.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.558191,'amu*angstrom^2'), symmetry=1, barrier=(12.8339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558758,'amu*angstrom^2'), symmetry=1, barrier=(12.847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45922,'amu*angstrom^2'), symmetry=1, barrier=(33.5503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.530907,0.109525,-0.000178214,1.55807e-07,-5.39664e-11,-106641,34.3785], Tmin=(100,'K'), Tmax=(797.467,'K')), NASAPolynomial(coeffs=[11.3839,0.0374703,-1.95621e-05,3.84995e-09,-2.69712e-13,-108150,-17.9564], Tmin=(797.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-887.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(398)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.101,381.695],'cm^-1')),
        HinderedRotor(inertia=(0.482775,'amu*angstrom^2'), symmetry=1, barrier=(49.7784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80088,0.0290439,-3.02832e-05,1.66615e-08,-3.81631e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.22,'K')), NASAPolynomial(coeffs=[6.95396,0.0128879,-6.71487e-06,1.38086e-09,-1.01081e-13,-26608.7,-6.32755], Tmin=(1028.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
)

species(
    label = 'O=CC(=O)C(F)(F)[CH]C(F)F(7476)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,D} {7,S} {11,S}
11 C u0 p0 c0 {6,D} {10,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-946.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,217.967,240.02,240.139],'cm^-1')),
        HinderedRotor(inertia=(0.132661,'amu*angstrom^2'), symmetry=1, barrier=(5.98871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300849,'amu*angstrom^2'), symmetry=1, barrier=(10.5532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00476683,'amu*angstrom^2'), symmetry=1, barrier=(0.169386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06153,'amu*angstrom^2'), symmetry=1, barrier=(35.4712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.306132,0.103465,-0.000151857,1.20263e-07,-3.86214e-11,-113730,34.7332], Tmin=(100,'K'), Tmax=(758.548,'K')), NASAPolynomial(coeffs=[12.3911,0.0365142,-1.94733e-05,3.92211e-09,-2.80816e-13,-115657,-23.0157], Tmin=(758.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-946.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O)"""),
)

species(
    label = 'F[C](F)[CH]C(F)F(1426)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u1 p0 c0 {5,S} {7,S} {9,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-554.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,190,488,555,1236,1407,180,2250.74],'cm^-1')),
        HinderedRotor(inertia=(0.397589,'amu*angstrom^2'), symmetry=1, barrier=(9.14135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399336,'amu*angstrom^2'), symmetry=1, barrier=(9.18152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4542,0.0670926,-0.000128678,1.28357e-07,-4.81479e-11,-66597.3,23.7219], Tmin=(100,'K'), Tmax=(836.732,'K')), NASAPolynomial(coeffs=[4.34924,0.0286339,-1.55993e-05,3.09731e-09,-2.16322e-13,-66219.9,15.4209], Tmin=(836.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C(F)C(F)(F)[CH]C(F)F(6995)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1092.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,486,617,768,1157,1926,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.140717,'amu*angstrom^2'), symmetry=1, barrier=(3.23535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141398,'amu*angstrom^2'), symmetry=1, barrier=(3.25102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.630071,'amu*angstrom^2'), symmetry=1, barrier=(14.4866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374945,0.0868619,-0.000124867,9.51828e-08,-2.92827e-11,-131313,30.7811], Tmin=(100,'K'), Tmax=(791.989,'K')), NASAPolynomial(coeffs=[11.8067,0.0291239,-1.55118e-05,3.12962e-09,-2.24593e-13,-133124,-21.7043], Tmin=(791.989,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1092.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(COCsFO) + radical(CCJCC=O)"""),
)

species(
    label = '[O]C(F)(C=O)C(F)(F)C=CF(7477)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u0 p0 c0 {6,D} {8,S} {13,S}
11 C u0 p0 c0 {4,S} {9,D} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-881.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.660758,'amu*angstrom^2'), symmetry=1, barrier=(15.1921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.659404,'amu*angstrom^2'), symmetry=1, barrier=(15.161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25788,'amu*angstrom^2'), symmetry=1, barrier=(28.9212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.359882,0.106106,-0.000154777,1.13925e-07,-3.06665e-11,-105920,32.8003], Tmin=(100,'K'), Tmax=(633.185,'K')), NASAPolynomial(coeffs=[12.9634,0.0352625,-1.85128e-05,3.68564e-09,-2.61431e-13,-107875,-27.4971], Tmin=(633.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-881.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ)"""),
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
    label = '[O]C(F)(C=O)C(F)(F)C=C(F)F(7478)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
9  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u0 p0 c0 {7,D} {9,S} {14,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1082.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.595698,'amu*angstrom^2'), symmetry=1, barrier=(13.6963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595672,'amu*angstrom^2'), symmetry=1, barrier=(13.6957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44502,'amu*angstrom^2'), symmetry=1, barrier=(33.2239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.9825,0.119345,-0.000194288,1.6577e-07,-5.6165e-11,-130076,36.0028], Tmin=(100,'K'), Tmax=(772.63,'K')), NASAPolynomial(coeffs=[13.9675,0.0355154,-1.90518e-05,3.79221e-09,-2.67571e-13,-132194,-31.0234], Tmin=(772.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1082.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ)"""),
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
    label = '[O]C=C([O])C(F)=CC(F)F(7479)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {9,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {8,D} {12,S}
8  C u0 p0 c0 {3,S} {7,D} {9,S}
9  C u0 p0 c0 {4,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-674.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,2995,3025,975,1000,1300,1375,400,500,1630,1680,280,518,736,852,873,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.853712,'amu*angstrom^2'), symmetry=1, barrier=(19.6285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852507,'amu*angstrom^2'), symmetry=1, barrier=(19.6008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.371098,0.092358,-0.000110216,6.33668e-08,-1.40713e-11,-80972.4,30.8393], Tmin=(100,'K'), Tmax=(1109.61,'K')), NASAPolynomial(coeffs=[20.8173,0.0159771,-6.96313e-06,1.33159e-09,-9.45996e-14,-85674.6,-73.5865], Tmin=(1109.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-674.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
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
    label = 'O=[C]C(=O)C(F)(F)[CH]C(F)F(7480)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,D} {7,S} {11,S}
11 C u1 p0 c0 {6,D} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-786.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,375,552.5,462.5,1710,1855,455,950,198.26,199.35,200.475],'cm^-1')),
        HinderedRotor(inertia=(0.280368,'amu*angstrom^2'), symmetry=1, barrier=(7.90038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00421403,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282716,'amu*angstrom^2'), symmetry=1, barrier=(7.88335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880989,'amu*angstrom^2'), symmetry=1, barrier=(24.8622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.283037,0.104756,-0.000157541,1.14693e-07,-2.90976e-11,-94494,35.4978], Tmin=(100,'K'), Tmax=(622.399,'K')), NASAPolynomial(coeffs=[13.5717,0.0310173,-1.67056e-05,3.34099e-09,-2.36897e-13,-96515,-27.1546], Tmin=(622.399,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-786.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)(C=O)C(F)(F)C[C](F)F(7481)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {8,S} {11,S} {13,S} {14,S}
11 C u1 p0 c0 {4,S} {5,S} {10,S}
12 C u0 p0 c0 {7,D} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1044.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53036,0.135635,-0.000237306,2.12689e-07,-7.38054e-11,-125421,38.6012], Tmin=(100,'K'), Tmax=(826.061,'K')), NASAPolynomial(coeffs=[13.2536,0.0405099,-2.18327e-05,4.29984e-09,-2.99322e-13,-127060,-25.0371], Tmin=(826.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1044.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'O=[C]C(O)(F)C(F)(F)[CH]C(F)F(7482)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {15,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {13,S}
11 C u1 p0 c0 {8,S} {10,S} {14,S}
12 C u1 p0 c0 {7,D} {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1121.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61925,0.134012,-0.000222935,1.88535e-07,-6.25307e-11,-134640,40.6974], Tmin=(100,'K'), Tmax=(797.015,'K')), NASAPolynomial(coeffs=[16.6365,0.0337642,-1.80288e-05,3.55863e-09,-2.49003e-13,-137276,-41.5152], Tmin=(797.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1121.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)C(F)(F)CC(F)F(7483)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {9,S} {15,S}
12 C u1 p0 c0 {7,D} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1084.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66108,0.135869,-0.000229859,1.97275e-07,-6.5972e-11,-130294,38.53], Tmin=(100,'K'), Tmax=(817.945,'K')), NASAPolynomial(coeffs=[16.0625,0.0348991,-1.84777e-05,3.62065e-09,-2.5167e-13,-132715,-40.4908], Tmin=(817.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1084.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC(O)(F)C(F)(F)[CH][C](F)F(7484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {15,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
10 C u1 p0 c0 {9,S} {12,S} {13,S}
11 C u0 p0 c0 {7,D} {8,S} {14,S}
12 C u1 p0 c0 {4,S} {5,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1080.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50038,0.133933,-0.000230992,2.04848e-07,-7.07995e-11,-129766,40.8103], Tmin=(100,'K'), Tmax=(813.738,'K')), NASAPolynomial(coeffs=[13.8588,0.0393186,-2.135e-05,4.22954e-09,-2.95953e-13,-131633,-26.2348], Tmin=(813.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1080.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'O=CC(F)(OF)[C](F)[CH]C(F)F(7485)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {11,S} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {14,S}
12 C u0 p0 c0 {7,D} {8,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-731.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3177,0.13004,-0.000221692,1.99884e-07,-7.10363e-11,-87806.7,41.2994], Tmin=(100,'K'), Tmax=(790.286,'K')), NASAPolynomial(coeffs=[12.214,0.0437806,-2.42412e-05,4.8675e-09,-3.44346e-13,-89390.6,-17.2879], Tmin=(790.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-731.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C(F)(C=O)[C](F)C(F)C(F)F(7486)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {14,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 C u0 p0 c0 {7,D} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1005.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,2782.5,750,1395,475,1775,1000,186.747,186.762,186.763,186.77],'cm^-1')),
        HinderedRotor(inertia=(0.489951,'amu*angstrom^2'), symmetry=1, barrier=(12.1306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489831,'amu*angstrom^2'), symmetry=1, barrier=(12.1305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36393,'amu*angstrom^2'), symmetry=1, barrier=(33.7583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36367,'amu*angstrom^2'), symmetry=1, barrier=(33.7596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39206,0.130493,-0.000219737,1.92349e-07,-6.595e-11,-120739,39.6349], Tmin=(100,'K'), Tmax=(810.604,'K')), NASAPolynomial(coeffs=[13.8225,0.0392106,-2.08343e-05,4.10284e-09,-2.86562e-13,-122673,-27.288], Tmin=(810.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1005.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C=C(OF)C(F)(F)[CH]C(F)F(7487)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u1 p2 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {6,S} {8,S} {12,D}
12 C u0 p0 c0 {7,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-786.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,225.506,225.508,225.508,225.509],'cm^-1')),
        HinderedRotor(inertia=(0.361238,'amu*angstrom^2'), symmetry=1, barrier=(13.0361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15067,'amu*angstrom^2'), symmetry=1, barrier=(5.43759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150667,'amu*angstrom^2'), symmetry=1, barrier=(5.43756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825414,'amu*angstrom^2'), symmetry=1, barrier=(29.7879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49536,0.131792,-0.000217595,1.8651e-07,-6.32891e-11,-94450.1,38.888], Tmin=(100,'K'), Tmax=(774.578,'K')), NASAPolynomial(coeffs=[15.2434,0.0376845,-2.05031e-05,4.09732e-09,-2.89479e-13,-96813.1,-36.1062], Tmin=(774.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-786.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S)"""),
)

species(
    label = '[O]C=C([O])C(F)(F)C(F)C(F)F(7488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {11,S}
7  O u1 p2 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
11 C u0 p0 c0 {6,S} {9,S} {12,D}
12 C u0 p0 c0 {7,S} {11,D} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1173.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39185,0.123452,-0.000171213,1.16976e-07,-3.12148e-11,-140998,35.8441], Tmin=(100,'K'), Tmax=(921.652,'K')), NASAPolynomial(coeffs=[20.9997,0.0262724,-1.30536e-05,2.57317e-09,-1.83087e-13,-145126,-70.3552], Tmin=(921.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1173.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=CC(F)(OF)C(F)(F)[CH][CH]F(7489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
10 C u1 p0 c0 {9,S} {12,S} {13,S}
11 C u0 p0 c0 {7,D} {8,S} {14,S}
12 C u1 p0 c0 {4,S} {10,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-738.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,334,575,1197,1424,3202,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74476,0.140785,-0.000247724,2.22412e-07,-7.74818e-11,-88582.8,40.1643], Tmin=(100,'K'), Tmax=(817.123,'K')), NASAPolynomial(coeffs=[13.931,0.0411329,-2.27238e-05,4.51952e-09,-3.16509e-13,-90379.6,-27.615], Tmin=(817.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsCsHH) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O]C(F)(C=O)C(F)(F)C(F)[CH]F(7490)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {8,S} {11,S} {13,S}
11 C u1 p0 c0 {5,S} {10,S} {14,S}
12 C u0 p0 c0 {7,D} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1008.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.871475,'amu*angstrom^2'), symmetry=1, barrier=(20.0369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873775,'amu*angstrom^2'), symmetry=1, barrier=(20.0898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33239,'amu*angstrom^2'), symmetry=1, barrier=(30.6343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87466,'amu*angstrom^2'), symmetry=1, barrier=(20.1102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74478,0.138239,-0.00023481,2.02849e-07,-6.83542e-11,-121141,38.4855], Tmin=(100,'K'), Tmax=(813.917,'K')), NASAPolynomial(coeffs=[15.9172,0.0364273,-1.95133e-05,3.84239e-09,-2.67913e-13,-123519,-40.0324], Tmin=(813.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1008.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    E0 = (-359.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-32.9486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-202.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (86.5294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-226.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (58.7936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (4.91432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-351.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-334.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-102.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-301.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-257.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-292.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-84.2785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-280.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-178.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-320.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-273.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-79.9265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-191.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-87.9206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (21.3948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-110.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-231.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-284.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-296.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-290.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (50.0899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-160.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-12.7077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-194.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (20.4568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-163.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=CC(=O)F(4234)', 'FC(F)=CC(F)F(344)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C(F)(F)[CH]C(F)F(6773)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(88.9,'m^3/(mol*s)'), n=1.51, Ea=(338.188,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C(F)(C=O)C([C](F)F)C(F)F(6758)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=COF)C(F)(F)[CH]C(F)F(7470)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(199.087,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(6759)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(132.515,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(F)F-2(957)', '[O]C(F)(C=O)[C](F)F(4474)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', '[O]C=C(F)C(F)(F)[CH]C(F)F(7471)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=CC1(F)OC(C(F)F)C1(F)F(6766)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=CC(O)(F)C(F)(F)C=C(F)F(6772)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['FC(F)[CH]C(F)(F)C1(F)[CH]OO1(7472)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C1(F)[CH]OC(C(F)F)C1(F)F(7412)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C1OC1(F)C(F)(F)[CH]C(F)F(7473)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C1C(C(F)F)C(F)(F)C1([O])F(7474)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(67.0567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 63.5 to 67.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]C(F)(C=O)C(F)=CC(F)F(7475)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(52.8062,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(398)', 'FC(F)=CC(F)F(344)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(48.8108,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=CC(=O)C(F)(F)[CH]C(F)F(7476)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(17.388,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=CC(=O)F(4234)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.13223e-11,'m^3/(mol*s)'), n=4.48095, Ea=(38.9319,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HCO(15)', 'O=C(F)C(F)(F)[CH]C(F)F(6995)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(108.759,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F(37)', '[O]C(F)(C=O)C(F)(F)C=CF(7477)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.1183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(F)(C=O)C(F)(F)C=C(F)F(7478)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.5227,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][CH]C(=O)F(398)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(2.98938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F2(78)', '[O]C=C([O])C(F)=CC(F)F(7479)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(26.8179,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=[C]C(=O)C(F)(F)[CH]C(F)F(7480)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(279.94,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C(F)(C=O)C(F)(F)C[C](F)F(7481)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=[C]C(O)(F)C(F)(F)[CH]C(F)F(7482)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C(F)([C]=O)C(F)(F)CC(F)F(7483)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.44586,'s^-1'), n=3.07722, Ea=(63.2598,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['O=CC(O)(F)C(F)(F)[CH][C](F)F(7484)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.54863e+08,'s^-1'), n=0.802517, Ea=(68.5723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_noH] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=CC(F)(OF)[C](F)[CH]C(F)F(7485)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.683,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)(C=O)[C](F)C(F)C(F)F(7486)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(166.592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C=C(OF)C(F)(F)[CH]C(F)F(7487)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(96.1923,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    products = ['[O]C=C([O])C(F)(F)C(F)C(F)F(7488)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(165.248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=CC(F)(OF)C(F)(F)[CH][CH]F(7489)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(80.6192,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(F)(C=O)C(F)(F)C(F)[CH]F(7490)'],
    products = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2073',
    isomers = [
        '[O]C(F)(C=O)C(F)(F)[CH]C(F)F(6757)',
    ],
    reactants = [
        ('O=CC(=O)F(4234)', 'FC(F)=CC(F)F(344)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2073',
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

