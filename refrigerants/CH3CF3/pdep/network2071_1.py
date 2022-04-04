species(
    label = '[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1078.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,1380,1390,370,380,2900,435,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.694023,0.112303,-0.000163658,1.26382e-07,-3.9302e-11,-129493,41.7439], Tmin=(100,'K'), Tmax=(784.554,'K')), NASAPolynomial(coeffs=[14.0869,0.0369402,-1.95635e-05,3.93338e-09,-2.81497e-13,-131813,-25.9785], Tmin=(784.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1078.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
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
    label = '[O]C(C(=O)F)C([C](F)F)C(F)F(6756)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {8,S} {15,S}
11 C u1 p0 c0 {3,S} {4,S} {8,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1061.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4125.71,'J/mol'), sigma=(6.24249,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=644.43 K, Pc=38.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.548789,0.107475,-0.00014503,1.01968e-07,-2.87629e-11,-127543,41.3798], Tmin=(100,'K'), Tmax=(863.359,'K')), NASAPolynomial(coeffs=[15.2927,0.03408,-1.75145e-05,3.50337e-09,-2.51069e-13,-130278,-32.719], Tmin=(863.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1061.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'O=C[C](F)OC(F)(F)[CH]C(F)F(6761)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u0 p0 c0 {7,D} {11,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1114.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,280,501,1494,1531,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3622.93,'J/mol'), sigma=(6.02473,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.89 K, Pc=37.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07835,0.124761,-0.000210358,1.9162e-07,-6.8967e-11,-133835,37.6305], Tmin=(100,'K'), Tmax=(792.149,'K')), NASAPolynomial(coeffs=[10.6065,0.0464049,-2.53374e-05,5.06736e-09,-3.57973e-13,-135079,-12.1867], Tmin=(792.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1114.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C([C](F)F)C(=O)F(4487)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-591.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,190,488,555,1236,1407,486,617,768,1157,1926,378.555,378.585,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0876046,'amu*angstrom^2'), symmetry=1, barrier=(8.91333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259265,'amu*angstrom^2'), symmetry=1, barrier=(26.381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43774,0.0608417,-8.50182e-05,6.14589e-08,-1.77788e-11,-71047.5,27.7795], Tmin=(100,'K'), Tmax=(843.24,'K')), NASAPolynomial(coeffs=[10.3543,0.0185442,-9.77556e-06,1.97086e-09,-1.41659e-13,-72551.2,-13.7173], Tmin=(843.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-591.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(CsCsF1sF1s)"""),
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
    label = 'O=C(F)[CH]C(F)(F)[CH]C(F)F(7509)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u1 p0 c0 {7,S} {11,S} {14,S}
11 C u0 p0 c0 {5,S} {6,D} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-989.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,522,611,926,1093,1137,1374,1416,3112,3000,3050,390,425,1340,1360,335,370,611,648,830,1210,1753,268.482,268.482,268.482,2178.19],'cm^-1')),
        HinderedRotor(inertia=(0.141984,'amu*angstrom^2'), symmetry=1, barrier=(7.26265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141984,'amu*angstrom^2'), symmetry=1, barrier=(7.26266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668735,'amu*angstrom^2'), symmetry=1, barrier=(34.2066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668735,'amu*angstrom^2'), symmetry=1, barrier=(34.2066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368339,0.104953,-0.000154034,1.21698e-07,-3.89678e-11,-118888,35.3635], Tmin=(100,'K'), Tmax=(760.901,'K')), NASAPolynomial(coeffs=[12.5977,0.0367896,-1.9658e-05,3.96132e-09,-2.83702e-13,-120861,-23.6468], Tmin=(760.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-989.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(COCsFO) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCJC=O)"""),
)

species(
    label = 'O=C(F)C1OC(C(F)F)C1(F)F(6764)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {6,S} {9,S} {12,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1338.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.501212,0.0858352,-7.54074e-05,3.11698e-08,-5.01663e-12,-160770,35.9666], Tmin=(100,'K'), Tmax=(1501.06,'K')), NASAPolynomial(coeffs=[24.312,0.0197129,-9.33141e-06,1.82331e-09,-1.28971e-13,-168219,-93.8212], Tmin=(1501.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1338.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + ring(Cs-O2s-Cs-Cs(F)(F))"""),
)

species(
    label = 'O=C(F)C(=O)C(F)(F)CC(F)F(7510)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,D}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
11 C u0 p0 c0 {6,D} {9,S} {12,S}
12 C u0 p0 c0 {5,S} {7,D} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1410.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.816187,0.117011,-0.000174296,1.33627e-07,-3.91375e-11,-169448,35.0924], Tmin=(100,'K'), Tmax=(643.678,'K')), NASAPolynomial(coeffs=[13.699,0.0390019,-2.09208e-05,4.20314e-09,-2.99736e-13,-171569,-30.5021], Tmin=(643.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1410.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO)"""),
)

species(
    label = 'O=C(F)C(O)C(F)(F)C=C(F)F(6771)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {15,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {9,S} {12,D} {14,S}
11 C u0 p0 c0 {3,S} {7,D} {8,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1367.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.38371,0.103658,-0.000132091,8.69695e-08,-2.30687e-11,-164306,37.4851], Tmin=(100,'K'), Tmax=(913.609,'K')), NASAPolynomial(coeffs=[15.5065,0.0340885,-1.78704e-05,3.6235e-09,-2.6231e-13,-167210,-37.7405], Tmin=(913.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1367.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF)"""),
)

species(
    label = 'F[C]1OOC1C(F)(F)[CH]C(F)F(7511)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
11 C u1 p0 c0 {9,S} {10,S} {15,S}
12 C u1 p0 c0 {5,S} {7,S} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-774.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.847852,0.116595,-0.000185055,1.59512e-07,-5.52081e-11,-93035.8,39.1562], Tmin=(100,'K'), Tmax=(762.731,'K')), NASAPolynomial(coeffs=[12.2838,0.040583,-2.15136e-05,4.28434e-09,-3.03118e-13,-94831.1,-19.2769], Tmin=(762.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-774.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1[C](F)OC(C(F)F)C1(F)F(7340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {12,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1062.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.792993,0.104064,-0.000125122,7.72301e-08,-1.87072e-11,-127656,30.4491], Tmin=(100,'K'), Tmax=(1014.06,'K')), NASAPolynomial(coeffs=[18.8571,0.0265536,-1.04694e-05,1.85498e-09,-1.24745e-13,-131641,-64.6259], Tmin=(1014.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1062.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(CsCsFFH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)OC1C(F)(F)[CH]C(F)F(7512)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
11 C u0 p0 c0 {4,S} {5,S} {12,S} {14,S}
12 C u1 p0 c0 {9,S} {11,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-955.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.484421,0.107719,-0.00014868,1.1108e-07,-3.38438e-11,-114800,36.191], Tmin=(100,'K'), Tmax=(796.386,'K')), NASAPolynomial(coeffs=[13.0048,0.0399651,-2.10598e-05,4.24353e-09,-3.04811e-13,-116948,-25.8151], Tmin=(796.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-955.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFOO) + group(CsCsFFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C1C([O])(F)C(C(F)F)C1(F)F(7457)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u1 p2 c0 {11,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
11 C u0 p0 c0 {7,S} {9,S} {10,S} {14,S}
12 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-993.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.269937,0.0965548,-0.000103697,5.59253e-08,-1.21141e-11,-119379,32.288], Tmin=(100,'K'), Tmax=(1108.78,'K')), NASAPolynomial(coeffs=[17.6856,0.0317793,-1.60663e-05,3.23655e-09,-2.34282e-13,-123361,-56.1911], Tmin=(1108.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-993.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCsFFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(O2sj(Cs-F1sCsCs)) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C([C](O)F)C(F)(F)[CH]C(F)F(7513)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {12,S} {15,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {7,D} {8,S} {12,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1128.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,375,552.5,462.5,1710,280,501,1494,1531,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09362,0.122593,-0.000200096,1.73072e-07,-5.92385e-11,-135540,39.0076], Tmin=(100,'K'), Tmax=(795.015,'K')), NASAPolynomial(coeffs=[13.1845,0.0390899,-2.05381e-05,4.04668e-09,-2.83564e-13,-137442,-24.2823], Tmin=(795.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1128.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C(C(=O)F)C(F)=CC(F)F(7514)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 C u0 p0 c0 {4,S} {6,D} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-919.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,281.287,281.287,1718.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.184558,'amu*angstrom^2'), symmetry=1, barrier=(10.3625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18456,'amu*angstrom^2'), symmetry=1, barrier=(10.3625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440595,'amu*angstrom^2'), symmetry=1, barrier=(24.7383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162553,0.0919603,-0.000122997,8.97081e-08,-2.6767e-11,-110409,36.423], Tmin=(100,'K'), Tmax=(811.67,'K')), NASAPolynomial(coeffs=[11.6671,0.0352655,-1.82242e-05,3.65382e-09,-2.62053e-13,-112277,-16.6792], Tmin=(811.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-919.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(COCsFO) + radical(C=OCOJ)"""),
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
    label = 'CFO(51)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC(F)(F)[CH]C(F)F(6652)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u0 p0 c0 {5,D} {6,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-838.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,230.583,230.583],'cm^-1')),
        HinderedRotor(inertia=(0.135692,'amu*angstrom^2'), symmetry=1, barrier=(5.1196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135691,'amu*angstrom^2'), symmetry=1, barrier=(5.1196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375337,'amu*angstrom^2'), symmetry=1, barrier=(14.1614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851186,0.0753436,-9.67733e-05,6.71648e-08,-1.90871e-11,-100781,28.1501], Tmin=(100,'K'), Tmax=(849.898,'K')), NASAPolynomial(coeffs=[10.6965,0.0290068,-1.49919e-05,3.01427e-09,-2.16881e-13,-102454,-17.7465], Tmin=(849.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCC=O)"""),
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
    label = 'O=C(F)C(=O)C(F)(F)[CH]C(F)F(7515)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,D}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {6,D} {8,S} {12,S}
12 C u0 p0 c0 {5,S} {7,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1210.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,375,552.5,462.5,1710,286,619,818,1246,1924,199.056,199.057,199.057,199.058],'cm^-1')),
        HinderedRotor(inertia=(0.00425446,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00425446,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293575,'amu*angstrom^2'), symmetry=1, barrier=(8.25475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890586,'amu*angstrom^2'), symmetry=1, barrier=(25.0416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.739894,0.113542,-0.000173907,1.39216e-07,-4.46034e-11,-145406,37.2955], Tmin=(100,'K'), Tmax=(763.693,'K')), NASAPolynomial(coeffs=[14.3402,0.0345568,-1.87708e-05,3.79083e-09,-2.71193e-13,-147710,-31.3919], Tmin=(763.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1210.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO) + radical(CCJCC=O)"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C=CF(7516)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 C u0 p0 c0 {4,S} {9,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-922.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,194,682,905,1196,1383,3221,180,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.580893,'amu*angstrom^2'), symmetry=1, barrier=(13.3559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27417,'amu*angstrom^2'), symmetry=1, barrier=(29.2956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254791,'amu*angstrom^2'), symmetry=1, barrier=(29.3921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101818,0.0926599,-0.000118596,8.02341e-08,-2.20294e-11,-110830,35.514], Tmin=(100,'K'), Tmax=(881.475,'K')), NASAPolynomial(coeffs=[13.2524,0.0329849,-1.70475e-05,3.43255e-09,-2.47355e-13,-113149,-26.2708], Tmin=(881.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-922.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C=C(F)F(7517)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {6,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {8,S} {12,D} {14,S}
11 C u0 p0 c0 {3,S} {7,D} {9,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1123.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,182,240,577,636,1210,1413,297.1,297.102,297.102,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00190983,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0877576,'amu*angstrom^2'), symmetry=1, barrier=(5.49692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387177,'amu*angstrom^2'), symmetry=1, barrier=(24.2517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288349,0.102472,-0.000142096,1.0342e-07,-3.03425e-11,-134996,37.9216], Tmin=(100,'K'), Tmax=(829.304,'K')), NASAPolynomial(coeffs=[13.9405,0.0338427,-1.79652e-05,3.63396e-09,-2.61708e-13,-137356,-28.0613], Tmin=(829.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1123.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)C(=O)[C](F)[CH]C(F)F(7518)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u1 p0 c0 {7,S} {9,S} {13,S}
9  C u1 p0 c0 {3,S} {8,S} {10,S}
10 C u0 p0 c0 {5,D} {9,S} {11,S}
11 C u0 p0 c0 {4,S} {6,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-838.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,234,347,1316,1464,375,552.5,462.5,1710,286,619,818,1246,1924,327.389,327.404,327.406,327.413],'cm^-1')),
        HinderedRotor(inertia=(0.532341,'amu*angstrom^2'), symmetry=1, barrier=(40.4954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00157235,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230557,'amu*angstrom^2'), symmetry=1, barrier=(17.5385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859558,'amu*angstrom^2'), symmetry=1, barrier=(6.53876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0482577,0.0966635,-0.000136305,1.00653e-07,-2.99006e-11,-100696,34.8302], Tmin=(100,'K'), Tmax=(820.152,'K')), NASAPolynomial(coeffs=[13.3709,0.0312157,-1.66046e-05,3.35273e-09,-2.40973e-13,-102897,-27.2487], Tmin=(820.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO) + radical(CCJCC=O) + radical(CsCOCsF1s)"""),
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
    label = 'O=C(F)[C](O)C(F)(F)[CH]C(F)F(7519)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,S} {15,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u1 p0 c0 {6,S} {8,S} {12,S}
12 C u0 p0 c0 {5,S} {7,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1145.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.845401,0.120028,-0.000176252,1.18275e-07,-2.24479e-11,-137563,40.9878], Tmin=(100,'K'), Tmax=(592.661,'K')), NASAPolynomial(coeffs=[14.7521,0.0364601,-1.96723e-05,3.93081e-09,-2.78234e-13,-139793,-29.3171], Tmin=(592.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1145.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(C2CsJOH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C[C](F)F(7520)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
11 C u1 p0 c0 {3,S} {4,S} {10,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1085.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07893,0.121963,-0.000198122,1.70356e-07,-5.82372e-11,-130331,41.3701], Tmin=(100,'K'), Tmax=(776.517,'K')), NASAPolynomial(coeffs=[13.58,0.0381947,-2.0358e-05,4.0468e-09,-2.85402e-13,-132358,-24.0402], Tmin=(776.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1085.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)=C([O])C(F)(F)CC(F)F(7521)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {11,S}
7  O u1 p2 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
11 C u0 p0 c0 {6,S} {9,S} {12,D}
12 C u0 p0 c0 {5,S} {7,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1159.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33379,0.132711,-0.000234945,2.16907e-07,-7.75202e-11,-139317,36.7975], Tmin=(100,'K'), Tmax=(820.426,'K')), NASAPolynomial(coeffs=[10.8222,0.0456429,-2.4929e-05,4.94749e-09,-3.46319e-13,-140375,-13.7399], Tmin=(820.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1159.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C(F)C(O)C(F)(F)[CH][C](F)F(7522)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {15,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u1 p0 c0 {9,S} {12,S} {14,S}
11 C u0 p0 c0 {3,S} {7,D} {8,S}
12 C u1 p0 c0 {4,S} {5,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1121.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.760415,0.116132,-0.000173094,1.29789e-07,-3.59264e-11,-134688,42.5852], Tmin=(100,'K'), Tmax=(631.002,'K')), NASAPolynomial(coeffs=[13.8773,0.037586,-2.0237e-05,4.06666e-09,-2.8981e-13,-136819,-23.5406], Tmin=(631.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1121.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(F)C(OF)[C](F)[CH]C(F)F(7523)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {11,S} {14,S}
10 C u1 p0 c0 {3,S} {8,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {15,S}
12 C u0 p0 c0 {4,S} {7,D} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-767.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,3025,407.5,1350,352.5,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.759719,0.114933,-0.000176314,1.47427e-07,-5.01378e-11,-92093.4,43.7954], Tmin=(100,'K'), Tmax=(716.958,'K')), NASAPolynomial(coeffs=[12.3897,0.0415562,-2.27652e-05,4.62026e-09,-3.31579e-13,-93978.5,-15.2651], Tmin=(716.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O]C([C](F)C(F)C(F)F)C(=O)F(7524)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {6,S} {11,S} {12,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {15,S}
11 C u1 p0 c0 {4,S} {8,S} {9,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1040.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.842852,0.115508,-0.000174861,1.40621e-07,-4.53692e-11,-125025,42.1612], Tmin=(100,'K'), Tmax=(758.222,'K')), NASAPolynomial(coeffs=[14.0666,0.0368591,-1.92802e-05,3.83629e-09,-2.72137e-13,-127286,-25.6426], Tmin=(758.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1040.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(COCsFO) + radical(C=OCOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = 'O=C(F)C(OF)C(F)(F)[CH][CH]F(7525)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u1 p0 c0 {9,S} {12,S} {14,S}
11 C u0 p0 c0 {3,S} {7,D} {8,S}
12 C u1 p0 c0 {4,S} {10,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-778.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,486,617,768,1157,1926,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28035,0.126947,-0.000207905,1.79197e-07,-6.15263e-11,-93493.1,42.8876], Tmin=(100,'K'), Tmax=(761.086,'K')), NASAPolynomial(coeffs=[14.2071,0.0389089,-2.13041e-05,4.27992e-09,-3.0373e-13,-95658.2,-26.3387], Tmin=(761.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-778.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(COCsFO) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C(F)[CH]F(7526)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
11 C u1 p0 c0 {4,S} {10,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1049.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1672,0.12288,-0.000188671,1.49692e-07,-4.72027e-11,-126056,40.8139], Tmin=(100,'K'), Tmax=(777.74,'K')), NASAPolynomial(coeffs=[16.0083,0.0345443,-1.83015e-05,3.65394e-09,-2.59504e-13,-128728,-37.7304], Tmin=(777.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1049.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'O=[C]C(OF)C(F)(F)[CH]C(F)F(7527)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u1 p0 c0 {7,D} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-773.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,1380,1390,370,380,2900,435,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01432,0.119098,-0.000173619,1.29284e-07,-3.82815e-11,-92831.1,42.1817], Tmin=(100,'K'), Tmax=(826.351,'K')), NASAPolynomial(coeffs=[16.7154,0.0332768,-1.78346e-05,3.60384e-09,-2.58788e-13,-95761.3,-39.9721], Tmin=(826.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-773.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)C(F)(F)C(F)C(F)F(7528)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {6,S} {8,S} {12,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {9,S} {15,S}
12 C u1 p0 c0 {7,D} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1044.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,1855,455,950,347.035,347.04,347.042,347.042],'cm^-1')),
        HinderedRotor(inertia=(0.00139971,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130045,'amu*angstrom^2'), symmetry=1, barrier=(11.1141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130043,'amu*angstrom^2'), symmetry=1, barrier=(11.1141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130042,'amu*angstrom^2'), symmetry=1, barrier=(11.1141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03609,0.11812,-0.000170344,1.24518e-07,-3.59544e-11,-125409,40.6011], Tmin=(100,'K'), Tmax=(849.423,'K')), NASAPolynomial(coeffs=[17.5663,0.03052,-1.565e-05,3.10755e-09,-2.21075e-13,-128569,-46.1087], Tmin=(849.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1044.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (-355.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (44.0221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-179.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-211.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (62.7145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-24.0952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-347.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-292.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-330.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-52.1658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-297.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-233.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-271.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-158.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-73.3273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-266.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-293.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-232.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-259.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-76.0055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-187.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-46.2715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-149.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-86.7938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-201.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-227.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-211.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-286.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (57.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-154.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (24.3777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-159.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (53.0219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-96.2312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=CC(=O)F(4234)', 'FC(F)=CC(F)F(344)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C(F)(F)[CH]C(F)F(6773)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(370.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(C(=O)F)C([C](F)F)C(F)F(6756)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C[C](F)OC(F)(F)[CH]C(F)F(6761)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(144.017,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(F)F-2(957)', '[O]C([C](F)F)C(=O)F(4487)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', 'O=C(F)[CH]C(F)(F)[CH]C(F)F(7509)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C(F)C1OC(C(F)F)C1(F)F(6764)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C(F)C(=O)C(F)(F)CC(F)F(7510)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C(F)C(O)C(F)(F)C=C(F)F(6771)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['F[C]1OOC1C(F)(F)[CH]C(F)F(7511)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['[O]C1[C](F)OC(C(F)F)C1(F)F(7340)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['[O]C1(F)OC1C(F)(F)[CH]C(F)F(7512)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(122.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 122.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['[O]C1C([O])(F)C(C(F)F)C1(F)F(7457)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(84.1695,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_csHNd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 82.7 to 84.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C([C](O)F)C(F)(F)[CH]C(F)F(7513)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(246.927,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(C(=O)F)C(F)=CC(F)F(7514)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(50.2687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][CH]C(=O)F(398)', 'FC(F)=CC(F)F(344)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(18.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=CC(=O)F(4234)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.13223e-11,'m^3/(mol*s)'), n=4.48095, Ea=(21.3059,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(51)', 'O=CC(F)(F)[CH]C(F)F(6652)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(74.2259,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C(F)C(=O)C(F)(F)[CH]C(F)F(7515)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(16.0861,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F(37)', '[O]C(C(=O)F)C(F)(F)C=CF(7516)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.1183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]C(C(=O)F)C(F)(F)C=C(F)F(7517)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.5227,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][CH]C(=O)F(398)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=C(F)C(=O)[C](F)[CH]C(F)F(7518)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(247.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=[C]C(=O)C(F)(F)[CH]C(F)F(7480)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(258.587,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C(F)[C](O)C(F)(F)[CH]C(F)F(7519)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['[O]C(C(=O)F)C(F)(F)C[C](F)F(7520)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['[O]C(F)=C([O])C(F)(F)CC(F)F(7521)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    products = ['O=C(F)C(O)C(F)(F)[CH][C](F)F(7522)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.54863e+08,'s^-1'), n=0.802517, Ea=(68.5723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_noH] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=C(F)C(OF)[C](F)[CH]C(F)F(7523)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(102.054,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C([C](F)C(F)C(F)F)C(=O)F(7524)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(164.127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=C(F)C(OF)C(F)(F)[CH][CH]F(7525)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(80.6192,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C(C(=O)F)C(F)(F)C(F)[CH]F(7526)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=[C]C(OF)C(F)(F)[CH]C(F)F(7527)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C([C]=O)C(F)(F)C(F)C(F)F(7528)'],
    products = ['[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(225.329,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

network(
    label = 'PDepNetwork #2071',
    isomers = [
        '[O]C(C(=O)F)C(F)(F)[CH]C(F)F(6755)',
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
    label = 'PDepNetwork #2071',
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

