species(
    label = 'O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31452,0.129622,-0.000200363,1.51022e-07,-4.13613e-11,-131915,38.5594], Tmin=(100,'K'), Tmax=(632.916,'K')), NASAPolynomial(coeffs=[16.1812,0.0362794,-1.99736e-05,4.02367e-09,-2.86276e-13,-134475,-40.5707], Tmin=(632.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1098.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CC(=O)F(335)',
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
    label = 'FC(F)=CC(F)F(2834)',
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
    label = '[CH]C(=O)F(8932)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.05,1655.52],'cm^-1')),
        HinderedRotor(inertia=(0.0191536,'amu*angstrom^2'), symmetry=1, barrier=(5.31745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10758e-05,3.91582e-09,-5.6142e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.41,'K')), NASAPolynomial(coeffs=[6.53552,0.00717479,-3.36979e-06,6.65147e-10,-4.72043e-14,-1600.48,-4.84321], Tmin=(1580.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C(F)(F)[CH]C(F)F(2594)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-700.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0263176,'amu*angstrom^2'), symmetry=1, barrier=(12.676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292203,'amu*angstrom^2'), symmetry=1, barrier=(6.71831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15905,0.0685176,-9.91319e-05,7.75163e-08,-2.45999e-11,-84188.1,25.196], Tmin=(100,'K'), Tmax=(766.908,'K')), NASAPolynomial(coeffs=[9.60201,0.0244801,-1.29966e-05,2.63779e-09,-1.90099e-13,-85483.1,-13.2955], Tmin=(766.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-700.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[CH]C(F)F-2(3493)',
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
    label = 'O=C(F)[CH]O[C](F)F(1755)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {4,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-610.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,493,600,700,1144,1293,244.153,244.166,244.198,244.198],'cm^-1')),
        HinderedRotor(inertia=(0.237006,'amu*angstrom^2'), symmetry=1, barrier=(10.0272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235805,'amu*angstrom^2'), symmetry=1, barrier=(10.0269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645282,'amu*angstrom^2'), symmetry=1, barrier=(27.298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888445,0.0741621,-0.000124068,1.048e-07,-3.46432e-11,-73370.2,26.6211], Tmin=(100,'K'), Tmax=(800.21,'K')), NASAPolynomial(coeffs=[11.1627,0.0179605,-9.63843e-06,1.90235e-09,-1.32985e-13,-74859.4,-19.6874], Tmin=(800.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=C(F)C1OC(F)(F)C1C(F)F(8981)',
    structure = adjacencyList("""1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1392.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.384122,0.0845571,-7.44393e-05,3.11943e-08,-5.11092e-12,-167315,36.2417], Tmin=(100,'K'), Tmax=(1472.24,'K')), NASAPolynomial(coeffs=[23.0064,0.021006,-9.68947e-06,1.87384e-09,-1.31999e-13,-174202,-85.6508], Tmin=(1472.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1392.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(CsCsFFH) + group(COCsFO) + ring(Cs-Cs(F)(F)-O2s-Cs)"""),
)

species(
    label = 'O=C(F)COC(F)(F)C=C(F)F(8987)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {6,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {12,D} {15,S}
11 C u0 p0 c0 {3,S} {7,D} {9,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1344.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.619691,0.110787,-0.000155975,1.17108e-07,-3.56201e-11,-161527,35.0149], Tmin=(100,'K'), Tmax=(799.53,'K')), NASAPolynomial(coeffs=[13.8188,0.0385569,-2.04746e-05,4.13258e-09,-2.9692e-13,-163836,-31.413], Tmin=(799.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1344.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF)"""),
)

species(
    label = 'F[C]1OC1OC(F)(F)[CH]C(F)F(9213)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
11 C u1 p0 c0 {9,S} {10,S} {15,S}
12 C u1 p0 c0 {5,S} {7,S} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-948.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.689113,0.113454,-0.000160009,1.21014e-07,-3.7321e-11,-113884,35.7909], Tmin=(100,'K'), Tmax=(786.72,'K')), NASAPolynomial(coeffs=[13.4739,0.0414373,-2.2685e-05,4.6358e-09,-3.35426e-13,-116112,-29.1388], Tmin=(786.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-948.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsCFHO) + group(CsCsFFH) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'F[C]1[CH]OC(F)(F)C(C(F)F)O1(9030)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {7,S} {11,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1096.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56655,0.117226,-0.000155015,1.02405e-07,-2.5407e-11,-131717,31.0796], Tmin=(100,'K'), Tmax=(1135.83,'K')), NASAPolynomial(coeffs=[21.4601,0.019154,-3.07499e-06,6.31587e-11,1.57069e-14,-135853,-78.1221], Tmin=(1135.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1096.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(1,4-Dioxane) + radical(CsCsF1sO2s) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)(F)C1C(F)F(9214)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
12 C u1 p0 c0 {6,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1104.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0136,0.103853,-0.000116693,6.26727e-08,-1.20997e-11,-132600,30.2274], Tmin=(100,'K'), Tmax=(943.704,'K')), NASAPolynomial(coeffs=[22.1496,0.0201056,-6.51916e-06,1.04826e-09,-6.77923e-14,-137615,-83.5847], Tmin=(943.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1104.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(CsCFFO) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Tetrahydrofuran) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOCs)"""),
)

species(
    label = 'F[C](F)[CH]C(F)F(4845)',
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
    label = 'O=C(F)[CH]OC(F)=CC(F)F(9215)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 C u1 p0 c0 {5,S} {11,S} {14,S}
11 C u0 p0 c0 {4,S} {6,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-927.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.565244,0.105572,-0.000138334,8.97284e-08,-2.29653e-11,-111375,33.6755], Tmin=(100,'K'), Tmax=(954.375,'K')), NASAPolynomial(coeffs=[18.3262,0.0264038,-1.39201e-05,2.8319e-09,-2.05571e-13,-114982,-56.585], Tmin=(954.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CsCFFH) + group(Cds-CdsCsH) + group(CdCFO) + group(COCsFO) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]C(=O)F(509)',
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
    label = 'O=C(F)[CH]OC(F)(F)C=CF(9216)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u1 p0 c0 {5,S} {11,S} {13,S}
10 C u0 p0 c0 {4,S} {8,D} {14,S}
11 C u0 p0 c0 {3,S} {6,D} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-962.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.632421,0.109242,-0.000156232,1.13826e-07,-3.2902e-11,-115642,33.4953], Tmin=(100,'K'), Tmax=(846.87,'K')), NASAPolynomial(coeffs=[16.177,0.0298479,-1.56094e-05,3.12765e-09,-2.23766e-13,-118489,-44.8066], Tmin=(846.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-962.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(CCsJOCs)"""),
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
    label = 'O=C(F)[CH]OC(F)(F)C=C(F)F(9217)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {11,D} {13,S}
10 C u1 p0 c0 {6,S} {12,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {9,D}
12 C u0 p0 c0 {3,S} {7,D} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1163.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,182,240,577,636,1210,1413,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0262,0.119133,-0.000180185,1.37861e-07,-4.16994e-11,-139808,35.9136], Tmin=(100,'K'), Tmax=(811.278,'K')), NASAPolynomial(coeffs=[16.9353,0.0305736,-1.64453e-05,3.30874e-09,-2.36369e-13,-142722,-46.9838], Tmin=(811.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1163.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C=COC(F)(F)[CH]C(F)F(9218)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-897.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983578,0.111943,-0.000147532,9.37306e-08,-2.31487e-11,-107781,33.1205], Tmin=(100,'K'), Tmax=(995.507,'K')), NASAPolynomial(coeffs=[21.5955,0.0212181,-1.08285e-05,2.18215e-09,-1.58017e-13,-112276,-75.7083], Tmin=(995.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-897.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'O=C(F)[CH]OC(F)(F)C[C](F)F(9219)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {9,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u1 p0 c0 {6,S} {12,S} {15,S}
11 C u1 p0 c0 {3,S} {4,S} {8,S}
12 C u0 p0 c0 {5,S} {7,D} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1107.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25788,0.127924,-0.000196359,1.46539e-07,-3.94572e-11,-133023,38.192], Tmin=(100,'K'), Tmax=(636.522,'K')), NASAPolynomial(coeffs=[16.2518,0.0352927,-1.90785e-05,3.8141e-09,-2.70145e-13,-135605,-41.1414], Tmin=(636.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1107.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(F)COC(F)(F)[CH][C](F)F(9220)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {6,S} {11,S} {13,S} {14,S}
10 C u1 p0 c0 {8,S} {12,S} {15,S}
11 C u0 p0 c0 {3,S} {7,D} {9,S}
12 C u1 p0 c0 {4,S} {5,S} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1078.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,486,617,768,1157,1926,190,488,555,1236,1407,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29176,0.130149,-0.000225493,2.0513e-07,-7.30419e-11,-129501,40.3424], Tmin=(100,'K'), Tmax=(801.972,'K')), NASAPolynomial(coeffs=[11.801,0.0436726,-2.41473e-05,4.83424e-09,-3.40893e-13,-130920,-15.6883], Tmin=(801.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1078.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(Cs-(Cds-O2d)OsHH) + group(CsCsFFH) + group(COCsFO) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(F)[CH]O[C](F)C(F)C(F)F(9221)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {14,S}
10 C u1 p0 c0 {4,S} {6,S} {8,S}
11 C u1 p0 c0 {6,S} {12,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1056.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,395,473,707,1436,3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58986,0.132513,-0.000211732,1.71939e-07,-5.50604e-11,-126875,39.3145], Tmin=(100,'K'), Tmax=(768.226,'K')), NASAPolynomial(coeffs=[17.5029,0.0331017,-1.76298e-05,3.49964e-09,-2.46742e-13,-129808,-47.7633], Tmin=(768.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1056.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(CsCsF1sO2s) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(F)C(F)O[C](F)[CH]C(F)F(9222)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {14,S}
10 C u1 p0 c0 {9,S} {11,S} {15,S}
11 C u1 p0 c0 {4,S} {6,S} {10,S}
12 C u0 p0 c0 {5,S} {7,D} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1080.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,395,473,707,1436,486,617,768,1157,1926,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14596,0.12432,-0.000200451,1.74406e-07,-6.12199e-11,-129824,39.223], Tmin=(100,'K'), Tmax=(733.205,'K')), NASAPolynomial(coeffs=[12.9152,0.0427931,-2.38079e-05,4.83319e-09,-3.45844e-13,-131756,-23.3674], Tmin=(733.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1080.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(COCsFO) + radical(Csj(Cs-F1sO2sH)(Cs-F1sF1sH)(H)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C(F)[CH]OC(F)(F)C(F)[CH]F(9223)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u1 p0 c0 {6,S} {12,S} {14,S}
11 C u1 p0 c0 {4,S} {8,S} {15,S}
12 C u0 p0 c0 {5,S} {7,D} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1071.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,3025,407.5,1350,352.5,334,575,1197,1424,3202,611,648,830,1210,1753,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62511,0.132895,-0.000205649,1.59321e-07,-4.8531e-11,-128737,38.5934], Tmin=(100,'K'), Tmax=(807.434,'K')), NASAPolynomial(coeffs=[18.8983,0.0312198,-1.67558e-05,3.35421e-09,-2.38424e-13,-132051,-56.0294], Tmin=(807.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1071.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(F)C(F)OC(F)(F)[CH][CH]F(9224)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {11,S} {13,S}
10 C u1 p0 c0 {8,S} {12,S} {14,S}
11 C u0 p0 c0 {4,S} {7,D} {9,S}
12 C u1 p0 c0 {5,S} {10,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1099.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48561,0.131754,-0.000218078,1.87809e-07,-6.39936e-11,-132077,40.618], Tmin=(100,'K'), Tmax=(779.409,'K')), NASAPolynomial(coeffs=[14.9456,0.0383486,-2.08438e-05,4.15997e-09,-2.93584e-13,-134363,-32.7892], Tmin=(779.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1099.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + group(COCsFO) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O][C]=COC(F)(F)C(F)C(F)F(9225)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u1 p2 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
11 C u0 p0 c0 {6,S} {12,D} {15,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1080.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05351,0.121533,-0.000147478,8.29539e-08,-1.76963e-11,-129682,39.9909], Tmin=(100,'K'), Tmax=(1165.7,'K')), NASAPolynomial(coeffs=[30.5352,0.00970821,-3.5842e-06,6.6088e-10,-4.75347e-14,-137280,-122.227], Tmin=(1165.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1080.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]C(F)OC(F)(F)[CH]C(F)F(9099)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {11,S} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {12,S} {14,S}
11 C u1 p0 c0 {8,S} {9,S} {15,S}
12 C u1 p0 c0 {7,D} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1094.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,355,410,600,1181,1341,1420,3056,3025,407.5,1350,352.5,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68443,0.139958,-0.000249043,2.24803e-07,-7.83129e-11,-131492,39.2667], Tmin=(100,'K'), Tmax=(826.678,'K')), NASAPolynomial(coeffs=[13.5283,0.0406831,-2.23373e-05,4.42147e-09,-3.08319e-13,-133130,-25.926], Tmin=(826.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1094.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(F)(C=O)C(F)(F)[CH]C(F)F(8975)',
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
    label = 'F[C]=COC(F)(F)[CH]C(F)F(9226)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-776.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00344,0.118331,-0.000174847,1.29557e-07,-3.7905e-11,-93199,34.3613], Tmin=(100,'K'), Tmax=(838.185,'K')), NASAPolynomial(coeffs=[17.594,0.0295839,-1.60344e-05,3.24824e-09,-2.33435e-13,-96316.8,-52.0786], Tmin=(838.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-776.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FC1=COC(F)(F)C(C(F)F)O1(9043)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {7,S} {11,D} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1337.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26538,0.100755,-0.000109224,5.7867e-08,-1.16788e-11,-160685,26.6412], Tmin=(100,'K'), Tmax=(1310.99,'K')), NASAPolynomial(coeffs=[24.8127,0.014822,-3.61895e-06,4.6097e-10,-2.5432e-14,-166976,-104.145], Tmin=(1310.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1337.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(23dihydro14dioxin)"""),
)

species(
    label = 'OC(F)=COC(F)(F)C=C(F)F(9227)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,S} {15,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {12,D} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 C u0 p0 c0 {4,S} {5,S} {9,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1303.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44126,0.118518,-0.000148399,8.85871e-08,-2.04548e-11,-156527,34.3871], Tmin=(100,'K'), Tmax=(1065.81,'K')), NASAPolynomial(coeffs=[24.7649,0.0201651,-9.97822e-06,2.00449e-09,-1.45619e-13,-162114,-93.7128], Tmin=(1065.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1303.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF)"""),
)

species(
    label = '[O][C](F)C1OC(F)(F)C1C(F)F(9228)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u1 p2 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {6,S} {8,S} {12,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
12 C u1 p0 c0 {5,S} {7,S} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1022.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.233785,0.100073,-0.000125696,8.40789e-08,-2.2856e-11,-122860,34.8731], Tmin=(100,'K'), Tmax=(889.865,'K')), NASAPolynomial(coeffs=[13.9328,0.0363925,-1.83519e-05,3.65806e-09,-2.62224e-13,-125381,-31.8192], Tmin=(889.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1022.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + group(CsCsFFH) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'O=C=[C]OC(F)(F)[CH]C(F)F(9110)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u1 p0 c0 {5,S} {11,D}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-657.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.859857,0.11469,-0.000177803,1.36989e-07,-4.14054e-11,-78956.2,35.807], Tmin=(100,'K'), Tmax=(814.183,'K')), NASAPolynomial(coeffs=[17.2455,0.0257372,-1.39186e-05,2.79431e-09,-1.98949e-13,-81904.4,-47.8187], Tmin=(814.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-657.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)=[C]OC(F)(F)[CH]C(F)F(9229)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u0 p2 c0 {11,S} {15,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {5,S} {7,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-997.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,293,496,537,1218,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53889,0.131313,-0.000199795,1.52146e-07,-4.56764e-11,-119812,39.4627], Tmin=(100,'K'), Tmax=(817.889,'K')), NASAPolynomial(coeffs=[18.7877,0.0319045,-1.74822e-05,3.54286e-09,-2.54111e-13,-123137,-54.5152], Tmin=(817.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-997.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)=[C]OC(F)(F)CC(F)F(9230)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u1 p2 c0 {11,S}
8  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {15,S}
11 C u0 p0 c0 {5,S} {7,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1066.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,293,496,537,1218,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.961672,0.117882,-0.000162545,1.14382e-07,-3.21623e-11,-128044,37.2827], Tmin=(100,'K'), Tmax=(866.741,'K')), NASAPolynomial(coeffs=[17.0341,0.0348332,-1.88222e-05,3.83729e-09,-2.77788e-13,-131164,-46.9633], Tmin=(866.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1066.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)=COC(F)(F)[CH][C](F)F(9231)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,S} {15,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 C u1 p0 c0 {4,S} {5,S} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1036.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,326,540,652,719,1357,190,488,555,1236,1407,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59806,0.131112,-0.000190837,1.36257e-07,-3.81297e-11,-124524,37.9105], Tmin=(100,'K'), Tmax=(877.812,'K')), NASAPolynomial(coeffs=[20.904,0.0285737,-1.56184e-05,3.18309e-09,-2.29928e-13,-128474,-67.7163], Tmin=(877.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1036.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FO[C](F)[CH]OC(F)=CC(F)F(9232)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u0 p2 c0 {5,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {13,S}
9  C u0 p0 c0 {8,S} {10,D} {14,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 C u1 p0 c0 {6,S} {12,S} {15,S}
12 C u1 p0 c0 {3,S} {7,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-663.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3025,407.5,1350,352.5,395,473,707,1436,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03857,0.143913,-0.000243295,2.05535e-07,-6.78929e-11,-79563.8,39.6582], Tmin=(100,'K'), Tmax=(791.081,'K')), NASAPolynomial(coeffs=[18.4658,0.0326677,-1.80102e-05,3.58955e-09,-2.52218e-13,-82571.2,-52.9619], Tmin=(791.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-663.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFFH) + group(Cds-CdsCsH) + group(CdCFO) + radical(CCsJOC(O)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[CH][CH]C(F)(F)OC=C(F)OF(9233)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 C u1 p0 c0 {4,S} {9,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-672.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,326,540,652,719,1357,334,575,1197,1424,3202,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68954,0.134933,-0.000210687,1.65478e-07,-5.1221e-11,-80627.7,40.15], Tmin=(100,'K'), Tmax=(794.072,'K')), NASAPolynomial(coeffs=[18.5795,0.0328376,-1.78427e-05,3.58563e-09,-2.55284e-13,-83847,-52.964], Tmin=(794.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-672.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'FO[C]=COC(F)(F)[CH]C(F)F(9234)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {5,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {6,S} {12,D} {15,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-681.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57022,0.128273,-0.000181226,1.246e-07,-3.34017e-11,-81775.2,40.0237], Tmin=(100,'K'), Tmax=(917.698,'K')), NASAPolynomial(coeffs=[21.9818,0.0256162,-1.34304e-05,2.70361e-09,-1.94488e-13,-86097.9,-71.5781], Tmin=(917.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-681.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(H)) + radical(C=CJO)"""),
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
    E0 = (-355.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (37.0395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (63.5477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-347.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-316.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-171.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-329.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-331.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-271.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-64.3082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-264.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-90.3145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-205.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-41.8929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-25.9742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-227.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-275.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-151.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-118.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-159.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-158.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-147.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-165.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-161.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (209.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-347.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-330.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-279.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (94.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-104.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-238.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-213.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (127.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (120.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (113.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['O=CC(=O)F(335)', 'FC(F)=CC(F)F(2834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=O)F(8932)', '[O]C(F)(F)[CH]C(F)F(2594)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(F)F-2(3493)', 'O=C(F)[CH]O[C](F)F(1755)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['O=C(F)C1OC(F)(F)C1C(F)F(8981)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['O=C(F)COC(F)(F)C=C(F)F(8987)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['F[C]1OC1OC(F)(F)[CH]C(F)F(9213)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['F[C]1[CH]OC(F)(F)C(C(F)F)O1(9030)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['[O]C1(F)[CH]OC(F)(F)C1C(F)F(9214)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(13013.2,'s^-1'), n=1.81618, Ea=(24.1174,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_csHNd] for rate rule [R6;carbonylbond_intra;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(=O)F(335)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.06973e-08,'m^3/(mol*s)'), n=3.60774, Ea=(23.3001,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5239778944948545, var=2.087050032983542, Tref=1000.0, N=387, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O=C(F)[CH]OC(F)=CC(F)F(9215)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.2617,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]C(=O)F(509)', 'FC(F)=CC(F)F(2834)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'O=C(F)[CH]OC(F)(F)C=CF(9216)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(56.7354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=C(F)[CH]OC(F)(F)C=C(F)F(9217)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(3.29466,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C=COC(F)(F)[CH]C(F)F(9218)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.9385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(509)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['O=C(F)[CH]OC(F)(F)C[C](F)F(9219)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C(F)COC(F)(F)[CH][C](F)F(9220)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(61692.6,'s^-1'), n=1.86119, Ea=(60.1503,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_H/CO]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C(F)[CH]O[C](F)C(F)C(F)F(9221)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(162.072,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C(F)C(F)O[C](F)[CH]C(F)F(9222)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(219.645,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C(F)[CH]OC(F)(F)C(F)[CH]F(9223)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.234,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['O=C(F)C(F)OC(F)(F)[CH][CH]F(9224)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(197.338,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C]=COC(F)(F)C(F)C(F)F(9225)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(190.021,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]C(F)OC(F)(F)[CH]C(F)F(9099)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(C=O)C(F)(F)[CH]C(F)F(8975)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(132.515,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(6)', 'F[C]=COC(F)(F)[CH]C(F)F(9226)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['FC1=COC(F)(F)C(C(F)F)O1(9043)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['OC(F)=COC(F)(F)C=C(F)F(9227)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    products = ['[O][C](F)C1OC(F)(F)C1C(F)F(9228)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(75.5624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 71.5 to 75.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'O=C=[C]OC(F)(F)[CH]C(F)F(9110)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(290.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction30',
    reactants = ['OC(F)=[C]OC(F)(F)[CH]C(F)F(9229)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)=[C]OC(F)(F)CC(F)F(9230)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(9.51492e+06,'s^-1'), n=1.58625, Ea=(84.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS_OCs;Y_rad_out;Cs_H_out_H/NonDeC] + [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OC(F)=COC(F)(F)[CH][C](F)F(9231)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(72286.1,'s^-1'), n=1.57898, Ea=(80.9523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_single;XH_out] for rate rule [R7HJ_1;C_rad_out_noH;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FO[C](F)[CH]OC(F)=CC(F)F(9232)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(47.5589,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction34',
    reactants = ['F[CH][CH]C(F)(F)OC=C(F)OF(9233)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(49.5154,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction35',
    reactants = ['FO[C]=COC(F)(F)[CH]C(F)F(9234)'],
    products = ['O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(52.0458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #2384',
    isomers = [
        'O=C(F)[CH]OC(F)(F)[CH]C(F)F(8977)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'FC(F)=CC(F)F(2834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2384',
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

