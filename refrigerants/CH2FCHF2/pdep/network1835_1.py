species(
    label = 'O=C(F)[CH]OO[C](F)CF(6213)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-610.105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.389635,0.105786,-0.000169482,1.46908e-07,-5.11987e-11,-73229.4,33.593], Tmin=(100,'K'), Tmax=(754.463,'K')), NASAPolynomial(coeffs=[11.4958,0.0366597,-1.98958e-05,3.99074e-09,-2.83406e-13,-74848.9,-19.2456], Tmin=(754.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(CsCsF1sO2s) + radical(OCJC=O)"""),
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
    label = 'O=C(F)CF(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]O[CH]C(=O)F(1948)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u1 p0 c0 {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-215.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(0.2824,'amu*angstrom^2'), symmetry=1, barrier=(6.49293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4819,'amu*angstrom^2'), symmetry=1, barrier=(57.0637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06408,0.0468221,-7.59984e-05,6.61717e-08,-2.2668e-11,-25811.9,19.1781], Tmin=(100,'K'), Tmax=(821.51,'K')), NASAPolynomial(coeffs=[7.05431,0.0160915,-8.14163e-06,1.57328e-09,-1.08835e-13,-26414.8,-2.59461], Tmin=(821.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(ROOJ) + radical(OCJC=O)"""),
)

species(
    label = 'F[C]CF-2(3008)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u2 p0 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-90.9927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([90,150,250,247,323,377,431,433,1065,1274,4000],'cm^-1')),
        HinderedRotor(inertia=(1.0536e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30451,0.015538,-7.22211e-06,9.68417e-10,5.6631e-14,-10919.1,11.642], Tmin=(100,'K'), Tmax=(1819.26,'K')), NASAPolynomial(coeffs=[8.39313,0.00776393,-3.62739e-06,6.82728e-10,-4.58702e-14,-13335.6,-17.5057], Tmin=(1819.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]C(=O)F(323)',
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
    label = '[O]O[C](F)CF(188)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-233.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180],'cm^-1')),
        HinderedRotor(inertia=(0.397531,'amu*angstrom^2'), symmetry=1, barrier=(9.14001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396473,'amu*angstrom^2'), symmetry=1, barrier=(9.1157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7711,0.058036,-0.000110671,1.05783e-07,-3.76048e-11,-27970.8,20.3112], Tmin=(100,'K'), Tmax=(874.045,'K')), NASAPolynomial(coeffs=[5.51302,0.0201688,-1.00878e-05,1.91304e-09,-1.2907e-13,-27832.6,7.29488], Tmin=(874.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C(F)C1OOC1(F)CF(6216)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-845.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.373911,0.0703799,-6.1386e-05,2.5411e-08,-4.118e-12,-101515,30.9863], Tmin=(100,'K'), Tmax=(1484.44,'K')), NASAPolynomial(coeffs=[19.8008,0.0180306,-8.48673e-06,1.6532e-09,-1.16776e-13,-107282,-70.4113], Tmin=(1484.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-845.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s)"""),
)

species(
    label = 'O=C(F)COOC(F)=CF(6220)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-752.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79501,0.0835197,-9.01896e-05,9.4965e-09,4.04576e-11,-90394.2,30.3174], Tmin=(100,'K'), Tmax=(501.329,'K')), NASAPolynomial(coeffs=[9.18952,0.0384212,-2.07177e-05,4.16779e-09,-2.97193e-13,-91510.8,-7.12715], Tmin=(501.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-752.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FC[C](F)OOC1O[C]1F(6510)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {10,S} {12,S} {13,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u1 p0 c0 {2,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-418.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342056,0.0883063,-0.000108247,7.04843e-08,-1.89027e-11,-50154.4,30.1784], Tmin=(100,'K'), Tmax=(894.428,'K')), NASAPolynomial(coeffs=[12.2602,0.0350069,-1.88617e-05,3.86051e-09,-2.8081e-13,-52286.4,-25.99], Tmin=(894.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-418.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'FCC1(F)OO[CH][C](F)O1(6487)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {10,S}
10 C u1 p0 c0 {6,S} {9,S} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-587.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.927811,0.0933359,-0.000107637,6.07626e-08,-1.27676e-11,-70482.4,28.3376], Tmin=(100,'K'), Tmax=(1352.75,'K')), NASAPolynomial(coeffs=[21.7082,0.0114416,-2.38843e-07,-3.41724e-10,3.59929e-14,-75237.7,-82.6475], Tmin=(1352.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(124trioxane) + radical(CsCsF1sO2s) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1(F)CF(6511)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
10 C u1 p0 c0 {5,S} {8,S} {13,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-557.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.36777,0.0869262,-9.71145e-05,5.44164e-08,-1.17333e-11,-66834.9,28.1573], Tmin=(100,'K'), Tmax=(1197.04,'K')), NASAPolynomial(coeffs=[19.8365,0.0162278,-4.53255e-06,6.32425e-10,-3.64419e-14,-71443.9,-71.9973], Tmin=(1197.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-557.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(12dioxolane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[O][C](F)CF(882)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-241.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.77787,'amu*angstrom^2'), symmetry=1, barrier=(17.8848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43603,0.0332746,-3.37965e-05,1.70828e-08,-3.40461e-12,-28931.2,16.3061], Tmin=(100,'K'), Tmax=(1215.45,'K')), NASAPolynomial(coeffs=[9.73799,0.00924406,-4.14003e-06,8.16317e-10,-5.88299e-14,-30706.2,-20.3463], Tmin=(1215.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'C=C(F)OO[CH]C(=O)F(6512)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {3,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-447.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,326,540,652,719,1357,611,648,830,1210,1753,2950,3100,1380,975,1025,1650,180,2064.18],'cm^-1')),
        HinderedRotor(inertia=(0.289696,'amu*angstrom^2'), symmetry=1, barrier=(6.66068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95574,'amu*angstrom^2'), symmetry=1, barrier=(44.9663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95389,'amu*angstrom^2'), symmetry=1, barrier=(44.9237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95439,'amu*angstrom^2'), symmetry=1, barrier=(44.9353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28137,0.072054,-6.50193e-05,-2.52121e-08,6.3264e-11,-53685.4,27.5796], Tmin=(100,'K'), Tmax=(478.989,'K')), NASAPolynomial(coeffs=[8.11853,0.0362143,-1.93514e-05,3.87686e-09,-2.75741e-13,-54584.2,-2.91846], Tmin=(478.989,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-447.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + group(COCsFO) + group(Cds-CdsHH) + radical(OCJC=O)"""),
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
    label = 'O=C(F)[CH]OOC(F)=CF(6513)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {4,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {10,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-613.945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,326,540,652,719,1357,611,648,830,1210,1753,194,682,905,1196,1383,3221,180,940.818],'cm^-1')),
        HinderedRotor(inertia=(0.225696,'amu*angstrom^2'), symmetry=1, barrier=(5.18919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226559,'amu*angstrom^2'), symmetry=1, barrier=(5.20905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25827,'amu*angstrom^2'), symmetry=1, barrier=(51.9221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26194,'amu*angstrom^2'), symmetry=1, barrier=(52.0066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0200604,0.0982857,-0.000167402,1.54057e-07,-5.58458e-11,-73707.7,32.6992], Tmin=(100,'K'), Tmax=(795.903,'K')), NASAPolynomial(coeffs=[8.74215,0.0373996,-2.05176e-05,4.10538e-09,-2.89918e-13,-74556,-3.99605], Tmin=(795.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COCsFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(OCJC=O)"""),
)

species(
    label = 'O=C=COO[C](F)CF(6514)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {3,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-312.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06561,'amu*angstrom^2'), symmetry=1, barrier=(24.5004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06088,'amu*angstrom^2'), symmetry=1, barrier=(24.3918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0688,'amu*angstrom^2'), symmetry=1, barrier=(24.5739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06034,'amu*angstrom^2'), symmetry=1, barrier=(24.3793,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198364,0.0902907,-0.000132845,1.02114e-07,-3.13641e-11,-37443.1,29.9507], Tmin=(100,'K'), Tmax=(796.519,'K')), NASAPolynomial(coeffs=[12.7379,0.027314,-1.42381e-05,2.83508e-09,-2.01575e-13,-39440.5,-27.6916], Tmin=(796.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[C]CF(126)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p1 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-105.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,734,1109,1255,1358,2983,3011,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.0460479,'amu*angstrom^2'), symmetry=1, barrier=(35.5232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4418.31,'J/mol'), sigma=(4.687e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8659,0.0107238,1.79401e-05,-3.81644e-08,1.95617e-11,-12729.5,8.48308], Tmin=(10,'K'), Tmax=(672.698,'K')), NASAPolynomial(coeffs=[3.4045,0.0188633,-1.22417e-05,3.67097e-09,-4.17395e-13,-12789.5,9.6187], Tmin=(672.698,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-105.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[CH]OOC(F)[CH]F(6515)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {5,S} {10,S} {12,S}
9  C u1 p0 c0 {2,S} {7,S} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-606.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,334,575,1197,1424,3202,611,648,830,1210,1753,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.264948,0.102449,-0.000148304,1.13398e-07,-3.50157e-11,-72741.4,32.716], Tmin=(100,'K'), Tmax=(788.844,'K')), NASAPolynomial(coeffs=[13.1869,0.0342361,-1.85904e-05,3.77037e-09,-2.71227e-13,-74863.6,-28.9901], Tmin=(788.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(OCJC=O)"""),
)

species(
    label = 'O=C(F)COO[C](F)[CH]F(6516)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u1 p0 c0 {3,S} {8,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-549.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,486,617,768,1157,1926,334,575,1197,1424,3202,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.598434,0.110971,-0.000184486,1.60955e-07,-5.55844e-11,-65987.4,34.936], Tmin=(100,'K'), Tmax=(783.247,'K')), NASAPolynomial(coeffs=[12.4881,0.0342651,-1.86765e-05,3.72991e-09,-2.63299e-13,-67734.5,-23.0682], Tmin=(783.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-549.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)(F)OO[CH]C(=O)F(6517)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
8  C u1 p0 c0 {5,S} {10,S} {11,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-651.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.185174,0.100434,-0.000142748,1.0655e-07,-3.20942e-11,-78182.9,32.9743], Tmin=(100,'K'), Tmax=(808.093,'K')), NASAPolynomial(coeffs=[13.3874,0.0332505,-1.80385e-05,3.66502e-09,-2.64202e-13,-80376.4,-29.6134], Tmin=(808.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-651.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFFO) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(CJCOOH) + radical(OCJC=O)"""),
)

species(
    label = '[CH2][C](F)OOC(F)C(=O)F(6518)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-589.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,395,473,707,1436,486,617,768,1157,1926,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.20335,'amu*angstrom^2'), symmetry=1, barrier=(4.67541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73229,'amu*angstrom^2'), symmetry=1, barrier=(39.8287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73217,'amu*angstrom^2'), symmetry=1, barrier=(39.826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547544,'amu*angstrom^2'), symmetry=1, barrier=(12.5891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73232,'amu*angstrom^2'), symmetry=1, barrier=(39.8294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.477161,0.107863,-0.000175434,1.52048e-07,-5.27355e-11,-70701.7,34.9603], Tmin=(100,'K'), Tmax=(756.044,'K')), NASAPolynomial(coeffs=[12.1341,0.0351501,-1.92868e-05,3.8803e-09,-2.75794e-13,-72437.4,-21.2225], Tmin=(756.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-589.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(COCsFO) + radical(CsCsF1sO2s) + radical(CJCOOH)"""),
)

species(
    label = 'O=[C][CH]OOC(F)(F)CF(6519)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-625.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0342698,0.0970083,-0.000132265,9.52576e-08,-2.78597e-11,-75090.3,31.6002], Tmin=(100,'K'), Tmax=(829.102,'K')), NASAPolynomial(coeffs=[13.003,0.0341116,-1.84758e-05,3.76389e-09,-2.72217e-13,-77252.1,-28.8541], Tmin=(829.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]C(F)OO[C](F)CF(6352)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {4,S} {10,S} {13,S}
9  C u1 p0 c0 {3,S} {5,S} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-564.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,355,410,600,1181,1341,1420,3056,395,473,707,1436,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.71679,'amu*angstrom^2'), symmetry=1, barrier=(39.4723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71529,'amu*angstrom^2'), symmetry=1, barrier=(39.438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830144,'amu*angstrom^2'), symmetry=1, barrier=(19.0866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.829818,'amu*angstrom^2'), symmetry=1, barrier=(19.0791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830797,'amu*angstrom^2'), symmetry=1, barrier=(19.1017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.743092,0.11585,-0.000200844,1.78906e-07,-6.19432e-11,-67758.6,33.9776], Tmin=(100,'K'), Tmax=(819.95,'K')), NASAPolynomial(coeffs=[12.2022,0.0344802,-1.86616e-05,3.68408e-09,-2.57043e-13,-69269.1,-22.1719], Tmin=(819.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(CsCsF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(F)(C=O)O[C](F)CF(6211)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-738.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.71017,0.115723,-0.000204607,1.82786e-07,-6.27203e-11,-88641.7,33.2009], Tmin=(100,'K'), Tmax=(842.947,'K')), NASAPolynomial(coeffs=[12.1091,0.0326716,-1.72773e-05,3.35897e-09,-2.31377e-13,-90013.4,-21.7715], Tmin=(842.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs))"""),
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
    label = 'F[C]=COO[C](F)CF(6520)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-191.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.615656,'amu*angstrom^2'), symmetry=1, barrier=(14.1551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72887,'amu*angstrom^2'), symmetry=1, barrier=(39.7502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72784,'amu*angstrom^2'), symmetry=1, barrier=(39.7264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615448,'amu*angstrom^2'), symmetry=1, barrier=(14.1504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0428325,0.0996848,-0.000172631,1.57095e-07,-5.57162e-11,-22851.8,31.9603], Tmin=(100,'K'), Tmax=(814.223,'K')), NASAPolynomial(coeffs=[9.65694,0.0340088,-1.84344e-05,3.65446e-09,-2.56006e-13,-23833.9,-9.17301], Tmin=(814.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCsF1sO2s) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FCC1(F)OOC=C(F)O1(6501)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-784.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0762432,0.0795957,-7.98058e-05,3.92144e-08,-7.48125e-12,-94215.6,25.0815], Tmin=(100,'K'), Tmax=(1285.69,'K')), NASAPolynomial(coeffs=[20.1527,0.0166601,-6.37996e-06,1.1412e-09,-7.80163e-14,-99417.2,-77.5952], Tmin=(1285.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-784.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFOO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
)

species(
    label = 'FC1=COO1(285)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-52.3112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,353.677,1138.56,1138.58,1138.63,1726.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95335,0.00258575,6.07332e-05,-1.08656e-07,5.80728e-11,-6290.59,8.5661], Tmin=(10,'K'), Tmax=(614.296,'K')), NASAPolynomial(coeffs=[2.56312,0.0230016,-1.68658e-05,5.67097e-09,-7.10066e-13,-6334.19,12.8506], Tmin=(614.296,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-52.3112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC1DCOO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)=COOC(F)=CF(6521)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u0 p0 c0 {2,S} {6,S} {7,D}
10 C u0 p0 c0 {3,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-656.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22587,0.101483,-0.000161317,1.3789e-07,-4.73622e-11,-78770.8,33.2537], Tmin=(100,'K'), Tmax=(748.709,'K')), NASAPolynomial(coeffs=[11.7116,0.033882,-1.82174e-05,3.64776e-09,-2.58943e-13,-80451.2,-20.1667], Tmin=(748.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = '[O][C](F)C1OOC1(F)CF(6522)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
10 C u1 p0 c0 {3,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-475.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267778,0.0889683,-0.000123633,9.29721e-08,-2.83617e-11,-57049.1,30.5346], Tmin=(100,'K'), Tmax=(797.465,'K')), NASAPolynomial(coeffs=[11.6447,0.0319028,-1.62945e-05,3.23873e-09,-2.30719e-13,-58863.6,-21.7775], Tmin=(797.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)-O2s-O2s) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'O=C=[C]OO[C](F)CF(6361)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {3,S} {6,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-72.6628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.83744,'amu*angstrom^2'), symmetry=1, barrier=(42.2464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459004,'amu*angstrom^2'), symmetry=1, barrier=(10.5534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458586,'amu*angstrom^2'), symmetry=1, barrier=(10.5438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47785,'amu*angstrom^2'), symmetry=1, barrier=(33.9787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201999,0.0947319,-0.000170382,1.56725e-07,-5.53293e-11,-8613.34,33.0501], Tmin=(100,'K'), Tmax=(835.662,'K')), NASAPolynomial(coeffs=[9.15298,0.030441,-1.64856e-05,3.24104e-09,-2.24946e-13,-9360.52,-4.0459], Tmin=(835.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.6628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCsF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)=[C]OO[C](F)CF(6523)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-412.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,293,496,537,1218,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.565215,0.112513,-0.000197035,1.78955e-07,-6.31566e-11,-49465.7,37.0148], Tmin=(100,'K'), Tmax=(816.445,'K')), NASAPolynomial(coeffs=[10.8605,0.0363082,-1.98682e-05,3.94541e-09,-2.76354e-13,-50657.2,-11.6622], Tmin=(816.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-412.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCsF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)[C]OOC(F)CF(6524)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {6,D} {10,S}
10 C u2 p0 c0 {5,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-482.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,486,617,768,1157,1926,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0245853,0.0949005,-0.000117191,7.21263e-08,-1.77866e-11,-57850.2,30.4783], Tmin=(100,'K'), Tmax=(980.486,'K')), NASAPolynomial(coeffs=[16.1881,0.0287591,-1.60049e-05,3.32607e-09,-2.44222e-13,-61029.4,-47.4191], Tmin=(980.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(CH2_triplet)"""),
)

species(
    label = 'OC(F)=COO[C](F)[CH]F(6525)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {7,D}
10 C u1 p0 c0 {3,S} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-453.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3010,987.5,1337.5,450,1655,395,473,707,1436,326,540,652,719,1357,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.674996,'amu*angstrom^2'), symmetry=1, barrier=(15.5195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.67495,'amu*angstrom^2'), symmetry=1, barrier=(15.5184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68612,'amu*angstrom^2'), symmetry=1, barrier=(38.7672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674925,'amu*angstrom^2'), symmetry=1, barrier=(15.5179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68564,'amu*angstrom^2'), symmetry=1, barrier=(38.7562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.9351,0.118004,-0.000198065,1.68696e-07,-5.63669e-11,-54393.3,35.569], Tmin=(100,'K'), Tmax=(792.973,'K')), NASAPolynomial(coeffs=[14.9799,0.0298011,-1.62306e-05,3.22474e-09,-2.26415e-13,-56668.2,-35.9494], Tmin=(792.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[CH2][C](F)OOC=C(F)OF(6526)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {7,D}
10 C u1 p0 c0 {8,S} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-106.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,3010,987.5,1337.5,450,1655,395,473,707,1436,326,540,652,719,1357,3000,3100,440,815,1455,1000,180,2929.24],'cm^-1')),
        HinderedRotor(inertia=(0.315806,'amu*angstrom^2'), symmetry=1, barrier=(7.26101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73749,'amu*angstrom^2'), symmetry=1, barrier=(39.9484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317124,'amu*angstrom^2'), symmetry=1, barrier=(7.2913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73462,'amu*angstrom^2'), symmetry=1, barrier=(39.8824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73391,'amu*angstrom^2'), symmetry=1, barrier=(39.8659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.591212,0.115598,-0.000210453,1.97767e-07,-7.14091e-11,-12641.3,37.1316], Tmin=(100,'K'), Tmax=(825.322,'K')), NASAPolynomial(coeffs=[9.12039,0.0403282,-2.23962e-05,4.45769e-09,-3.12003e-13,-13283.9,-2.03794], Tmin=(825.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCsF1sO2s) + radical(CJCOOH)"""),
)

species(
    label = 'FC[C](F)OOC=[C]OF(6527)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {3,S} {10,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-96.3324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3010,987.5,1337.5,450,1655,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.242482,0.104209,-0.000153536,1.05867e-07,-2.28482e-11,-11443.7,36.3671], Tmin=(100,'K'), Tmax=(609.206,'K')), NASAPolynomial(coeffs=[13.615,0.0308541,-1.63355e-05,3.2357e-09,-2.2792e-13,-13459.3,-26.3048], Tmin=(609.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.3324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CsCsF1sO2s) + radical(C=CJO)"""),
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
    E0 = (-246.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (57.8207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (125.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-237.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-206.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-54.0753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-175.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-184.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-246.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-231.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (39.0747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-29.5131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (146.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-73.6724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (45.4288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-97.7424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-125.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-112.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-37.3582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-55.7478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-34.8484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-74.2792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (415.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-238.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (71.7726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-221.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-111.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (260.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (101.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-73.8749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-8.73011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (292.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (301.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['O=CC(=O)F(335)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]O[CH]C(=O)F(1948)', 'F[C]CF-2(3008)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(=O)F(323)', '[O]O[C](F)CF(188)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['O=C(F)C1OOC1(F)CF(6216)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['FC[C](F)OOC1O[C]1F(6510)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(192.058,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 190.6 to 192.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['FCC1(F)OO[CH][C](F)O1(6487)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['[O]C1(F)[CH]OOC1(F)CF(6511)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(61.9573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]C(=O)F(509)', 'O=C(F)CF(879)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(216.316,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Ea raised from 215.8 to 216.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=CC(=O)F(335)', '[O][C](F)CF(882)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(129.146,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C=C(F)OO[CH]C(=O)F(6512)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(49.2862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C(F)[CH]OOC(F)=CF(6513)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(50675.7,'m^3/(mol*s)'), n=0.856706, Ea=(8.65598,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3149223058909152, var=1.98898808922157, Tref=1000.0, N=74, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=C=COO[C](F)CF(6514)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(21.7215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]C(=O)F(509)', '[O][C](F)CF(882)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(17.8562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]O[CH]C(=O)F(1948)', 'F[C]CF(126)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(2.46435,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C(F)[CH]OOC(F)[CH]F(6515)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(34257.4,'s^-1'), n=2.45724, Ea=(144.304,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out] for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C(F)COO[C](F)[CH]F(6516)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(61692.6,'s^-1'), n=1.86119, Ea=(60.1503,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;C_rad_out_1H;Cs_H_out_H/CO]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['[CH2]C(F)(F)OO[CH]C(=O)F(6517)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(133.813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C](F)OOC(F)C(=O)F(6518)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(187.783,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['O=[C][CH]OOC(F)(F)CF(6519)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(190.385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=[C]C(F)OO[C](F)CF(6352)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(165.885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['[O]C(F)(C=O)O[C](F)CF(6211)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(171.854,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(6)', 'F[C]=COO[C](F)CF(6520)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['FCC1(F)OOC=C(F)O1(6501)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['FC1=COO1(285)', '[O][C](F)CF(882)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(317.906,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['OC(F)=COOC(F)=CF(6521)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['[O][C](F)C1OOC1(F)CF(6522)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(134.706,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 132.5 to 134.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', 'O=C=[C]OO[C](F)CF(6361)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(249.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OC(F)=[C]OO[C](F)CF(6523)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=C(F)[C]OOC(F)CF(6524)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['OC(F)=COO[C](F)[CH]F(6525)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(72286.1,'s^-1'), n=1.57898, Ea=(80.9523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_single;XH_out] for rate rule [R7HJ_1;C_rad_out_1H;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C](F)OOC=C(F)OF(6526)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(34.9661,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FC[C](F)OOC=[C]OF(6527)'],
    products = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(33.7202,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1835',
    isomers = [
        'O=C(F)[CH]OO[C](F)CF(6213)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1835',
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

