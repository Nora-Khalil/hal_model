species(
    label = 'F[CH]C(F)C(F)(F)[CH]C(F)F(6699)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u1 p0 c0 {6,S} {8,S} {15,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-989.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,334,575,1197,1424,3202,266.207,266.207,266.207,266.207],'cm^-1')),
        HinderedRotor(inertia=(0.142114,'amu*angstrom^2'), symmetry=1, barrier=(7.14669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142115,'amu*angstrom^2'), symmetry=1, barrier=(7.14668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.676711,'amu*angstrom^2'), symmetry=1, barrier=(34.0306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142114,'amu*angstrom^2'), symmetry=1, barrier=(7.14669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10795,0.122673,-0.000199594,1.71717e-07,-5.86065e-11,-118819,36.8582], Tmin=(100,'K'), Tmax=(784.646,'K')), NASAPolynomial(coeffs=[13.5826,0.0383295,-2.02858e-05,4.01606e-09,-2.8241e-13,-120833,-28.5983], Tmin=(784.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-989.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'CHFCHF[Z](59)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
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
    label = 'F[CH]C(F)C([C](F)F)C(F)F(6700)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {14,S}
10 C u1 p0 c0 {4,S} {5,S} {7,S}
11 C u1 p0 c0 {6,S} {8,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-993.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,334,575,1197,1424,3202,327.929,327.929,327.929,327.929],'cm^-1')),
        HinderedRotor(inertia=(0.124288,'amu*angstrom^2'), symmetry=1, barrier=(9.48451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124288,'amu*angstrom^2'), symmetry=1, barrier=(9.4845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124288,'amu*angstrom^2'), symmetry=1, barrier=(9.4845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453634,'amu*angstrom^2'), symmetry=1, barrier=(34.6167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3132.16,'J/mol'), sigma=(5.21338,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.24 K, Pc=50.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08835,0.121604,-0.000195107,1.65339e-07,-5.56791e-11,-119257,36.4754], Tmin=(100,'K'), Tmax=(779.739,'K')), NASAPolynomial(coeffs=[14.0287,0.0371742,-1.94498e-05,3.83704e-09,-2.69487e-13,-121406,-31.3531], Tmin=(779.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-993.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)[C](F)F(818)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 C u1 p0 c0 {4,S} {5,S} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-506.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.165161,'amu*angstrom^2'), symmetry=1, barrier=(3.79738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165272,'amu*angstrom^2'), symmetry=1, barrier=(3.79994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23835,0.0729054,-0.000142741,1.40496e-07,-5.15307e-11,-60798.1,23.0029], Tmin=(100,'K'), Tmax=(853.903,'K')), NASAPolynomial(coeffs=[5.03358,0.0276124,-1.48436e-05,2.90815e-09,-2.00625e-13,-60443.2,11.1663], Tmin=(853.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = '[CH]F(804)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]C(F)F(1054)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u1 p0 c0 {5,S} {6,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-764.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,334,575,1197,1424,3202,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.242395,'amu*angstrom^2'), symmetry=1, barrier=(5.57315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242202,'amu*angstrom^2'), symmetry=1, barrier=(5.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242309,'amu*angstrom^2'), symmetry=1, barrier=(5.57116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00640479,0.101749,-0.00018755,1.77353e-07,-6.38397e-11,-91814.8,29.977], Tmin=(100,'K'), Tmax=(838.537,'K')), NASAPolynomial(coeffs=[8.04264,0.0354661,-1.90946e-05,3.764e-09,-2.61651e-13,-92184.2,-1.59146], Tmin=(838.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-764.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FC(F)C1C(F)C(F)C1(F)F(6963)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
10 C u0 p0 c0 {4,S} {8,S} {9,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1245.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.034012,0.0830291,-7.55483e-05,3.41326e-08,-6.16602e-12,-149626,29.2454], Tmin=(100,'K'), Tmax=(1320.47,'K')), NASAPolynomial(coeffs=[18.2433,0.0278691,-1.28888e-05,2.4977e-09,-1.76704e-13,-154435,-63.6661], Tmin=(1320.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1245.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsFFH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FC=C(F)C(F)(F)CC(F)F(7058)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {14,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1236.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.426979,0.10914,-0.000151692,1.00666e-07,-1.98506e-11,-148623,30.6742], Tmin=(100,'K'), Tmax=(598.547,'K')), NASAPolynomial(coeffs=[12.7668,0.0388152,-2.01789e-05,4.00218e-09,-2.8326e-13,-150522,-28.8768], Tmin=(598.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1236.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FCC(F)C(F)(F)C=C(F)F(6701)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {11,D} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1227.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.264405,0.101788,-0.000133186,9.36827e-08,-2.68609e-11,-147459,30.9601], Tmin=(100,'K'), Tmax=(844.348,'K')), NASAPolynomial(coeffs=[13.2474,0.0377801,-1.94787e-05,3.90701e-09,-2.8057e-13,-149741,-31.9408], Tmin=(844.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1227.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF)"""),
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
    label = 'F[CH]C(F)C(F)=CC(F)F(7059)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 C u1 p0 c0 {5,S} {6,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-824.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.200998,'amu*angstrom^2'), symmetry=1, barrier=(4.62135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200933,'amu*angstrom^2'), symmetry=1, barrier=(4.61984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00781,'amu*angstrom^2'), symmetry=1, barrier=(23.1716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.181842,0.0999583,-0.000152021,1.27674e-07,-4.33804e-11,-99081.3,32.021], Tmin=(100,'K'), Tmax=(754.683,'K')), NASAPolynomial(coeffs=[11.1532,0.0363647,-1.86365e-05,3.67332e-09,-2.58906e-13,-100692,-18.8105], Tmin=(754.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-824.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]F(811)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-71.0739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,406,528,622,1148,1246,1368,1480,3164,3240,1663.86],'cm^-1')),
        HinderedRotor(inertia=(0.36711,'amu*angstrom^2'), symmetry=1, barrier=(8.44058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88781,0.0292781,-5.20579e-05,5.24786e-08,-2.00705e-11,-8512.84,14.1702], Tmin=(100,'K'), Tmax=(832.942,'K')), NASAPolynomial(coeffs=[3.49542,0.0155533,-7.88019e-06,1.5434e-09,-1.07577e-13,-8239.17,13.6003], Tmin=(832.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.0739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sH)"""),
)

species(
    label = 'FC=CC(F)(F)[CH]C(F)F(7060)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {5,S} {9,D} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-875.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.27528,'amu*angstrom^2'), symmetry=1, barrier=(6.32923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275401,'amu*angstrom^2'), symmetry=1, barrier=(6.33202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275386,'amu*angstrom^2'), symmetry=1, barrier=(6.33166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.552536,0.088226,-0.000102994,3.59636e-08,1.89077e-11,-105183,29.7234], Tmin=(100,'K'), Tmax=(526.72,'K')), NASAPolynomial(coeffs=[9.63449,0.0393049,-2.07708e-05,4.15907e-09,-2.96252e-13,-106418,-10.9094], Tmin=(526.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-875.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S)"""),
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
    label = 'FC=C(F)C(F)(F)[CH]C(F)F(7032)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1042.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.25879,'amu*angstrom^2'), symmetry=1, barrier=(5.95008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258728,'amu*angstrom^2'), symmetry=1, barrier=(5.94866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258848,'amu*angstrom^2'), symmetry=1, barrier=(5.95142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.682193,0.114403,-0.000193158,1.73109e-07,-6.0967e-11,-125216,34.0774], Tmin=(100,'K'), Tmax=(803.365,'K')), NASAPolynomial(coeffs=[11.1016,0.0389392,-2.09038e-05,4.147e-09,-2.91263e-13,-126567,-16.82], Tmin=(803.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1042.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_S)"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C=CF(7061)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {10,D} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {13,S}
10 C u0 p0 c0 {5,S} {8,D} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-833.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,194,682,905,1196,1383,3221,205.155,205.864,210.163],'cm^-1')),
        HinderedRotor(inertia=(0.771949,'amu*angstrom^2'), symmetry=1, barrier=(25.4063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202118,'amu*angstrom^2'), symmetry=1, barrier=(6.36773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200504,'amu*angstrom^2'), symmetry=1, barrier=(6.39298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.311749,0.102908,-0.000153514,1.23343e-07,-3.99532e-11,-100156,30.6341], Tmin=(100,'K'), Tmax=(754.235,'K')), NASAPolynomial(coeffs=[12.5325,0.0347859,-1.80272e-05,3.57968e-09,-2.53848e-13,-102093,-27.7085], Tmin=(754.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-833.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C=C(F)F(7062)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.197426,'amu*angstrom^2'), symmetry=1, barrier=(4.53922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196654,'amu*angstrom^2'), symmetry=1, barrier=(4.52146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14824,'amu*angstrom^2'), symmetry=1, barrier=(26.4003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.419719,0.108593,-0.000157594,1.10709e-07,-2.5783e-11,-124334,32.073], Tmin=(100,'K'), Tmax=(609.51,'K')), NASAPolynomial(coeffs=[13.2279,0.0356428,-1.89495e-05,3.78311e-09,-2.6843e-13,-126306,-29.5433], Tmin=(609.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    label = 'F[CH]C=C(F)[CH]C(F)F(7063)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u1 p0 c0 {5,S} {7,S} {11,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u1 p0 c0 {4,S} {8,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-574.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,271,519,563,612,1379,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,215.116,215.156],'cm^-1')),
        HinderedRotor(inertia=(0.287025,'amu*angstrom^2'), symmetry=1, barrier=(9.4257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0140665,'amu*angstrom^2'), symmetry=1, barrier=(35.6293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08475,'amu*angstrom^2'), symmetry=1, barrier=(35.6293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867279,0.0736678,-7.76319e-05,4.3325e-08,-1.00038e-11,-68989.3,24.9782], Tmin=(100,'K'), Tmax=(1027.35,'K')), NASAPolynomial(coeffs=[11.7583,0.0312633,-1.57181e-05,3.14766e-09,-2.26814e-13,-71227.1,-27.8583], Tmin=(1027.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-574.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Csj(Cd-CdH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)=C(F)[CH]C(F)F(6936)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {6,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 C u1 p0 c0 {5,S} {9,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-734.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,206,336,431,607,515,611,528,696,1312,1446,234,589,736,816,1240,3237,180,2453.42],'cm^-1')),
        HinderedRotor(inertia=(1.92252,'amu*angstrom^2'), symmetry=1, barrier=(44.2025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305872,'amu*angstrom^2'), symmetry=1, barrier=(7.03259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9252,'amu*angstrom^2'), symmetry=1, barrier=(44.2641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156496,0.0913922,-0.000122724,8.76715e-08,-2.5356e-11,-88235.5,27.8256], Tmin=(100,'K'), Tmax=(839.854,'K')), NASAPolynomial(coeffs=[12.6541,0.0318685,-1.64113e-05,3.28046e-09,-2.34874e-13,-90334.7,-30.2867], Tmin=(839.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-734.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Allyl_S) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'FC[C](F)C(F)(F)[CH]C(F)F(6851)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {10,S} {13,S} {14,S}
10 C u1 p0 c0 {6,S} {7,S} {9,S}
11 C u1 p0 c0 {7,S} {8,S} {15,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-990.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,522,611,926,1093,1137,1374,1416,3112,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,3025,407.5,1350,352.5,228.622,228.624,228.697,2231.53],'cm^-1')),
        HinderedRotor(inertia=(0.783163,'amu*angstrom^2'), symmetry=1, barrier=(29.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232709,'amu*angstrom^2'), symmetry=1, barrier=(8.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232787,'amu*angstrom^2'), symmetry=1, barrier=(8.62543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782538,'amu*angstrom^2'), symmetry=1, barrier=(29.0206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.707247,0.113312,-0.000177303,1.52664e-07,-5.3183e-11,-118992,37.1055], Tmin=(100,'K'), Tmax=(746.537,'K')), NASAPolynomial(coeffs=[11.7081,0.0414232,-2.20773e-05,4.41666e-09,-3.13727e-13,-120697,-18.1607], Tmin=(746.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-990.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C[C](F)F(6827)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
10 C u1 p0 c0 {4,S} {9,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-996.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,190,488,555,1236,1407,266.97,266.979,266.979,266.983],'cm^-1')),
        HinderedRotor(inertia=(0.144175,'amu*angstrom^2'), symmetry=1, barrier=(7.29181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144162,'amu*angstrom^2'), symmetry=1, barrier=(7.29144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.592431,'amu*angstrom^2'), symmetry=1, barrier=(29.9647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144162,'amu*angstrom^2'), symmetry=1, barrier=(7.2916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31044,0.129892,-0.000223983,1.99977e-07,-6.94096e-11,-119664,35.8475], Tmin=(100,'K'), Tmax=(821.965,'K')), NASAPolynomial(coeffs=[12.7548,0.0401736,-2.14392e-05,4.21764e-09,-2.9384e-13,-121258,-24.8806], Tmin=(821.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-996.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][C](F)C(F)(F)CC(F)F(7064)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {14,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {6,S} {10,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1006.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1709,0.125805,-0.000211951,1.88271e-07,-6.57891e-11,-120824,36.1953], Tmin=(100,'K'), Tmax=(799.998,'K')), NASAPolynomial(coeffs=[12.5166,0.040951,-2.20681e-05,4.38384e-09,-3.0815e-13,-122489,-23.5014], Tmin=(799.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1006.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'FCC(F)C(F)(F)[CH][C](F)F(7065)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
10 C u1 p0 c0 {8,S} {11,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-981.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,215,315,519,588,595,1205,1248,528,1116,1182,1331,1402,1494,3075,3110,3025,407.5,1350,352.5,190,488,555,1236,1407,299.512,299.518,299.528,1260.38],'cm^-1')),
        HinderedRotor(inertia=(0.130754,'amu*angstrom^2'), symmetry=1, barrier=(8.32429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130758,'amu*angstrom^2'), symmetry=1, barrier=(8.32426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130762,'amu*angstrom^2'), symmetry=1, barrier=(8.32426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.578008,'amu*angstrom^2'), symmetry=1, barrier=(36.7956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.878981,0.117816,-0.000190974,1.66731e-07,-5.78966e-11,-117831,36.8708], Tmin=(100,'K'), Tmax=(786.806,'K')), NASAPolynomial(coeffs=[12.0533,0.0404509,-2.13304e-05,4.22157e-09,-2.96957e-13,-119506,-20.1336], Tmin=(786.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-981.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(F)[C](F)C(F)C(F)F(7066)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {14,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 C u1 p0 c0 {6,S} {8,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-957.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([207,311,436,622,447,691,1074,1182,1276,1366,1311,1469,3045,3235,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,334,575,1197,1424,3202,180,180,180,2343.83],'cm^-1')),
        HinderedRotor(inertia=(0.552077,'amu*angstrom^2'), symmetry=1, barrier=(12.6933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41961,'amu*angstrom^2'), symmetry=1, barrier=(32.6396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551867,'amu*angstrom^2'), symmetry=1, barrier=(12.6885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41905,'amu*angstrom^2'), symmetry=1, barrier=(32.6268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16854,0.124704,-0.000206237,1.79386e-07,-6.1439e-11,-114982,36.8685], Tmin=(100,'K'), Tmax=(803.848,'K')), NASAPolynomial(coeffs=[13.3118,0.0388958,-2.04536e-05,4.02377e-09,-2.81345e-13,-116866,-27.0653], Tmin=(803.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-957.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[C]([CH]C(F)F)C(F)C(F)F(7067)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {13,S}
9  C u0 p0 c0 {4,S} {5,S} {11,S} {14,S}
10 C u1 p0 c0 {6,S} {7,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-982.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,3025,407.5,1350,352.5,180,184.347,188.253,2585.61],'cm^-1')),
        HinderedRotor(inertia=(0.216316,'amu*angstrom^2'), symmetry=1, barrier=(5.06224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218427,'amu*angstrom^2'), symmetry=1, barrier=(5.08235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39676,'amu*angstrom^2'), symmetry=1, barrier=(32.4897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40154,'amu*angstrom^2'), symmetry=1, barrier=(32.4956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.651683,0.112808,-0.000180349,1.59553e-07,-5.66938e-11,-118064,37.914], Tmin=(100,'K'), Tmax=(771.438,'K')), NASAPolynomial(coeffs=[10.6035,0.0431595,-2.29716e-05,4.57993e-09,-3.24127e-13,-119465,-11.288], Tmin=(771.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-982.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'F[CH][CH]C(F)(F)C(F)C(F)F(7068)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {13,S}
10 C u1 p0 c0 {8,S} {11,S} {14,S}
11 C u1 p0 c0 {6,S} {10,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-989.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09158,0.123718,-0.000207013,1.82969e-07,-6.35363e-11,-118839,36.8242], Tmin=(100,'K'), Tmax=(807.393,'K')), NASAPolynomial(coeffs=[12.3658,0.0404298,-2.14047e-05,4.21986e-09,-2.95263e-13,-120471,-21.8664], Tmin=(807.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-989.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'FC(F)[CH]C(F)(F)[CH]C(F)F(6822)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {11,S} {13,S}
9  C u0 p0 c0 {5,S} {6,S} {10,S} {12,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u1 p0 c0 {7,S} {8,S} {15,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1017.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,463,581,545,677,848,1004,1037,1149,1097,1177,1301,1447,1377,1455,3039,3185,3000,3050,390,425,1340,1360,335,370,203.261,203.404,203.646,2008.88],'cm^-1')),
        HinderedRotor(inertia=(0.236002,'amu*angstrom^2'), symmetry=1, barrier=(6.93571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236222,'amu*angstrom^2'), symmetry=1, barrier=(6.93426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00809,'amu*angstrom^2'), symmetry=1, barrier=(29.5431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236417,'amu*angstrom^2'), symmetry=1, barrier=(6.93507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.688581,0.112819,-0.000177591,1.53172e-07,-5.32228e-11,-122254,37.2828], Tmin=(100,'K'), Tmax=(758.947,'K')), NASAPolynomial(coeffs=[11.7287,0.0406014,-2.14722e-05,4.27693e-09,-3.02818e-13,-123943,-17.9132], Tmin=(758.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1017.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C(F)[CH]F(7069)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {10,S} {13,S}
10 C u1 p0 c0 {5,S} {9,S} {14,S}
11 C u1 p0 c0 {6,S} {8,S} {15,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-960.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,207,311,436,622,447,691,1074,1182,1276,1366,1311,1469,3045,3235,262,406,528,622,1148,1246,1368,1480,3164,3240,236.454,236.576,236.726,237.747],'cm^-1')),
        HinderedRotor(inertia=(0.236432,'amu*angstrom^2'), symmetry=1, barrier=(9.42517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236426,'amu*angstrom^2'), symmetry=1, barrier=(9.42524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607936,'amu*angstrom^2'), symmetry=1, barrier=(24.1019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833017,'amu*angstrom^2'), symmetry=1, barrier=(33.2461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52118,0.132448,-0.000221303,1.8987e-07,-6.38328e-11,-115385,35.0257], Tmin=(100,'K'), Tmax=(806.902,'K')), NASAPolynomial(coeffs=[15.4075,0.0361108,-1.91316e-05,3.76308e-09,-2.62675e-13,-117712,-40.5082], Tmin=(806.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-960.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    E0 = (-394.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-236.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (20.4526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (45.5608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-385.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-330.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-369.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-108.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-262.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-144.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-266.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-235.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-114.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-226.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-30.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (30.7571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-180.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-30.6113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-257.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-266.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-250.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-343.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-195.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-162.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-166.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-226.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-198.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['CHFCHF[Z](59)', 'FC(F)=CC(F)F(344)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['F[CH]C(F)C([C](F)F)C(F)F(6700)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(F)F-2(957)', 'F[CH]C(F)[C](F)F(818)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(804)', 'F[CH]C(F)(F)[CH]C(F)F(1054)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['FC(F)C1C(F)C(F)C1(F)F(6963)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['FC=C(F)C(F)(F)CC(F)F(7058)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['FCC(F)C(F)(F)C=C(F)F(6701)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'F[CH]C(F)C(F)=CC(F)F(7059)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(48.9307,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][CH]F(811)', 'FC(F)=CC(F)F(344)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(6.59488,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'FC=CC(F)(F)[CH]C(F)F(7060)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(63.2913,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CHFCHF[Z](59)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.26446e-11,'m^3/(mol*s)'), n=4.48095, Ea=(3.25025,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'FC=C(F)C(F)(F)[CH]C(F)F(7032)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'F[CH]C(F)C(F)(F)C=CF(7061)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.1183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[CH]C(F)C(F)(F)C=C(F)F(7062)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.5227,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][CH]F(811)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -6.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', 'F[CH]C=C(F)[CH]C(F)F(7063)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(18.9625,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[CH]C(F)=C(F)[CH]C(F)F(6936)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(240.455,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'F[CH]C(F)(F)[CH]C(F)F(1054)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['FC[C](F)C(F)(F)[CH]C(F)F(6851)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.30951e+09,'s^-1'), n=1.14834, Ea=(136.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['F[CH]C(F)C(F)(F)C[C](F)F(6827)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['F[CH][C](F)C(F)(F)CC(F)F(7064)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FCC(F)C(F)(F)[CH][C](F)F(7065)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(112.053,'s^-1'), n=2.49718, Ea=(42.3838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(F)[C](F)C(F)C(F)F(7066)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(166.592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]([CH]C(F)F)C(F)C(F)F(7067)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(225.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['F[CH][CH]C(F)(F)C(F)C(F)F(7068)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(227.576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    products = ['FC(F)[CH]C(F)(F)[CH]C(F)F(6822)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(F)C(F)(F)C(F)[CH]F(7069)'],
    products = ['F[CH]C(F)C(F)(F)[CH]C(F)F(6699)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #2019',
    isomers = [
        'F[CH]C(F)C(F)(F)[CH]C(F)F(6699)',
    ],
    reactants = [
        ('CHFCHF[Z](59)', 'FC(F)=CC(F)F(344)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2019',
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

