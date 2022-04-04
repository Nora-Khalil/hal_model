species(
    label = 'FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {15,S}
11 C u1 p0 c0 {8,S} {10,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1021.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,180,2241.89],'cm^-1')),
        HinderedRotor(inertia=(0.203212,'amu*angstrom^2'), symmetry=1, barrier=(4.67223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203187,'amu*angstrom^2'), symmetry=1, barrier=(4.67166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204125,'amu*angstrom^2'), symmetry=1, barrier=(4.69323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55209,'amu*angstrom^2'), symmetry=1, barrier=(35.6856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32681,0.129436,-0.000211359,1.88309e-07,-6.72684e-11,-122624,42.6331], Tmin=(100,'K'), Tmax=(763.254,'K')), NASAPolynomial(coeffs=[11.9671,0.0474391,-2.59866e-05,5.23497e-09,-3.72467e-13,-124294,-15.5587], Tmin=(763.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1021.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC=C=C(F)F(1325)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C(C(F)F)C(F)[C]=C(F)F(6744)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {13,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {16,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1024.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,562,600,623,1070,1265,1685,370,198.829,198.908,198.993,1741.43],'cm^-1')),
        HinderedRotor(inertia=(0.174321,'amu*angstrom^2'), symmetry=1, barrier=(4.89806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174391,'amu*angstrom^2'), symmetry=1, barrier=(4.89899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1745,'amu*angstrom^2'), symmetry=1, barrier=(4.8984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16607,'amu*angstrom^2'), symmetry=1, barrier=(32.7679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3374.93,'J/mol'), sigma=(5.29414,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.16 K, Pc=51.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30365,0.128321,-0.00020668,1.81642e-07,-6.41986e-11,-123062,42.2379], Tmin=(100,'K'), Tmax=(757.094,'K')), NASAPolynomial(coeffs=[12.4038,0.046301,-2.51612e-05,5.05854e-09,-3.59765e-13,-124863,-18.2616], Tmin=(757.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1024.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {8,S} {9,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1145.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50681,0.133506,-0.000219857,1.94354e-07,-6.84139e-11,-137560,41.5557], Tmin=(100,'K'), Tmax=(775.942,'K')), NASAPolynomial(coeffs=[12.9424,0.0458445,-2.49268e-05,4.99304e-09,-3.53626e-13,-139406,-21.9322], Tmin=(775.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1145.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Cs_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[C](F)C(F)[C]=C(F)F(3149)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {3,S} {6,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-545.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,190,488,555,1236,1407,562,600,623,1070,1265,1685,370,180,549.592],'cm^-1')),
        HinderedRotor(inertia=(0.178029,'amu*angstrom^2'), symmetry=1, barrier=(4.09324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184635,'amu*angstrom^2'), symmetry=1, barrier=(4.24512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930086,0.0756727,-0.000127379,1.17982e-07,-4.34838e-11,-65525.6,28.6742], Tmin=(100,'K'), Tmax=(771.236,'K')), NASAPolynomial(coeffs=[7.62,0.0298118,-1.64702e-05,3.34227e-09,-2.38708e-13,-66225.5,0.289503], Tmin=(771.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-545.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[C]=C(F)F(1205)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u2 p0 c0 {3,D}
"""),
    E0 = (195.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([182,240,577,636,1210,1413],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41427e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.13,'K')), NASAPolynomial(coeffs=[7.28631,0.0032378,-1.64877e-06,3.54597e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
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
    label = 'FC(F)=C1C(F)C(F)(F)C1C(F)F(7289)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {11,S} {12,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {9,S} {12,S} {15,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {16,S}
12 C u0 p0 c0 {8,S} {10,S} {13,D}
13 C u0 p0 c0 {6,S} {7,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1339.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.511264,0.108287,-0.000131759,8.56852e-08,-2.29223e-11,-160971,31.0428], Tmin=(100,'K'), Tmax=(897.543,'K')), NASAPolynomial(coeffs=[14.1482,0.0429547,-2.25715e-05,4.58358e-09,-3.32144e-13,-163602,-38.0961], Tmin=(897.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1339.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FC(F)=C=C(F)C(F)(F)CC(F)F(7290)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {16,S}
11 C u0 p0 c0 {5,S} {9,S} {13,D}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u0 p0 c0 {11,D} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1280.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60419,0.134228,-0.000216725,1.84927e-07,-6.26968e-11,-153834,37.8334], Tmin=(100,'K'), Tmax=(779.57,'K')), NASAPolynomial(coeffs=[14.828,0.0414325,-2.18554e-05,4.32441e-09,-3.04163e-13,-156138,-35.6975], Tmin=(779.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1280.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = 'FC(F)=CC(F)C(F)(F)C=C(F)F(6749)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {15,S}
11 C u0 p0 c0 {9,S} {13,D} {16,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1329.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13808,0.123411,-0.000186177,1.52614e-07,-5.07154e-11,-159758,37.1396], Tmin=(100,'K'), Tmax=(734.203,'K')), NASAPolynomial(coeffs=[13.4936,0.0436854,-2.32715e-05,4.67302e-09,-3.33463e-13,-161907,-28.9268], Tmin=(734.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1329.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF)"""),
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
    label = 'FC(F)=[C]C(F)C(F)=CC(F)F(7291)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {14,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-896.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.376538,'amu*angstrom^2'), symmetry=1, barrier=(8.65734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375587,'amu*angstrom^2'), symmetry=1, barrier=(8.63549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937984,'amu*angstrom^2'), symmetry=1, barrier=(21.5661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.892797,0.120453,-0.000203896,1.85278e-07,-6.59798e-11,-107712,36.7154], Tmin=(100,'K'), Tmax=(810.724,'K')), NASAPolynomial(coeffs=[10.2761,0.0441383,-2.34568e-05,4.63302e-09,-3.2448e-13,-108826,-10.5258], Tmin=(810.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-896.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Cds_S)"""),
)

species(
    label = 'F[CH][C]=C(F)F(2842)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-200.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.793],'cm^-1')),
        HinderedRotor(inertia=(0.35565,'amu*angstrom^2'), symmetry=1, barrier=(8.17708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377452,-4.40203e-05,2.68135e-08,-6.58286e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.929,'K')), NASAPolynomial(coeffs=[8.46193,0.0130101,-6.31222e-06,1.26456e-09,-9.14162e-14,-25275.2,-12.1012], Tmin=(983.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'FC(F)=C=CC(F)(F)[CH]C(F)F(7292)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u1 p0 c0 {7,S} {8,S} {14,S}
10 C u0 p0 c0 {7,S} {12,D} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {10,D} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-911.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.236997,'amu*angstrom^2'), symmetry=1, barrier=(5.44902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237108,'amu*angstrom^2'), symmetry=1, barrier=(5.45159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237296,'amu*angstrom^2'), symmetry=1, barrier=(5.45591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.932433,0.11955,-0.000195221,1.72438e-07,-6.07607e-11,-109494,37.5891], Tmin=(100,'K'), Tmax=(775.741,'K')), NASAPolynomial(coeffs=[11.7854,0.0421517,-2.27063e-05,4.53833e-09,-3.21205e-13,-111112,-18.2453], Tmin=(775.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-911.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(Cs_S)"""),
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
    label = 'FC(F)=C=C(F)C(F)(F)[CH]C(F)F(7293)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u1 p0 c0 {8,S} {9,S} {15,S}
11 C u0 p0 c0 {5,S} {8,S} {13,D}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u0 p0 c0 {11,D} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1086.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,145,326,398,834,1303,94,120,354,641,825,1294,540,610,2055,180,2028.4],'cm^-1')),
        HinderedRotor(inertia=(0.303204,'amu*angstrom^2'), symmetry=1, barrier=(6.97126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303936,'amu*angstrom^2'), symmetry=1, barrier=(6.9881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304124,'amu*angstrom^2'), symmetry=1, barrier=(6.99242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46194,0.133645,-0.000230658,2.06898e-07,-7.2435e-11,-130444,39.875], Tmin=(100,'K'), Tmax=(810.542,'K')), NASAPolynomial(coeffs=[12.9118,0.0420418,-2.28861e-05,4.54616e-09,-3.1885e-13,-132095,-22.2627], Tmin=(810.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1086.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + radical(Cs_S)"""),
)

species(
    label = 'FC=CC(F)(F)C(F)[C]=C(F)F(7294)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {4,S} {9,D} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {8,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-865.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,164,312,561,654,898,1207,1299,3167,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,562,600,623,1070,1265,1685,370,180,180,1125.7],'cm^-1')),
        HinderedRotor(inertia=(0.234681,'amu*angstrom^2'), symmetry=1, barrier=(5.39577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234927,'amu*angstrom^2'), symmetry=1, barrier=(5.40144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235498,'amu*angstrom^2'), symmetry=1, barrier=(5.41457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.173221,0.10424,-0.000138624,8.88369e-08,-1.53449e-11,-103976,35.1933], Tmin=(100,'K'), Tmax=(571.899,'K')), NASAPolynomial(coeffs=[10.722,0.0442816,-2.39752e-05,4.86143e-09,-3.49406e-13,-105488,-13.6046], Tmin=(571.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-865.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=[C]C(F)C(F)(F)C=C(F)F(7295)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {11,D} {15,S}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1066.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,164,312,561,654,898,1207,1299,3167,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,562,600,623,1070,1265,1685,370,180,180,1177.04],'cm^-1')),
        HinderedRotor(inertia=(0.200274,'amu*angstrom^2'), symmetry=1, barrier=(4.6047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205172,'amu*angstrom^2'), symmetry=1, barrier=(4.71732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203106,'amu*angstrom^2'), symmetry=1, barrier=(4.66981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.329723,0.110323,-0.000142333,6.98349e-08,7.26334e-12,-128152,36.813], Tmin=(100,'K'), Tmax=(537.222,'K')), NASAPolynomial(coeffs=[11.6271,0.0447351,-2.46443e-05,5.00133e-09,-3.58478e-13,-129775,-16.5903], Tmin=(537.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1066.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC#CC(F)C(F)(F)[CH]C(F)F(7296)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {14,S}
10 C u1 p0 c0 {7,S} {9,S} {15,S}
11 C u0 p0 c0 {8,S} {12,T}
12 C u0 p0 c0 {6,S} {11,T}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-825.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,233,378,609,1068,1270,1314,3037,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,2175,525,239,401,1367,221.486,221.54,221.7,1540.67],'cm^-1')),
        HinderedRotor(inertia=(0.179464,'amu*angstrom^2'), symmetry=1, barrier=(6.24369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1793,'amu*angstrom^2'), symmetry=1, barrier=(6.24377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1793,'amu*angstrom^2'), symmetry=1, barrier=(6.24055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32103,'amu*angstrom^2'), symmetry=1, barrier=(46.0295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17674,0.126376,-0.000214271,1.9069e-07,-6.62576e-11,-99064.4,37.712], Tmin=(100,'K'), Tmax=(819.133,'K')), NASAPolynomial(coeffs=[12.2039,0.041117,-2.16687e-05,4.25052e-09,-2.96048e-13,-100588,-20.0934], Tmin=(819.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-825.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
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
    label = 'FC(F)=[C]C=C(F)[CH]C(F)F(7297)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
7  C u1 p0 c0 {6,S} {8,S} {13,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {11,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-604.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,271,519,563,612,1379,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.852582,'amu*angstrom^2'), symmetry=1, barrier=(19.6025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429435,'amu*angstrom^2'), symmetry=1, barrier=(61.3445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476588,'amu*angstrom^2'), symmetry=1, barrier=(10.9577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0152323,0.0927973,-0.000110832,6.85175e-08,-1.70407e-11,-72604.7,31.0922], Tmin=(100,'K'), Tmax=(973.305,'K')), NASAPolynomial(coeffs=[15.0828,0.0307482,-1.52049e-05,3.01667e-09,-2.16204e-13,-75543.7,-41.3383], Tmin=(973.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + radical(Allyl_S) + radical(Cdj(Cd-CdH)(Cd-F1sF1s))"""),
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
    label = 'FC(F)=[C]C(F)=C(F)[CH]C(F)F(7298)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
8  C u1 p0 c0 {7,S} {9,S} {14,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {10,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-778.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,271,519,563,612,1379,86,203,488,582,605,741,562,600,623,1070,1265,1685,370,201.475,203.817,207.931],'cm^-1')),
        HinderedRotor(inertia=(0.329081,'amu*angstrom^2'), symmetry=1, barrier=(10.0697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.042309,'amu*angstrom^2'), symmetry=1, barrier=(39.6785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0135222,'amu*angstrom^2'), symmetry=1, barrier=(39.6642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254114,0.100169,-0.000123486,7.81086e-08,-1.99197e-11,-93523.1,34.4675], Tmin=(100,'K'), Tmax=(948.01,'K')), NASAPolynomial(coeffs=[15.6277,0.0331578,-1.74583e-05,3.54754e-09,-2.57301e-13,-96534.4,-41.3053], Tmin=(948.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-778.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + radical(Allyl_S) + radical(Cdj(Cd-F1sCd)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C](F)C=C(F)C(F)(F)[CH]C(F)F(7299)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u1 p0 c0 {8,S} {9,S} {15,S}
11 C u0 p0 c0 {5,S} {8,S} {12,D}
12 C u0 p0 c0 {11,D} {13,S} {16,S}
13 C u1 p0 c0 {6,S} {7,S} {12,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1146.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60637,0.136686,-0.000230703,2.06133e-07,-7.25457e-11,-137740,41.3385], Tmin=(100,'K'), Tmax=(796.463,'K')), NASAPolynomial(coeffs=[12.8889,0.0456096,-2.47539e-05,4.93352e-09,-3.47542e-13,-139469,-21.6547], Tmin=(796.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1146.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC(F)(F)C(F)[C]=C(F)F(7300)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {3,S} {8,S} {13,S} {16,S}
11 C u1 p0 c0 {4,S} {5,S} {9,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {10,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1028.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53698,0.136745,-0.000236032,2.16828e-07,-7.80877e-11,-123468,41.65], Tmin=(100,'K'), Tmax=(800.312,'K')), NASAPolynomial(coeffs=[11.1971,0.049179,-2.70774e-05,5.42129e-09,-3.82602e-13,-124741,-12.162], Tmin=(800.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1028.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C]([C]=C(F)F)C(F)(F)CC(F)F(7301)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {16,S}
11 C u1 p0 c0 {5,S} {9,S} {13,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {11,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1090.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39175,0.129867,-0.000200684,1.6727e-07,-5.64289e-11,-130998,39.9635], Tmin=(100,'K'), Tmax=(723.794,'K')), NASAPolynomial(coeffs=[13.9863,0.0448657,-2.4493e-05,4.95576e-09,-3.54793e-13,-133224,-29.2529], Tmin=(723.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1090.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCdCsF1s) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C](F)[CH]C(F)(F)C(F)C=C(F)F(7302)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
10 C u1 p0 c0 {8,S} {12,S} {16,S}
11 C u0 p0 c0 {9,S} {13,D} {15,S}
12 C u1 p0 c0 {4,S} {5,S} {10,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1083.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59453,0.137401,-0.000235888,2.13463e-07,-7.55646e-11,-130137,42.4936], Tmin=(100,'K'), Tmax=(808.072,'K')), NASAPolynomial(coeffs=[12.1863,0.0465526,-2.52376e-05,5.01478e-09,-3.52114e-13,-131626,-16.4821], Tmin=(808.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1083.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](C(F)[C]=C(F)F)C(F)C(F)F(7303)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {11,S} {13,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {16,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-989.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,174,267,591,721,1107,1278,1348,3273,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,562,600,623,1070,1265,1685,370,180,180,180,1120.97],'cm^-1')),
        HinderedRotor(inertia=(0.226805,'amu*angstrom^2'), symmetry=1, barrier=(5.21469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225844,'amu*angstrom^2'), symmetry=1, barrier=(5.1926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225727,'amu*angstrom^2'), symmetry=1, barrier=(5.18991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86222,'amu*angstrom^2'), symmetry=1, barrier=(42.8162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38866,0.131479,-0.000218017,1.95926e-07,-7.0024e-11,-118787,42.6481], Tmin=(100,'K'), Tmax=(780.821,'K')), NASAPolynomial(coeffs=[11.7168,0.0479682,-2.61321e-05,5.23722e-09,-3.70937e-13,-120334,-14.1393], Tmin=(780.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-989.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCsCsF1s) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C]([CH]C(F)F)C(F)C(F)=C(F)F(7304)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {12,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {11,S} {15,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {16,S}
12 C u0 p0 c0 {4,S} {8,S} {13,D}
13 C u0 p0 c0 {6,S} {7,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1027.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54419,0.138016,-0.000242324,2.24456e-07,-8.07346e-11,-123452,43.3477], Tmin=(100,'K'), Tmax=(816.224,'K')), NASAPolynomial(coeffs=[10.5845,0.0497406,-2.71015e-05,5.38419e-09,-3.77522e-13,-124471,-6.81878], Tmin=(816.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1027.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'FC(F)=[C][CH]C(F)(F)C(F)C(F)F(7305)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
11 C u1 p0 c0 {9,S} {13,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {11,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1088.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10493,0.121557,-0.00016769,1.21401e-07,-3.53948e-11,-130788,38.1498], Tmin=(100,'K'), Tmax=(834.911,'K')), NASAPolynomial(coeffs=[15.9429,0.0398771,-2.09348e-05,4.21178e-09,-3.02397e-13,-133634,-41.0192], Tmin=(834.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1088.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Allyl_S) + radical(Cdj(Cs-CsHH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=C(F)[CH]C(F)(F)[CH]C(F)F(7306)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u1 p0 c0 {8,S} {9,S} {15,S}
11 C u1 p0 c0 {8,S} {12,S} {16,S}
12 C u0 p0 c0 {5,S} {11,S} {13,D}
13 C u0 p0 c0 {6,S} {7,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1131.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59477,0.134329,-0.000217211,1.86566e-07,-6.3865e-11,-135945,39.5962], Tmin=(100,'K'), Tmax=(770.929,'K')), NASAPolynomial(coeffs=[14.4677,0.0427086,-2.28337e-05,4.54575e-09,-3.21009e-13,-138176,-32.1215], Tmin=(770.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1131.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Allyl_S)"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C(F)[C]=C(F)F(7307)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {4,S} {8,S} {13,S} {15,S}
11 C u1 p0 c0 {5,S} {9,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {10,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-992.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,334,575,1197,1424,3202,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0997574,'amu*angstrom^2'), symmetry=1, barrier=(2.29362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989395,'amu*angstrom^2'), symmetry=1, barrier=(2.27481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229278,'amu*angstrom^2'), symmetry=1, barrier=(5.27154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7022,'amu*angstrom^2'), symmetry=1, barrier=(39.1369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73774,0.139175,-0.00023288,2.06085e-07,-7.22431e-11,-119189,41.486], Tmin=(100,'K'), Tmax=(781.358,'K')), NASAPolynomial(coeffs=[13.8089,0.0451899,-2.48142e-05,4.97753e-09,-3.52354e-13,-121179,-26.8693], Tmin=(781.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-992.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH][CH]C(F)(F)C(F)C(F)=C(F)F(7308)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
10 C u1 p0 c0 {8,S} {12,S} {15,S}
11 C u0 p0 c0 {4,S} {9,S} {13,D}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1034.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9605,0.148632,-0.000267901,2.46442e-07,-8.69989e-11,-124228,42.1745], Tmin=(100,'K'), Tmax=(834.658,'K')), NASAPolynomial(coeffs=[12.2428,0.0471979,-2.56468e-05,5.05141e-09,-3.50971e-13,-125437,-16.8188], Tmin=(834.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1034.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFHH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'F[C]=[C]C(F)C(F)(F)C(F)C(F)F(7309)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {4,S} {8,S} {12,S} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {16,S}
12 C u1 p0 c0 {10,S} {13,D}
13 C u1 p0 c0 {7,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-959.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,164,312,561,654,898,1207,1299,3167,235,523,627,1123,1142,1372,1406,3097,1685,370,167,640,1190,278.299,278.303,278.307],'cm^-1')),
        HinderedRotor(inertia=(0.136467,'amu*angstrom^2'), symmetry=1, barrier=(7.50107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136474,'amu*angstrom^2'), symmetry=1, barrier=(7.50099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136467,'amu*angstrom^2'), symmetry=1, barrier=(7.50093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636149,'amu*angstrom^2'), symmetry=1, barrier=(34.9592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90663,0.14297,-0.000241691,2.12055e-07,-7.30206e-11,-115218,40.8515], Tmin=(100,'K'), Tmax=(801.718,'K')), NASAPolynomial(coeffs=[14.9077,0.042589,-2.30275e-05,4.57069e-09,-3.20843e-13,-117384,-33.2469], Tmin=(801.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-959.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=C(F)C(F)C(F)(F)[CH]C(F)F(7310)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {12,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {15,S}
11 C u1 p0 c0 {8,S} {10,S} {16,S}
12 C u0 p0 c0 {6,S} {9,S} {13,D}
13 C u1 p0 c0 {7,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-995.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,284,328,853,1146,1135,1297,3239,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,246,474,533,1155,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82803,0.141979,-0.000242872,2.15693e-07,-7.48222e-11,-119562,42.1691], Tmin=(100,'K'), Tmax=(812.145,'K')), NASAPolynomial(coeffs=[13.9954,0.0436228,-2.34916e-05,4.64671e-09,-3.25147e-13,-121459,-26.7306], Tmin=(812.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-995.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Cdj(Cd-CsF1s)(F1s))"""),
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
    E0 = (-354.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-196.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-259.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (52.8239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (98.3792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-345.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-290.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-314.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-97.5264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-277.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-107.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-246.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-207.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-74.7309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-186.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-44.6875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-63.8861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (71.8143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-146.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-171.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-226.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-210.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-196.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-155.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-129.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-159.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-222.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-158.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-161.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-124.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-159.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC=C=C(F)F(1325)', 'FC(F)=CC(F)F(344)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C](F)C(C(F)F)C(F)[C]=C(F)F(6744)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(F)F-2(957)', 'F[C](F)C(F)[C]=C(F)F(3149)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=C(F)F(1205)', 'F[CH]C(F)(F)[CH]C(F)F(1054)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC(F)=C1C(F)C(F)(F)C1C(F)F(7289)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC(F)=C=C(F)C(F)(F)CC(F)F(7290)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC(F)=CC(F)C(F)(F)C=C(F)F(6749)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'FC(F)=[C]C(F)C(F)=CC(F)F(7291)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.5968,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH][C]=C(F)F(2842)', 'FC(F)=CC(F)F(344)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(49.4574,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC(F)=C=CC(F)(F)[CH]C(F)F(7292)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(64.0574,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC=C=C(F)F(1325)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.13223e-11,'m^3/(mol*s)'), n=4.48095, Ea=(6.07887,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC(F)=C=C(F)C(F)(F)[CH]C(F)F(7293)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'FC=CC(F)(F)C(F)[C]=C(F)F(7294)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.1183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'FC(F)=[C]C(F)C(F)(F)C=C(F)F(7295)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.5227,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'FC#CC(F)C(F)(F)[CH]C(F)F(7296)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(24.2568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(78)', 'FC(F)=[C]C=C(F)[CH]C(F)F(7297)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(18.5436,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'FC(F)=[C]C(F)=C(F)[CH]C(F)F(7298)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(246.477,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C](F)C=C(F)C(F)(F)[CH]C(F)F(7299)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C](F)CC(F)(F)C(F)[C]=C(F)F(7300)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C]([C]=C(F)F)C(F)(F)CC(F)F(7301)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C](F)[CH]C(F)(F)C(F)C=C(F)F(7302)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.01114e+10,'s^-1'), n=0.8095, Ea=(157.607,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C](C(F)[C]=C(F)F)C(F)C(F)F(7303)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(166.592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[C]([CH]C(F)F)C(F)C(F)=C(F)F(7304)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(224.656,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC(F)=[C][CH]C(F)(F)C(F)C(F)F(7305)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(195.088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['FC(F)=C(F)[CH]C(F)(F)[CH]C(F)F(7306)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(132.101,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(F)C(F)(F)C(F)[C]=C(F)F(7307)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction29',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[CH][CH]C(F)(F)C(F)C(F)=C(F)F(7308)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(192.72,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[C]=[C]C(F)C(F)(F)C(F)C(F)F(7309)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(168.063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]=C(F)C(F)C(F)(F)[CH]C(F)F(7310)'],
    products = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2059',
    isomers = [
        'FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)',
    ],
    reactants = [
        ('FC=C=C(F)F(1325)', 'FC(F)=CC(F)F(344)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2059',
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

