species(
    label = 'FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {4,S} {11,D} {14,S}
10 C u0 p0 c0 {5,S} {6,S} {12,D}
11 C u1 p0 c0 {7,S} {9,D}
12 C u1 p0 c0 {8,S} {10,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-606.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,164,312,561,654,898,1207,1299,3167,615,860,1140,1343,3152,562,600,623,1070,1265,1670,1700,300,440,180,180,1094.2],'cm^-1')),
        HinderedRotor(inertia=(0.245507,'amu*angstrom^2'), symmetry=1, barrier=(5.6447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242233,'amu*angstrom^2'), symmetry=1, barrier=(5.56941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244681,'amu*angstrom^2'), symmetry=1, barrier=(5.62569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279359,0.0980193,-0.00010752,1.40221e-09,6.1285e-11,-72818.6,35.2415], Tmin=(100,'K'), Tmax=(490.853,'K')), NASAPolynomial(coeffs=[10.6187,0.0421147,-2.33214e-05,4.71771e-09,-3.36438e-13,-74175.2,-10.7615], Tmin=(490.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC=C=C(F)F(5948)',
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
    label = 'F[CH]C(=C(F)F)C(F)[C]=C(F)F(7581)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {8,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-730.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.866452,0.119072,-0.00020217,1.81579e-07,-6.39712e-11,-87680,37.0293], Tmin=(100,'K'), Tmax=(806.176,'K')), NASAPolynomial(coeffs=[11.3357,0.0403048,-2.1706e-05,4.30431e-09,-3.02049e-13,-89055.3,-15.5372], Tmin=(806.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {3,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-756.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02382,0.122557,-0.000209261,1.86533e-07,-6.50473e-11,-90851.9,36.6752], Tmin=(100,'K'), Tmax=(809.493,'K')), NASAPolynomial(coeffs=[12.3022,0.0388206,-2.09506e-05,4.15048e-09,-2.90819e-13,-92423.3,-21.179], Tmin=(809.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = '[C]=CF(1054)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (404.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([194,682,905,1196,1383,3221],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.0070472,-1.17409e-06,-1.98501e-09,8.12281e-13,48709.9,8.54797], Tmin=(100,'K'), Tmax=(1284.59,'K')), NASAPolynomial(coeffs=[5.40185,0.00468,-2.11337e-06,4.24439e-10,-3.06832e-14,47991.3,-1.49755], Tmin=(1284.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'F[C](F)C(F)[C]=C(F)F(7771)',
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
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,190,488,555,1236,1407,562,600,623,1070,1265,1685,370,180,563.728],'cm^-1')),
        HinderedRotor(inertia=(0.188989,'amu*angstrom^2'), symmetry=1, barrier=(4.34523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180705,'amu*angstrom^2'), symmetry=1, barrier=(4.15475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930052,0.0756732,-0.000127381,1.17985e-07,-4.34855e-11,-65525.6,28.6743], Tmin=(100,'K'), Tmax=(771.224,'K')), NASAPolynomial(coeffs=[7.62005,0.0298117,-1.64701e-05,3.34225e-09,-2.38707e-13,-66225.5,0.289234], Tmin=(771.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-545.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[C]=C(F)F(1544)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41428e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.12,'K')), NASAPolynomial(coeffs=[7.28631,0.00323781,-1.64877e-06,3.54598e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'F[CH]C(F)(F)[C]=CF(6423)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u1 p0 c0 {3,S} {5,S} {9,S}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-348.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,334,575,1197,1424,3202,615,860,1140,1343,3152,1685,370,180,1014.91],'cm^-1')),
        HinderedRotor(inertia=(0.228324,'amu*angstrom^2'), symmetry=1, barrier=(5.24963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228576,'amu*angstrom^2'), symmetry=1, barrier=(5.25541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14737,0.0698491,-0.000114017,1.03508e-07,-3.74616e-11,-41866.5,27.3406], Tmin=(100,'K'), Tmax=(780.409,'K')), NASAPolynomial(coeffs=[7.49347,0.0277474,-1.46913e-05,2.93768e-09,-2.08232e-13,-42565.5,0.165871], Tmin=(780.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCd)(F1s)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'FC=C1C(=C(F)F)C(F)C1(F)F(8627)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {4,S} {9,D} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1002.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00333932,0.0942675,-0.000112104,6.95522e-08,-1.75415e-11,-120442,26.8292], Tmin=(100,'K'), Tmax=(954.574,'K')), NASAPolynomial(coeffs=[14.3535,0.0341344,-1.761e-05,3.55741e-09,-2.57358e-13,-123181,-41.7345], Tmin=(954.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1002.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCFH) + group(CdCFF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'FC=CC(F)(F)C(F)=C=C(F)F(8628)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {12,D}
10 C u0 p0 c0 {4,S} {8,D} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {9,D} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-935.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.876231,0.117239,-0.000192164,1.65775e-07,-5.64832e-11,-112342,34.155], Tmin=(100,'K'), Tmax=(796.596,'K')), NASAPolynomial(coeffs=[13.1134,0.0361344,-1.89986e-05,3.74307e-09,-2.62174e-13,-114226,-27.9932], Tmin=(796.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-935.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCFH) + group(CdCddFF) + group(Cdd-CdsCds)"""),
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
    label = 'FC=C=C(F)C(F)[C]=C(F)F(8629)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {2,S} {6,S} {11,D}
8  C u0 p0 c0 {4,S} {5,S} {10,D}
9  C u0 p0 c0 {3,S} {11,D} {13,S}
10 C u1 p0 c0 {6,S} {8,D}
11 C u0 p0 c0 {7,D} {9,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-482.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,562,600,623,1070,1265,113,247,382,1207,3490,1685,370,540,610,2055,180,180,180,522.447],'cm^-1')),
        HinderedRotor(inertia=(0.112847,'amu*angstrom^2'), symmetry=1, barrier=(2.59458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101175,'amu*angstrom^2'), symmetry=1, barrier=(2.32621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377775,0.108526,-0.000189261,1.7287e-07,-6.10904e-11,-57894.9,32.9244], Tmin=(100,'K'), Tmax=(829.209,'K')), NASAPolynomial(coeffs=[9.70847,0.0373577,-1.97952e-05,3.88084e-09,-2.697e-13,-58793.7,-9.1798], Tmin=(829.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCFF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'F[CH][C]=C(F)F(7343)',
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
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.805],'cm^-1')),
        HinderedRotor(inertia=(0.355657,'amu*angstrom^2'), symmetry=1, barrier=(8.17725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377451,-4.40201e-05,2.68133e-08,-6.58282e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.946,'K')), NASAPolynomial(coeffs=[8.46194,0.0130101,-6.31221e-06,1.26455e-09,-9.14158e-14,-25275.2,-12.1013], Tmin=(983.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'FC=[C]C(F)(F)C=C=C(F)F(8630)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {11,D} {12,S}
8  C u0 p0 c0 {3,S} {10,D} {13,S}
9  C u0 p0 c0 {4,S} {5,S} {11,D}
10 C u1 p0 c0 {6,S} {8,D}
11 C u0 p0 c0 {7,D} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-523.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,615,860,1140,1343,3152,94,120,354,641,825,1294,1685,370,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23306,'amu*angstrom^2'), symmetry=1, barrier=(28.3504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21916,'amu*angstrom^2'), symmetry=1, barrier=(28.0308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.448789,0.108593,-0.000184598,1.64588e-07,-5.73869e-11,-62786.1,32.7063], Tmin=(100,'K'), Tmax=(812.861,'K')), NASAPolynomial(coeffs=[11.1411,0.0351352,-1.87346e-05,3.69685e-09,-2.58489e-13,-64127.7,-17.4692], Tmin=(812.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(Cds_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,25474.2,-0.444973], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)=C=C(F)F(8631)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {3,S} {7,S} {12,D}
9  C u0 p0 c0 {4,S} {11,D} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {12,D}
11 C u1 p0 c0 {7,S} {9,D}
12 C u0 p0 c0 {8,D} {10,D}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-697.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,615,860,1140,1343,3152,94,120,354,641,825,1294,1685,370,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12592,'amu*angstrom^2'), symmetry=1, barrier=(25.8871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12713,'amu*angstrom^2'), symmetry=1, barrier=(25.915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.958739,0.122446,-0.000219148,1.979e-07,-6.86103e-11,-83736.3,34.923], Tmin=(100,'K'), Tmax=(838.104,'K')), NASAPolynomial(coeffs=[12.1766,0.0351888,-1.90124e-05,3.72851e-09,-2.58153e-13,-85075.3,-20.9802], Tmin=(838.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-697.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCFH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(F)(F)C(F)[C]=C(F)F(8632)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {3,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 C u0 p0 c0 {7,S} {11,T}
11 C u0 p0 c0 {10,T} {13,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-507.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([164,312,561,654,898,1207,1299,3167,154,355,414,641,686,1150,1196,562,600,623,1070,1265,1685,370,2175,525,750,770,3400,2100,180,2326.09],'cm^-1')),
        HinderedRotor(inertia=(0.519173,'amu*angstrom^2'), symmetry=1, barrier=(11.9368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80973,'amu*angstrom^2'), symmetry=1, barrier=(41.6093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80987,'amu*angstrom^2'), symmetry=1, barrier=(41.6124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.308339,0.103297,-0.000164592,1.40737e-07,-4.80411e-11,-60941.3,32.3664], Tmin=(100,'K'), Tmax=(776.937,'K')), NASAPolynomial(coeffs=[11.7949,0.0341493,-1.78927e-05,3.53465e-09,-2.48611e-13,-62615.7,-21.6419], Tmin=(776.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC#CC(F)(F)C(F)[C]=C(F)F(8633)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 C u0 p0 c0 {8,S} {12,T}
12 C u0 p0 c0 {6,S} {11,T}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-607.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([164,312,561,654,898,1207,1299,3167,154,355,414,641,686,1150,1196,562,600,623,1070,1265,1685,370,2175,525,239,401,1367,180,1717.77,1719.4],'cm^-1')),
        HinderedRotor(inertia=(0.347323,'amu*angstrom^2'), symmetry=1, barrier=(7.98564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04891,'amu*angstrom^2'), symmetry=1, barrier=(47.1084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346277,'amu*angstrom^2'), symmetry=1, barrier=(7.96158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.568701,0.112819,-0.000196131,1.79423e-07,-6.4053e-11,-72934.7,35.6185], Tmin=(100,'K'), Tmax=(808.989,'K')), NASAPolynomial(coeffs=[10.2527,0.0389354,-2.13529e-05,4.25566e-09,-2.99145e-13,-74018.7,-10.1731], Tmin=(808.989,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-607.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + group(Ct-CtCs) + group(CtCF) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC#CC(F)C(F)(F)[C]=CF(8634)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {12,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u1 p0 c0 {6,S} {8,D}
10 C u0 p0 c0 {7,S} {11,T}
11 C u0 p0 c0 {5,S} {10,T}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-410.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,233,378,609,1068,1270,1314,3037,615,860,1140,1343,3152,1685,370,2175,525,239,401,1367,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.311138,'amu*angstrom^2'), symmetry=1, barrier=(7.15367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311565,'amu*angstrom^2'), symmetry=1, barrier=(7.16349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0175153,'amu*angstrom^2'), symmetry=1, barrier=(45.3711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459656,0.109326,-0.000187272,1.68021e-07,-5.87232e-11,-49221.3,33.306], Tmin=(100,'K'), Tmax=(820.556,'K')), NASAPolynomial(coeffs=[10.7887,0.0359018,-1.9064e-05,3.74708e-09,-2.61133e-13,-50441.4,-14.9223], Tmin=(820.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFH) + group(Ct-CtCs) + group(CtCF) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = 'FC=C=C(F)[CH][C]=C(F)F(8635)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u1 p0 c0 {6,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {10,D}
7  C u0 p0 c0 {2,S} {3,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {12,S}
9  C u1 p0 c0 {5,S} {7,D}
10 C u0 p0 c0 {6,D} {8,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-188.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,132,421,945,1169,1268,562,600,623,1070,1265,113,247,382,1207,3490,1685,370,540,610,2055,543.596,543.687,543.735,543.912],'cm^-1')),
        HinderedRotor(inertia=(0.000570334,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000569963,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678635,0.075454,-8.80418e-05,5.36296e-08,-1.31167e-11,-22570.8,29.4835], Tmin=(100,'K'), Tmax=(991.245,'K')), NASAPolynomial(coeffs=[13.1887,0.0249724,-1.16514e-05,2.25335e-09,-1.59319e-13,-25050.9,-30.7603], Tmin=(991.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCFF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(Cds_S)"""),
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
    label = 'FC=[C]C(F)=C(F)[C]=C(F)F(8636)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {6,D} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {12,S}
10 C u1 p0 c0 {6,S} {8,D}
11 C u1 p0 c0 {7,S} {9,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-349.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([59,113,117,289,459,517,515,649,544,666,703,779,562,600,623,1070,1265,615,860,1140,1343,3152,1670,1700,300,440,286.988,288.657],'cm^-1')),
        HinderedRotor(inertia=(0.0019809,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00201108,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.200773,0.0904438,-0.000119712,8.19031e-08,-2.25755e-11,-41927.8,33.5529], Tmin=(100,'K'), Tmax=(880.469,'K')), NASAPolynomial(coeffs=[13.6288,0.0294406,-1.57859e-05,3.21379e-09,-2.32746e-13,-44292.4,-29.5201], Tmin=(880.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-349.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCFF) + radical(Cdj(Cd-F1sCd)(Cd-F1sH)) + radical(Cdj(Cd-F1sCd)(Cd-F1sF1s))"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)=C[C](F)F(8637)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {10,S} {13,S}
10 C u1 p0 c0 {4,S} {5,S} {9,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-758.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11055,0.125576,-0.000219507,1.97512e-07,-6.885e-11,-91032.2,36.4128], Tmin=(100,'K'), Tmax=(825.138,'K')), NASAPolynomial(coeffs=[12.1941,0.0386843,-2.08371e-05,4.10541e-09,-2.85961e-13,-92465.5,-20.5976], Tmin=(825.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-758.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=CC(F)(F)C(F)[C]=C(F)F(8638)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {7,S} {12,D} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 C u1 p0 c0 {6,S} {9,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-615.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.696005,0.114436,-0.000189869,1.70378e-07,-6.09645e-11,-73874,37.4821], Tmin=(100,'K'), Tmax=(771.715,'K')), NASAPolynomial(coeffs=[10.9861,0.0411531,-2.26805e-05,4.57038e-09,-3.24865e-13,-75297.9,-13.3935], Tmin=(771.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=CC(F)(F)C(F)=[C][C](F)F(8639)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {12,D}
10 C u0 p0 c0 {4,S} {8,D} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {9,D} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-737.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.951985,0.119884,-0.000199931,1.75766e-07,-6.10465e-11,-88473.8,36.695], Tmin=(100,'K'), Tmax=(790.56,'K')), NASAPolynomial(coeffs=[12.6636,0.0381969,-2.06587e-05,4.11426e-09,-2.89884e-13,-90226.7,-23.2633], Tmin=(790.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsF1s))"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(F)C=C(F)F(8640)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 C u1 p0 c0 {8,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-619.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.95459,0.120508,-0.000203711,1.81125e-07,-6.34478e-11,-74327.3,37.4951], Tmin=(100,'K'), Tmax=(793.427,'K')), NASAPolynomial(coeffs=[12.3222,0.0387653,-2.11776e-05,4.23035e-09,-2.98369e-13,-75968,-20.548], Tmin=(793.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[CH]C(F)=C(F)C(F)[C]=C(F)F(8641)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {10,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-707.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14691,0.126921,-0.000223229,2.01764e-07,-7.04327e-11,-84911.3,36.2746], Tmin=(100,'K'), Tmax=(830.545,'K')), NASAPolynomial(coeffs=[11.9333,0.0394882,-2.11859e-05,4.16072e-09,-2.89031e-13,-86241.2,-19.3272], Tmin=(830.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-707.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cd-CdF1s)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH][C]=C(F)C(F)C(F)=C(F)F(8642)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {3,S} {7,S} {12,D}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u1 p0 c0 {6,S} {12,S} {14,S}
12 C u1 p0 c0 {9,D} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-694.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12346,0.124831,-0.000213703,1.89626e-07,-6.56376e-11,-83319.1,36.2045], Tmin=(100,'K'), Tmax=(815.813,'K')), NASAPolynomial(coeffs=[12.7944,0.0381649,-2.04757e-05,4.03914e-09,-2.82057e-13,-84976.8,-24.35], Tmin=(815.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-694.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-CsF1s))"""),
)

species(
    label = 'FC=C(F)C(F)(F)C=[C][C](F)F(8643)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {7,S} {12,D} {13,S}
10 C u0 p0 c0 {4,S} {8,D} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {9,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-723.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03888,0.122039,-0.000205316,1.80516e-07,-6.24499e-11,-86862.7,36.9157], Tmin=(100,'K'), Tmax=(797.503,'K')), NASAPolynomial(coeffs=[13.0332,0.037572,-2.03245e-05,4.0401e-09,-2.84045e-13,-88665.6,-25.0211], Tmin=(797.503,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-723.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsH))"""),
)

species(
    label = 'FC=[C]C(F)(F)C=C(F)[C](F)F(8644)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {8,D} {10,S}
10 C u1 p0 c0 {4,S} {5,S} {9,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-743.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14935,0.127343,-0.000225572,2.05265e-07,-7.21361e-11,-89216.7,37.4481], Tmin=(100,'K'), Tmax=(827.007,'K')), NASAPolynomial(coeffs=[11.713,0.0400521,-2.17576e-05,4.29607e-09,-2.99405e-13,-90486.5,-16.9767], Tmin=(827.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-743.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C(F)[C]=C(F)F(8645)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {4,S} {7,S} {12,D}
10 C u0 p0 c0 {5,S} {6,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 C u1 p0 c0 {9,D} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-607.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,164,312,561,654,898,1207,1299,3167,246,474,533,1155,562,600,623,1070,1265,1685,370,3120,650,792.5,1650,180,180,1173.8],'cm^-1')),
        HinderedRotor(inertia=(0.189796,'amu*angstrom^2'), symmetry=1, barrier=(4.36377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189792,'amu*angstrom^2'), symmetry=1, barrier=(4.36368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189824,'amu*angstrom^2'), symmetry=1, barrier=(4.36442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.809838,0.116516,-0.000191898,1.69333e-07,-5.95705e-11,-72898.7,37.3005], Tmin=(100,'K'), Tmax=(770.868,'K')), NASAPolynomial(coeffs=[11.9871,0.0395673,-2.16458e-05,4.34752e-09,-3.08445e-13,-74558.3,-19.0746], Tmin=(770.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-607.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(F)C(F)=C(F)F(8646)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {5,S} {6,S} {9,D}
11 C u1 p0 c0 {8,S} {12,D}
12 C u1 p0 c0 {11,D} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-606.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,136,307,446,511,682,757,1180,1185,323,467,575,827,1418,182,240,577,636,1210,1413,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.441691,'amu*angstrom^2'), symmetry=1, barrier=(10.1553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44136,'amu*angstrom^2'), symmetry=1, barrier=(10.1477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441477,'amu*angstrom^2'), symmetry=1, barrier=(10.1504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26974,0.128293,-0.000222184,1.96977e-07,-6.79949e-11,-72727.2,37.4568], Tmin=(100,'K'), Tmax=(813.9,'K')), NASAPolynomial(coeffs=[13.6362,0.0368413,-2.01073e-05,3.98852e-09,-2.79186e-13,-74551,-27.6842], Tmin=(813.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = 'F[C]=[C]C(F)C(F)(F)C(F)=CF(8647)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {14,S}
11 C u1 p0 c0 {8,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-591.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,164,312,561,654,898,1207,1299,3167,323,467,575,827,1418,194,682,905,1196,1383,3221,1685,370,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.173122,'amu*angstrom^2'), symmetry=1, barrier=(3.98042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172775,'amu*angstrom^2'), symmetry=1, barrier=(3.97244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172678,'amu*angstrom^2'), symmetry=1, barrier=(3.97021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17056,0.12792,-0.000227614,2.07227e-07,-7.27868e-11,-70918.9,37.1406], Tmin=(100,'K'), Tmax=(827.489,'K')), NASAPolynomial(coeffs=[11.8656,0.0395849,-2.15914e-05,4.2678e-09,-2.97505e-13,-72209.5,-18.0449], Tmin=(827.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-591.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCCFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=C(F)C(F)C(F)(F)[C]=CF(8648)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {7,S} {12,D}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u1 p0 c0 {8,S} {10,D}
12 C u1 p0 c0 {6,S} {9,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-581.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,136,307,446,511,682,757,1180,1185,246,474,533,1155,615,860,1140,1343,3152,1685,370,167,640,1190,180,180,180,1347.74],'cm^-1')),
        HinderedRotor(inertia=(0.256081,'amu*angstrom^2'), symmetry=1, barrier=(5.88781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256284,'amu*angstrom^2'), symmetry=1, barrier=(5.89247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25631,'amu*angstrom^2'), symmetry=1, barrier=(5.89308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11166,0.124938,-0.000215909,1.93078e-07,-6.73142e-11,-69719,37.7656], Tmin=(100,'K'), Tmax=(812.826,'K')), NASAPolynomial(coeffs=[12.5819,0.0384043,-2.08849e-05,4.1428e-09,-2.90192e-13,-71312.6,-21.5693], Tmin=(812.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cd-CsF1s)(F1s))"""),
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
    E0 = (-222.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-96.9447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-128.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (242.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (230.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-214.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-129.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (33.3683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-137.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (5.38847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-100.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (15.9146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-12.2938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (86.6032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (29.0006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (204.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-8.66233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-40.2598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-18.7039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-80.8642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-43.9961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-86.1941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-35.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-48.9327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-100.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-41.4236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-25.5546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-17.8137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-28.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=C=C(F)F(5948)', 'FC=C=C(F)F(5948)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[CH]C(=C(F)F)C(F)[C]=C(F)F(7581)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=CF(1054)', 'F[C](F)C(F)[C]=C(F)F(7771)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=C(F)F(1544)', 'F[CH]C(F)(F)[C]=CF(6423)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=C1C(=C(F)F)C(F)C1(F)F(8627)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=CC(F)(F)C(F)=C=C(F)F(8628)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'FC=C=C(F)C(F)[C]=C(F)F(8629)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.4715,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC=C=C(F)F(5948)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(44.0711,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'FC=[C]C(F)(F)C=C=C(F)F(8630)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(72.1923,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'FC=[C]C(F)(F)C(F)=C=C(F)F(8631)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(1.94682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C#CC(F)(F)C(F)[C]=C(F)F(8632)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC#CC(F)(F)C(F)[C]=C(F)F(8633)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'FC#CC(F)C(F)(F)[C]=CF(8634)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C]=C(F)F(7343)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(46.7485,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', 'FC=C=C(F)[CH][C]=C(F)F(8635)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(18.3782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'FC=[C]C(F)=C(F)[C]=C(F)F(8636)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(238.563,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=[C]C(F)(F)C(F)=C[C](F)F(8637)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[C]=CC(F)(F)C(F)[C]=C(F)F(8638)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=CC(F)(F)C(F)=[C][C](F)F(8639)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[C]=[C]C(F)(F)C(F)C=C(F)F(8640)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.02951e+08,'s^-1'), n=1.2385, Ea=(178.845,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_single] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[CH]C(F)=C(F)C(F)[C]=C(F)F(8641)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(136.647,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['F[CH][C]=C(F)C(F)C(F)=C(F)F(8642)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(186.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=C(F)C(F)(F)C=[C][C](F)F(8643)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(173.909,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    products = ['FC=[C]C(F)(F)C=C(F)[C](F)F(8644)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(122.129,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(F)C(F)(F)C(F)[C]=C(F)F(8645)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]C(F)(F)C(F)C(F)=C(F)F(8646)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(197.029,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C]=[C]C(F)C(F)(F)C(F)=CF(8647)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(189.69,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[C]=C(F)C(F)C(F)(F)[C]=CF(8648)'],
    products = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2841',
    isomers = [
        'FC=[C]C(F)(F)C(F)[C]=C(F)F(7584)',
    ],
    reactants = [
        ('FC=C=C(F)F(5948)', 'FC=C=C(F)F(5948)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2841',
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

