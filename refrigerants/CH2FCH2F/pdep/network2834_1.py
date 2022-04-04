species(
    label = 'FC=[C]C(F)(F)C(F)[C](F)F(7577)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-781.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,136,307,446,511,682,757,1180,1185,190,488,555,1236,1407,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.202137,'amu*angstrom^2'), symmetry=1, barrier=(4.64753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202009,'amu*angstrom^2'), symmetry=1, barrier=(4.64459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786646,'amu*angstrom^2'), symmetry=1, barrier=(18.0865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.665923,0.113872,-0.000196018,1.74885e-07,-6.0928e-11,-93884.8,35.0989], Tmin=(100,'K'), Tmax=(810.044,'K')), NASAPolynomial(coeffs=[11.9235,0.0350032,-1.9044e-05,3.78094e-09,-2.6509e-13,-95376.4,-19.6032], Tmin=(810.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-781.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'CHFCF2(55)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-511.455,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(82.003,'amu')),
        NonlinearRotor(inertia=([47.283,130.848,178.131],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([226.38,309.801,488.171,585.738,628.102,776.692,950.625,1195.27,1295.73,1386.74,1841.68,3239.84],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93201,0.00413887,7.03233e-05,-1.49827e-07,9.42397e-11,-61509.2,9.50863], Tmin=(10,'K'), Tmax=(540.182,'K')), NASAPolynomial(coeffs=[3.82032,0.0197002,-1.38027e-05,4.49179e-09,-5.49698e-13,-61712.1,7.9889], Tmin=(540.182,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[CH]C(F)(F)C(F)(F)[C]=CF(7574)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {11,S}
9  C u1 p0 c0 {5,S} {7,S} {12,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-782.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,136,307,446,511,682,757,1180,1185,334,575,1197,1424,3202,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.223118,'amu*angstrom^2'), symmetry=1, barrier=(5.12992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223085,'amu*angstrom^2'), symmetry=1, barrier=(5.12916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223041,'amu*angstrom^2'), symmetry=1, barrier=(5.12815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.817149,0.119849,-0.000216609,1.98291e-07,-6.95912e-11,-93955.1,34.8401], Tmin=(100,'K'), Tmax=(836.504,'K')), NASAPolynomial(coeffs=[11.1403,0.0363893,-1.98252e-05,3.90285e-09,-2.70879e-13,-95036.2,-15.2172], Tmin=(836.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-782.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)[C](F)F(7575)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 C u1 p0 c0 {4,S} {8,S} {13,S}
11 C u0 p0 c0 {5,S} {6,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-865.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.241796,0.102091,-0.000159986,1.36705e-07,-4.73092e-11,-103940,34.4087], Tmin=(100,'K'), Tmax=(729.39,'K')), NASAPolynomial(coeffs=[11.3426,0.0361954,-1.96038e-05,3.9462e-09,-2.81373e-13,-105567,-17.3925], Tmin=(729.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-865.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-F1sCdH)(F1s)(F1s)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C]F(744)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
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
    label = 'F[C](F)C(F)[C](F)F(6178)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u1 p0 c0 {4,S} {5,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-711.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.130621,'amu*angstrom^2'), symmetry=1, barrier=(3.00323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131003,'amu*angstrom^2'), symmetry=1, barrier=(3.01201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03364,0.0738655,-0.000133101,1.21989e-07,-4.27831e-11,-85493.9,26.1787], Tmin=(100,'K'), Tmax=(838.764,'K')), NASAPolynomial(coeffs=[8.34905,0.0225506,-1.19531e-05,2.34629e-09,-1.62739e-13,-86143.2,-4.38292], Tmin=(838.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC=C1C(F)(F)C(F)C1(F)F(8345)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1096.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.324607,0.0888143,-0.000115769,8.21427e-08,-2.40106e-11,-131695,25.8074], Tmin=(100,'K'), Tmax=(824.857,'K')), NASAPolynomial(coeffs=[11.312,0.0355337,-1.88801e-05,3.83581e-09,-2.77435e-13,-133507,-25.0848], Tmin=(824.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1096.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FC=CC(F)(F)C(F)=C(F)F(8346)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1078.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.537225,0.109779,-0.000183909,1.60906e-07,-5.52154e-11,-129540,30.6475], Tmin=(100,'K'), Tmax=(809.45,'K')), NASAPolynomial(coeffs=[12.1484,0.0336087,-1.77726e-05,3.49756e-09,-2.44322e-13,-131152,-25.1434], Tmin=(809.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1078.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + group(CdCFH)"""),
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
    label = 'FC=[C]C(F)(F)C=C(F)F(7698)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-688.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2329,'amu*angstrom^2'), symmetry=1, barrier=(28.3469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23147,'amu*angstrom^2'), symmetry=1, barrier=(28.314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0644733,0.0962131,-0.000162631,1.45663e-07,-5.1156e-11,-82615.7,29.2123], Tmin=(100,'K'), Tmax=(809.074,'K')), NASAPolynomial(coeffs=[9.93846,0.0326351,-1.73913e-05,3.43652e-09,-2.40684e-13,-83730.3,-13.3462], Tmin=(809.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-688.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cds_S)"""),
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
    label = 'FC=[C]C(F)(F)C(F)=C(F)F(8347)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 C u0 p0 c0 {6,S} {11,D} {12,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-840.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,182,240,577,636,1210,1413,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17042,'amu*angstrom^2'), symmetry=1, barrier=(26.9102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17293,'amu*angstrom^2'), symmetry=1, barrier=(26.968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.60733,0.114828,-0.000210296,1.92199e-07,-6.69704e-11,-100935,31.3717], Tmin=(100,'K'), Tmax=(845.239,'K')), NASAPolynomial(coeffs=[11.1668,0.0327433,-1.78344e-05,3.49469e-09,-2.41289e-13,-101984,-17.8807], Tmin=(845.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-840.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + group(CdCFH) + radical(Cds_S)"""),
)

species(
    label = 'FC=C=C(F)C(F)[C](F)F(8348)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {10,D}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {7,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-617.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,145,326,398,834,1303,190,488,555,1236,1407,113,247,382,1207,3490,540,610,2055,1058.92,1059.13],'cm^-1')),
        HinderedRotor(inertia=(0.101153,'amu*angstrom^2'), symmetry=1, barrier=(2.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101182,'amu*angstrom^2'), symmetry=1, barrier=(2.32638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192615,0.0922558,-0.000149904,1.32158e-07,-4.64161e-11,-74152.4,30.494], Tmin=(100,'K'), Tmax=(783.734,'K')), NASAPolynomial(coeffs=[9.87984,0.0329482,-1.75109e-05,3.47816e-09,-2.45226e-13,-75367.8,-11.9475], Tmin=(783.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-617.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCddCF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Csj(Cs-F1sCdH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][C](F)F(3911)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1768.75],'cm^-1')),
        HinderedRotor(inertia=(0.512049,'amu*angstrom^2'), symmetry=1, barrier=(11.773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54622,0.0360267,-6.1862e-05,5.68999e-08,-2.03356e-11,-34259.7,17.1202], Tmin=(100,'K'), Tmax=(816.891,'K')), NASAPolynomial(coeffs=[5.79763,0.0130514,-6.72067e-06,1.32758e-09,-9.30761e-14,-34555.5,3.53231], Tmin=(816.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'C#CC(F)(F)C(F)[C](F)F(8349)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  C u0 p0 c0 {7,S} {10,T}
10 C u0 p0 c0 {9,T} {12,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-683.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,154,355,414,641,686,1150,1196,190,488,555,1236,1407,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.76025,'amu*angstrom^2'), symmetry=1, barrier=(63.4636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14489,'amu*angstrom^2'), symmetry=1, barrier=(26.3232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.477007,'amu*angstrom^2'), symmetry=1, barrier=(10.9673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.367884,0.104847,-0.000176628,1.50821e-07,-4.9828e-11,-82045.3,29.249], Tmin=(100,'K'), Tmax=(845.384,'K')), NASAPolynomial(coeffs=[13.1051,0.0270381,-1.36201e-05,2.59986e-09,-1.77483e-13,-83820.8,-30.5157], Tmin=(845.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC#CC(F)(F)C(F)[C](F)F(8350)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u0 p0 c0 {8,S} {11,T}
11 C u0 p0 c0 {6,S} {10,T}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-783.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,154,355,414,641,686,1150,1196,190,488,555,1236,1407,2175,525,239,401,1367,180,180,1800.47],'cm^-1')),
        HinderedRotor(inertia=(0.387514,'amu*angstrom^2'), symmetry=1, barrier=(8.90971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38784,'amu*angstrom^2'), symmetry=1, barrier=(8.91721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24101,'amu*angstrom^2'), symmetry=1, barrier=(51.5252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.59261,0.11391,-0.000206375,1.86916e-07,-6.46188e-11,-94040.2,32.3758], Tmin=(100,'K'), Tmax=(848.019,'K')), NASAPolynomial(coeffs=[11.4573,0.0320146,-1.71946e-05,3.34872e-09,-2.3038e-13,-95182.9,-18.4596], Tmin=(848.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = 'FC=[C]C(F)=C[C](F)F(7702)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,D} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {5,D} {9,S}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-336.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,86,203,488,582,605,741,161,297,490,584,780,1358,615,860,1140,1343,3152,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.00425008,'amu*angstrom^2'), symmetry=1, barrier=(4.40065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30238,'amu*angstrom^2'), symmetry=1, barrier=(52.9363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.872645,0.0751001,-0.000103232,7.74178e-08,-2.37407e-11,-40350.6,27.47], Tmin=(100,'K'), Tmax=(790.772,'K')), NASAPolynomial(coeffs=[10.0705,0.0285728,-1.49725e-05,3.0081e-09,-2.1569e-13,-41805.3,-14.745], Tmin=(790.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cd-F1sCd)(Cd-F1sH))"""),
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
    label = 'FC=[C]C(F)=C(F)[C](F)F(8351)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {11,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-488.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,86,203,488,582,605,741,161,297,490,584,780,1358,615,860,1140,1343,3152,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.201901,'amu*angstrom^2'), symmetry=1, barrier=(4.64211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0365344,'amu*angstrom^2'), symmetry=1, barrier=(54.3599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0382603,0.0958782,-0.000161085,1.41353e-07,-4.88298e-11,-58586.7,31.3887], Tmin=(100,'K'), Tmax=(797.057,'K')), NASAPolynomial(coeffs=[11.1567,0.0294167,-1.59409e-05,3.16778e-09,-2.22659e-13,-60020.4,-17.6046], Tmin=(797.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-488.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + radical(Cdj(Cd-F1sCd)(Cd-F1sH))"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'FC=[C]C(F)(F)[C](F)C(F)F(8352)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {5,S} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-792.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.646752,0.113876,-0.000197509,1.77507e-07,-6.21315e-11,-95197,35.7367], Tmin=(100,'K'), Tmax=(814.321,'K')), NASAPolynomial(coeffs=[11.5212,0.0356438,-1.93945e-05,3.84698e-09,-2.69415e-13,-96566.5,-16.7084], Tmin=(814.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-792.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(CsCsCsF1s) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'F[C]=CC(F)(F)C(F)[C](F)F(8353)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-791.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.742273,0.115797,-0.000201054,1.78983e-07,-6.18895e-11,-94978.5,34.3189], Tmin=(100,'K'), Tmax=(820.42,'K')), NASAPolynomial(coeffs=[12.3173,0.0340056,-1.83867e-05,3.63055e-09,-2.53318e-13,-96511.6,-22.385], Tmin=(820.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-791.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=CC(F)(F)[C](F)[C](F)F(8354)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-850.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.669572,0.11689,-0.000212061,1.96147e-07,-6.94802e-11,-102157,35.5763], Tmin=(100,'K'), Tmax=(836.417,'K')), NASAPolynomial(coeffs=[10.1425,0.037617,-2.04604e-05,4.02903e-09,-2.79788e-13,-103002,-8.8902], Tmin=(836.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-850.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(F)C(F)F(8355)',
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
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-733.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,136,307,446,511,682,757,1180,1185,235,523,627,1123,1142,1372,1406,3097,1685,370,167,640,1190,180.793,185.329],'cm^-1')),
        HinderedRotor(inertia=(0.336929,'amu*angstrom^2'), symmetry=1, barrier=(8.2115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12259,'amu*angstrom^2'), symmetry=1, barrier=(29.1392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35694,'amu*angstrom^2'), symmetry=1, barrier=(8.21005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.691044,0.112418,-0.000185088,1.58322e-07,-5.36055e-11,-88019.4,34.3794], Tmin=(100,'K'), Tmax=(776.797,'K')), NASAPolynomial(coeffs=[13.6044,0.0321977,-1.74202e-05,3.47269e-09,-2.44998e-13,-90041,-29.6937], Tmin=(776.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=[C]C(F)(F)[CH]C(F)(F)F(8356)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-863.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.321028,0.105593,-0.000177891,1.60269e-07,-5.71156e-11,-103673,34.8916], Tmin=(100,'K'), Tmax=(786.716,'K')), NASAPolynomial(coeffs=[10.432,0.0369175,-2.02534e-05,4.06252e-09,-2.87537e-13,-104932,-11.652], Tmin=(786.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-863.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'FC=C(F)C(F)(F)[CH][C](F)F(8357)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {7,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-841.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.646067,0.117361,-0.000216213,2.02984e-07,-7.27561e-11,-101101,35.3081], Tmin=(100,'K'), Tmax=(835.487,'K')), NASAPolynomial(coeffs=[9.32059,0.0393154,-2.16432e-05,4.28109e-09,-2.97987e-13,-101708,-4.64907], Tmin=(835.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-841.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][C]=C(F)C(F)C(F)(F)F(8358)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {11,S} {13,S}
11 C u1 p0 c0 {9,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-868.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107311,0.0978094,-0.000136354,9.92389e-08,-2.90243e-11,-104344,31.7706], Tmin=(100,'K'), Tmax=(832.984,'K')), NASAPolynomial(coeffs=[13.7512,0.0312607,-1.65158e-05,3.32858e-09,-2.39127e-13,-106653,-32.5561], Tmin=(832.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-868.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-CsF1s))"""),
)

species(
    label = 'FC=C(F)[C](F)C(F)[C](F)F(8359)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
8  C u1 p0 c0 {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {4,S} {5,S} {7,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-854.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.679207,0.112341,-0.000185314,1.59371e-07,-5.40845e-11,-102579,33.9563], Tmin=(100,'K'), Tmax=(793.184,'K')), NASAPolynomial(coeffs=[13.1809,0.0330515,-1.76043e-05,3.48057e-09,-2.44159e-13,-104482,-27.8366], Tmin=(793.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-854.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCdCsF1s) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(F)C(F)(F)F(8360)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {10,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-787.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.554257,0.107985,-0.000163911,1.26825e-07,-3.88564e-11,-94613.2,32.7346], Tmin=(100,'K'), Tmax=(800.883,'K')), NASAPolynomial(coeffs=[15.3406,0.0285963,-1.52161e-05,3.04621e-09,-2.16889e-13,-97159.1,-40.419], Tmin=(800.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-787.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C(F)[C](F)F(8361)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {11,D}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 C u1 p0 c0 {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-782.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,259,529,569,1128,1321,1390,3140,246,474,533,1155,190,488,555,1236,1407,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.116602,'amu*angstrom^2'), symmetry=1, barrier=(2.6809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116813,'amu*angstrom^2'), symmetry=1, barrier=(2.68577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.666824,'amu*angstrom^2'), symmetry=1, barrier=(15.3316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.862234,0.11796,-0.000203437,1.78509e-07,-6.08038e-11,-94002.9,34.1586], Tmin=(100,'K'), Tmax=(824.287,'K')), NASAPolynomial(coeffs=[13.3233,0.0324105,-1.73461e-05,3.40624e-09,-2.36772e-13,-95773.8,-28.0936], Tmin=(824.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-782.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + group(Cds-CdsHH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    E0 = (-311.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-154.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-185.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (155.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (163.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-302.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-222.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-75.1002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-234.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-154.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-25.0593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-173.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-72.3979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-100.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-0.454435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (140.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-78.3809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-81.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-177.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-107.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-169.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-229.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-167.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-112.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-125.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-162.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-116.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-129.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['CHFCF2(55)', 'FC=C=C(F)F(5948)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C]=CF(7574)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['F[CH]C(=C(F)F)C(F)[C](F)F(7575)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]F(744)', 'F[CH]C(F)(F)[C]=CF(6423)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CF(1054)', 'F[C](F)C(F)[C](F)F(6178)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=C1C(F)(F)C(F)C1(F)F(8345)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=CC(F)(F)C(F)=C(F)F(8346)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'FC=[C]C(F)(F)C=C(F)F(7698)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(69.2597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CHFCF2(55)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.36858e-17,'m^3/(mol*s)'), n=6.25044, Ea=(6.7149,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6646553062382831, var=1.4576930825215244, Tref=1000.0, N=48, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'FC=[C]C(F)(F)C(F)=C(F)F(8347)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(3.41699,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC=C=C(F)C(F)[C](F)F(8348)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(48.9034,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[CH][C](F)F(3911)', 'FC=C=C(F)F(5948)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(5.69886,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CC(F)(F)C(F)[C](F)F(8349)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'FC#CC(F)(F)C(F)[C](F)F(8350)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C](F)F(3911)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(14.7095,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', 'FC=[C]C(F)=C[C](F)F(7702)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(14.9135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'FC=[C]C(F)=C(F)[C](F)F(8351)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(220.211,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CF2(43)', 'F[CH]C(F)(F)[C]=CF(6423)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=[C]C(F)(F)[C](F)C(F)F(8352)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['F[C]=CC(F)(F)C(F)[C](F)F(8353)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=CC(F)(F)[C](F)[C](F)F(8354)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]=[C]C(F)(F)C(F)C(F)F(8355)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=[C]C(F)(F)[CH]C(F)(F)F(8356)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(143.75,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=C(F)C(F)(F)[CH][C](F)F(8357)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(198.893,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['F[CH][C]=C(F)C(F)C(F)(F)F(8358)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.958,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['FC=C(F)[C](F)C(F)[C](F)F(8359)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(148.191,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    products = ['[CH]=[C]C(F)(F)C(F)C(F)(F)F(8360)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(195.086,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(F)C(F)(F)C(F)[C](F)F(8361)'],
    products = ['FC=[C]C(F)(F)C(F)[C](F)F(7577)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2834',
    isomers = [
        'FC=[C]C(F)(F)C(F)[C](F)F(7577)',
    ],
    reactants = [
        ('CHFCF2(55)', 'FC=C=C(F)F(5948)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2834',
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

