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
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,326,540,652,719,1357,194,682,905,1196,1383,3221,180,1575.67],'cm^-1')),
        HinderedRotor(inertia=(0.335659,'amu*angstrom^2'), symmetry=1, barrier=(7.71747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333301,'amu*angstrom^2'), symmetry=1, barrier=(7.66324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79501,0.0835197,-9.01896e-05,9.4965e-09,4.04576e-11,-90394.2,30.3174], Tmin=(100,'K'), Tmax=(501.329,'K')), NASAPolynomial(coeffs=[9.18952,0.0384212,-2.07177e-05,4.16779e-09,-2.97193e-13,-91510.8,-7.12715], Tmin=(501.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-752.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
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
    label = 'FC=C(F)OOCF(3236)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-608.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,326,540,652,719,1357,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.226993,'amu*angstrom^2'), symmetry=1, barrier=(5.21902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226976,'amu*angstrom^2'), symmetry=1, barrier=(5.21861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24262,'amu*angstrom^2'), symmetry=1, barrier=(51.5623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3151.58,'J/mol'), sigma=(5.06004,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.27 K, Pc=55.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06891,0.072491,-0.000116297,1.08471e-07,-4.06571e-11,-73105,25.8181], Tmin=(100,'K'), Tmax=(777.321,'K')), NASAPolynomial(coeffs=[6.04495,0.034313,-1.83642e-05,3.6723e-09,-2.60394e-13,-73498.8,5.50817], Tmin=(777.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-608.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'C=C(F)OOOC(F)=CF(6925)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {4,S} {5,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-438.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,390,483,597,572,732,631,807,1275,1439,2950,3100,1380,975,1025,1650,194,682,905,1196,1383,3221,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858055,0.0802465,-9.26549e-05,3.47325e-08,1.40229e-11,-52596.2,33.1148], Tmin=(100,'K'), Tmax=(528.153,'K')), NASAPolynomial(coeffs=[8.72078,0.0376934,-2.00698e-05,4.03968e-09,-2.8905e-13,-53663.8,-2.04306], Tmin=(528.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CdCFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=C(F)COC(F)C(=O)F(6926)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {5,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-1143.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.531257,0.0831631,-0.000102451,6.753e-08,-1.82852e-11,-137374,29.8517], Tmin=(100,'K'), Tmax=(888.08,'K')), NASAPolynomial(coeffs=[11.7481,0.0326402,-1.7114e-05,3.46772e-09,-2.50881e-13,-139367,-22.9317], Tmin=(888.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1143.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(COCsFO)"""),
)

species(
    label = 'F[C]1COO[C](F)C(F)O1(6927)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
9  C u1 p0 c0 {2,S} {4,S} {8,S}
10 C u1 p0 c0 {3,S} {6,S} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-534.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,2750,3150,900,1100,353,437,368,578,616,798,1350,1522,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0536693,0.0482074,7.2211e-05,-1.65952e-07,7.75878e-11,-64115.2,30.1998], Tmin=(100,'K'), Tmax=(916.117,'K')), NASAPolynomial(coeffs=[37.9956,-0.0165072,1.28815e-05,-2.49439e-09,1.57615e-13,-75303.3,-172.644], Tmin=(916.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-534.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCFHO) + ring(Cycloheptane) + radical(CsCsF1sO2s) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)[CH]OOC(F)=CF(6928)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u1 p0 c0 {4,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,391,562,707,872,1109,1210,1289,3137,3025,407.5,1350,352.5,326,540,652,719,1357,194,682,905,1196,1383,3221,180,2509.27],'cm^-1')),
        HinderedRotor(inertia=(0.22901,'amu*angstrom^2'), symmetry=1, barrier=(5.2654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229162,'amu*angstrom^2'), symmetry=1, barrier=(5.2689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41421,'amu*angstrom^2'), symmetry=1, barrier=(55.5074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41894,'amu*angstrom^2'), symmetry=1, barrier=(55.6162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.251199,0.106675,-0.000188082,1.75939e-07,-6.36058e-11,-48494.2,34.8538], Tmin=(100,'K'), Tmax=(825.107,'K')), NASAPolynomial(coeffs=[8.13688,0.0406162,-2.18244e-05,4.30492e-09,-3.00383e-13,-49014,1.237], Tmin=(825.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-CsF1sH)) + radical(CCsJOOC)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3651.34,'J/mol'), sigma=(5.76317,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.33 K, Pc=43.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.389635,0.105786,-0.000169482,1.46908e-07,-5.11987e-11,-73229.4,33.593], Tmin=(100,'K'), Tmax=(754.463,'K')), NASAPolynomial(coeffs=[11.4958,0.0366597,-1.98958e-05,3.99074e-09,-2.83406e-13,-74848.9,-19.2456], Tmin=(754.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(CsCsF1sO2s) + radical(OCJC=O)"""),
)

species(
    label = '[O]C(F)COOC(F)=[C]F(6929)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-321.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,293,496,537,1218,167,640,1190,180,180,1711.68],'cm^-1')),
        HinderedRotor(inertia=(0.21495,'amu*angstrom^2'), symmetry=1, barrier=(4.94212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21388,'amu*angstrom^2'), symmetry=1, barrier=(4.91751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215102,'amu*angstrom^2'), symmetry=1, barrier=(4.94562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33307,'amu*angstrom^2'), symmetry=1, barrier=(53.6419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.229308,0.106758,-0.000190439,1.79885e-07,-6.54755e-11,-38500.1,35.5164], Tmin=(100,'K'), Tmax=(825.683,'K')), NASAPolynomial(coeffs=[7.67433,0.0413587,-2.23775e-05,4.42394e-09,-3.08961e-13,-38881.2,4.49603], Tmin=(825.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-CsF1sH)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O[C](F)COOC(F)=[C]F(6930)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-355.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,293,496,537,1218,167,640,1190,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.731726,0.12038,-0.000226463,2.13384e-07,-7.61422e-11,-42588.8,36.8868], Tmin=(100,'K'), Tmax=(846.219,'K')), NASAPolynomial(coeffs=[9.34787,0.0386822,-2.12844e-05,4.18628e-09,-2.89578e-13,-43075.4,-2.85436], Tmin=(846.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-O2sHH)(F1s)(O2s-H)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3010,987.5,1337.5,450,1655,262,390,483,597,572,732,631,807,1275,1439,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.31982,'amu*angstrom^2'), symmetry=1, barrier=(7.35329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625962,'amu*angstrom^2'), symmetry=1, barrier=(14.3921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321569,'amu*angstrom^2'), symmetry=1, barrier=(7.39351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08416,'amu*angstrom^2'), symmetry=1, barrier=(47.919,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22587,0.101483,-0.000161317,1.3789e-07,-4.73622e-11,-78770.8,33.2537], Tmin=(100,'K'), Tmax=(748.709,'K')), NASAPolynomial(coeffs=[11.7116,0.033882,-1.82174e-05,3.64776e-09,-2.58943e-13,-80451.2,-20.1667], Tmin=(748.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
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
    label = 'O=[C]COOC(F)=CF(6931)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,S} {8,D}
8  C u0 p0 c0 {2,S} {7,D} {12,S}
9  C u1 p0 c0 {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-343.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,326,540,652,719,1357,194,682,905,1196,1383,3221,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.366316,'amu*angstrom^2'), symmetry=1, barrier=(8.42233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76105,'amu*angstrom^2'), symmetry=1, barrier=(40.4901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00302348,'amu*angstrom^2'), symmetry=1, barrier=(8.38955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75896,'amu*angstrom^2'), symmetry=1, barrier=(40.4419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37069,0.070753,-5.7709e-05,-4.89708e-08,8.64431e-11,-41264.3,28.0587], Tmin=(100,'K'), Tmax=(464.267,'K')), NASAPolynomial(coeffs=[8.17325,0.0354737,-1.91013e-05,3.82727e-09,-2.71663e-13,-42147.4,-2.24793], Tmin=(464.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C(F)COO[C]=CF(6932)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-326.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,615,860,1140,1343,3152,1685,370,180,180,1303.71],'cm^-1')),
        HinderedRotor(inertia=(1.89683,'amu*angstrom^2'), symmetry=1, barrier=(43.6119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89718,'amu*angstrom^2'), symmetry=1, barrier=(43.6199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951805,'amu*angstrom^2'), symmetry=1, barrier=(21.8839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.430908,'amu*angstrom^2'), symmetry=1, barrier=(9.90742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671306,0.0805574,-0.000108944,8.12209e-08,-2.49927e-11,-39140.9,29.3269], Tmin=(100,'K'), Tmax=(784.998,'K')), NASAPolynomial(coeffs=[10.0654,0.0326904,-1.74797e-05,3.54578e-09,-2.55801e-13,-40615.8,-13.7201], Tmin=(784.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-F1sH)(O2s-O2s))"""),
)

species(
    label = '[CH]=C(F)OOCC(=O)F(6933)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-323.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,293,496,537,1218,3120,650,792.5,1650,180,1403.12],'cm^-1')),
        HinderedRotor(inertia=(0.335435,'amu*angstrom^2'), symmetry=1, barrier=(7.71232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335395,'amu*angstrom^2'), symmetry=1, barrier=(7.71139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68702,'amu*angstrom^2'), symmetry=1, barrier=(38.7879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68699,'amu*angstrom^2'), symmetry=1, barrier=(38.7872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.704425,0.081827,-0.000111333,7.49088e-08,-1.57393e-11,-38751.2,29.8819], Tmin=(100,'K'), Tmax=(585.031,'K')), NASAPolynomial(coeffs=[9.53365,0.0334522,-1.8051e-05,3.64956e-09,-2.61757e-13,-39989.5,-9.73482], Tmin=(585.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'O=C(F)[CH]F(215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-442.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(1.5365,'amu*angstrom^2'), symmetry=1, barrier=(35.3272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91969,0.0052474,7.61145e-05,-1.81228e-07,1.26791e-10,-53179.7,10.2212], Tmin=(10,'K'), Tmax=(489.215,'K')), NASAPolynomial(coeffs=[4.0714,0.0192239,-1.3397e-05,4.3331e-09,-5.27008e-13,-53376.6,7.73666], Tmin=(489.215,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-442.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC(=O)F(1554)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-363.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,180,475.353],'cm^-1')),
        HinderedRotor(inertia=(0.0192623,'amu*angstrom^2'), symmetry=1, barrier=(3.18239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86108,0.0121233,3.70793e-05,-8.49933e-08,5.10894e-11,-43750.9,11.7734], Tmin=(10,'K'), Tmax=(581.186,'K')), NASAPolynomial(coeffs=[3.98792,0.0215968,-1.40743e-05,4.31471e-09,-5.0294e-13,-43940.4,9.72707], Tmin=(581.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-363.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)=CF(214)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-249.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.358734,'amu*angstrom^2'), symmetry=1, barrier=(8.24799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3441.78,'J/mol'), sigma=(5.35945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.60 K, Pc=50.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75962,0.0203814,3.21936e-05,-1.12324e-07,8.28e-11,-30007.8,12.0802], Tmin=(10,'K'), Tmax=(521.492,'K')), NASAPolynomial(coeffs=[5.63612,0.021414,-1.5147e-05,4.91851e-09,-5.9756e-13,-30413.3,2.2378], Tmin=(521.492,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-249.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CFO(79)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-261.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,611,648,830,1210,1753],'cm^-1')),
        HinderedRotor(inertia=(0.112588,'amu*angstrom^2'), symmetry=1, barrier=(22.3442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95484,0.00266306,5.93929e-05,-1.151e-07,6.74789e-11,-31406.5,8.78105], Tmin=(10,'K'), Tmax=(558.359,'K')), NASAPolynomial(coeffs=[3.00171,0.0196207,-1.33752e-05,4.27391e-09,-5.17274e-13,-31457.9,11.4099], Tmin=(558.359,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-261.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[CH2]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHFCF[Z](72)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-58.4765,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([18.6764,93.2568,111.933],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([203.406,408.59,750.559,801.84,1025.92,1150.34,1361.98,1782.01,3238.53],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9644,0.00217073,4.76141e-05,-9.63411e-08,5.92621e-11,-7031.25,9.1333], Tmin=(10,'K'), Tmax=(529.917,'K')), NASAPolynomial(coeffs=[3.26037,0.0150322,-1.01557e-05,3.21347e-09,-3.84694e-13,-7062.61,11.0829], Tmin=(529.917,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-58.4765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""F[C]DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OCC(=O)F(2461)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-365.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.431789,'amu*angstrom^2'), symmetry=1, barrier=(9.92767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0958417,'amu*angstrom^2'), symmetry=1, barrier=(54.572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81447,0.024853,-7.75773e-06,-6.46438e-09,3.7492e-12,-43928,11.2993], Tmin=(10,'K'), Tmax=(1043.58,'K')), NASAPolynomial(coeffs=[7.79384,0.0179742,-9.90682e-06,2.59764e-09,-2.63701e-13,-45214.6,-10.2533], Tmin=(1043.58,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-365.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""[O]OCC(DO)F""", comment="""Thermo library: 2-BTP_G4"""),
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
    label = '[CH2]OOC(F)=CF(3295)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-200.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,194,682,905,1196,1383,3221,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.275753,'amu*angstrom^2'), symmetry=1, barrier=(6.3401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27446,'amu*angstrom^2'), symmetry=1, barrier=(52.2943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0152951,'amu*angstrom^2'), symmetry=1, barrier=(52.3822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40076,0.0641754,-0.000102073,9.45495e-08,-3.51759e-11,-24052.3,22.72], Tmin=(100,'K'), Tmax=(783.277,'K')), NASAPolynomial(coeffs=[5.83685,0.0303891,-1.60524e-05,3.19079e-09,-2.2544e-13,-24405.7,4.58193], Tmin=(783.277,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsJOOC)"""),
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
    label = 'O=C(F)COOC(F)=[C]F(6934)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-482.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,293,496,537,1218,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61565,'amu*angstrom^2'), symmetry=1, barrier=(37.1471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346848,'amu*angstrom^2'), symmetry=1, barrier=(7.97471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00845218,'amu*angstrom^2'), symmetry=1, barrier=(37.115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00845585,'amu*angstrom^2'), symmetry=1, barrier=(37.1198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.084484,0.0970966,-0.000166997,1.55187e-07,-5.67255e-11,-57947,33.7472], Tmin=(100,'K'), Tmax=(794.499,'K')), NASAPolynomial(coeffs=[8.36545,0.0375277,-2.07798e-05,4.17401e-09,-2.95385e-13,-58698.7,-0.747789], Tmin=(794.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O=C(F)CO[O+]=[C-]C(F)F(6935)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p1 c+1 {4,S} {10,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u0 p1 c-1 {5,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-629.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13903,0.105338,-0.00018809,1.79936e-07,-6.61952e-11,-75615.6,44.2099], Tmin=(100,'K'), Tmax=(827.604,'K')), NASAPolynomial(coeffs=[6.51365,0.043884,-2.36012e-05,4.65546e-09,-3.24777e-13,-75713.3,19.4356], Tmin=(827.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-629.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFFH) + group(COCsFO) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(F)(F)OOCC(=O)F(6936)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u0 p1 c0 {8,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-552.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2749,0.103284,-0.000158574,1.30125e-07,-4.32149e-11,-66315.6,31.6753], Tmin=(100,'K'), Tmax=(734.507,'K')), NASAPolynomial(coeffs=[12.2209,0.035236,-1.96112e-05,4.00105e-09,-2.87988e-13,-68151.3,-24.7547], Tmin=(734.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-552.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFFO) + group(COCsFO) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C(F)COOC(F)[C]F(6937)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p1 c0 {3,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-591.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,486,617,768,1157,1926,617,898,1187,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.26427,'amu*angstrom^2'), symmetry=1, barrier=(29.068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26427,'amu*angstrom^2'), symmetry=1, barrier=(29.0681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26427,'amu*angstrom^2'), symmetry=1, barrier=(29.068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26427,'amu*angstrom^2'), symmetry=1, barrier=(29.0681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32271,'amu*angstrom^2'), symmetry=1, barrier=(30.4117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279371,0.0895108,-0.000112991,7.51028e-08,-2.04081e-11,-70992.9,31.5998], Tmin=(100,'K'), Tmax=(886.42,'K')), NASAPolynomial(coeffs=[12.7146,0.0333977,-1.80386e-05,3.69198e-09,-2.68424e-13,-73197.5,-26.8942], Tmin=(886.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-591.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + group(COCsFO) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C=COOC(F)=CF(6938)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-316.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,194,682,905,1196,1383,3221,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.330755,'amu*angstrom^2'), symmetry=1, barrier=(7.60471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331792,'amu*angstrom^2'), symmetry=1, barrier=(7.62855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74675,'amu*angstrom^2'), symmetry=1, barrier=(40.1613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504755,0.0841757,-0.000136517,1.18333e-07,-4.07786e-11,-37917,29.4175], Tmin=(100,'K'), Tmax=(784.691,'K')), NASAPolynomial(coeffs=[10.1235,0.0277979,-1.47039e-05,2.91141e-09,-2.04817e-13,-39200.4,-13.214], Tmin=(784.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d)"""),
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
    label = 'C#COOCC(=O)F(6939)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-168.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.41761,'amu*angstrom^2'), symmetry=1, barrier=(32.5936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41635,'amu*angstrom^2'), symmetry=1, barrier=(32.5648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42027,'amu*angstrom^2'), symmetry=1, barrier=(32.6549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41915,'amu*angstrom^2'), symmetry=1, barrier=(32.6291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02555,0.0702279,-8.42228e-05,5.2448e-08,-1.32638e-11,-20203.3,24.7832], Tmin=(100,'K'), Tmax=(952.413,'K')), NASAPolynomial(coeffs=[11.7802,0.0250587,-1.3082e-05,2.64996e-09,-1.91908e-13,-22251.8,-26.5772], Tmin=(952.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'O=C(F)COOC#CF(6940)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {2,S} {8,T}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-268.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,2175,525,239,401,1367,180,1309.32],'cm^-1')),
        HinderedRotor(inertia=(1.36253,'amu*angstrom^2'), symmetry=1, barrier=(31.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63806,'amu*angstrom^2'), symmetry=1, barrier=(14.6703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35901,'amu*angstrom^2'), symmetry=1, barrier=(31.2464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6659,0.0808975,-0.000119739,9.64926e-08,-3.18011e-11,-32192.4,28.3932], Tmin=(100,'K'), Tmax=(737.215,'K')), NASAPolynomial(coeffs=[9.98709,0.0303242,-1.68422e-05,3.44615e-09,-2.48958e-13,-33566.8,-13.7347], Tmin=(737.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Ct-CtOs) + group(CtCF)"""),
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
    label = 'O=C(F)CO[O+]=[C-]F(2464)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p1 c+1 {3,S} {8,D}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p1 c-1 {2,S} {4,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-643.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52132,0.060256,-7.86318e-05,5.83454e-08,-1.81584e-11,-77310.8,24.0988], Tmin=(100,'K'), Tmax=(771.887,'K')), NASAPolynomial(coeffs=[7.89168,0.0272436,-1.4478e-05,2.93592e-09,-2.11966e-13,-78294.3,-4.98496], Tmin=(771.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CJ2_singlet-FO)"""),
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
    E0 = (-390.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-19.2243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1.83559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-334.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-199.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-46.7474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-231.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-235.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (39.0316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (4.9065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-144.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (64.4475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (81.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (85.1055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-400.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-171.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-88.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-55.7402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-66.8183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (64.2442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-132.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-85.2513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-169.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-68.5189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (157.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (57.2195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-157.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)COOC(F)=CF(6220)'],
    products = ['O=CC(=O)F(335)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(26.2384,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'FC=C(F)OOCF(3236)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(372.836,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(F)OOOC(F)=CF(6925)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(104.673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C(F)COOC(F)=CF(6220)'],
    products = ['O=C(F)COC(F)C(=O)F(6926)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(82.5351,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]1COO[C](F)C(F)O1(6927)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)[CH]OOC(F)=CF(6928)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)[CH]OOC(F)[CH]F(6515)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)[CH]OO[C](F)CF(6213)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)COOC(F)=[C]F(6929)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O[C](F)COOC(F)=[C]F(6930)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC(F)=COOC(F)=CF(6521)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(175.89,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'O=[C]COOC(F)=CF(6931)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=C(F)COO[C]=CF(6932)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[CH]=C(F)OOCC(=O)F(6933)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.11549e+09,'m^3/(mol*s)'), n=-0.68237, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.29431919638206655, var=1.0853977775937997, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C(F)[CH]F(215)', '[O]CC(=O)F(1554)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(70.0634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC(F)=CF(214)', 'CH2CFO(79)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(3.79219,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHFCF[Z](72)', '[O]OCC(=O)F(2461)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(51)', '[CH2]OOC(F)=CF(3295)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -6.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C(F)[CH]OOC(F)=CF(6513)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.43549e+07,'m^3/(mol*s)'), n=0.0910533, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0028590110630623746, var=4.845846562261507, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'O=C(F)COOC(F)=[C]F(6934)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_N-Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_N-Sp-4C-2C_Ext-4C-R
Ea raised from -14.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C(F)CO[O+]=[C-]C(F)F(6935)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.82066e+10,'s^-1'), n=0.827, Ea=(162.407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(F)(F)OOCC(=O)F(6936)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.82066e+10,'s^-1'), n=0.827, Ea=(132.014,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(F)COOC(F)[C]F(6937)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.33334e+12,'s^-1'), n=-1.40567e-07, Ea=(86.0109,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=C=COOC(F)=CF(6938)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(193.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F2(78)', 'C#COOCC(=O)F(6939)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'O=C(F)COOC#CF(6940)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(271.613,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CHF(40)', 'O=C(F)CO[O+]=[C-]F(2464)'],
    products = ['O=C(F)COOC(F)=CF(6220)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.1e+24,'cm^3/(mol*s)'), n=-3.8, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for CF
Exact match found for rate rule [CF]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #1842',
    isomers = [
        'O=C(F)COOC(F)=CF(6220)',
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
    label = 'PDepNetwork #1842',
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

