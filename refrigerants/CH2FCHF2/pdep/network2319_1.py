species(
    label = 'O=[C]CC([C](F)F)C(F)F(8661)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {14,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-721.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,1855,455,950,340.905,340.905,340.908],'cm^-1')),
        HinderedRotor(inertia=(0.0980926,'amu*angstrom^2'), symmetry=1, barrier=(8.09108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980991,'amu*angstrom^2'), symmetry=1, barrier=(8.09102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981082,'amu*angstrom^2'), symmetry=1, barrier=(8.09101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.098101,'amu*angstrom^2'), symmetry=1, barrier=(8.09087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.308506,0.103691,-0.000166617,1.43665e-07,-4.89064e-11,-86648.8,34.7938], Tmin=(100,'K'), Tmax=(813.272,'K')), NASAPolynomial(coeffs=[11.2667,0.0348124,-1.75406e-05,3.39853e-09,-2.35809e-13,-88136.4,-16.2283], Tmin=(813.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-721.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(CCCJ=O)"""),
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
    label = 'O=[C]CC[C](F)F(6187)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-276.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,190,488,555,1236,1407,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.17743,'amu*angstrom^2'), symmetry=1, barrier=(4.07946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177246,'amu*angstrom^2'), symmetry=1, barrier=(4.07524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177823,'amu*angstrom^2'), symmetry=1, barrier=(4.08851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3294.32,'J/mol'), sigma=(5.47565,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.57 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24589,0.0686908,-0.000113564,1.04502e-07,-3.73809e-11,-33217.2,25.4254], Tmin=(100,'K'), Tmax=(839.14,'K')), NASAPolynomial(coeffs=[6.01416,0.0296705,-1.46926e-05,2.81718e-09,-1.93804e-13,-33443.8,6.67503], Tmin=(839.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(CCCJ=O)"""),
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
    label = 'O=[C]CC(F)[C](F)F(7863)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-468.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.276161,'amu*angstrom^2'), symmetry=1, barrier=(6.34948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274614,'amu*angstrom^2'), symmetry=1, barrier=(6.31393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275592,'amu*angstrom^2'), symmetry=1, barrier=(6.33641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.47,'J/mol'), sigma=(5.394,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.40 K, Pc=47.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.713632,0.0822936,-0.000144918,1.33394e-07,-4.70241e-11,-56239,27.6427], Tmin=(100,'K'), Tmax=(848.531,'K')), NASAPolynomial(coeffs=[7.62081,0.0292651,-1.49936e-05,2.88885e-09,-1.98261e-13,-56674.4,-0.20403], Tmin=(848.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC(F)(F)[CH]C(F)F(8660)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-737.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.385676,0.107693,-0.000182207,1.63833e-07,-5.74391e-11,-88602.2,34.9153], Tmin=(100,'K'), Tmax=(823.858,'K')), NASAPolynomial(coeffs=[10.0725,0.0376396,-1.95645e-05,3.82149e-09,-2.65579e-13,-89671.3,-9.54284], Tmin=(823.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[CH]C(F)(F)C(F)F(8762)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {13,S}
9  C u1 p0 c0 {6,S} {7,S} {14,S}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-713.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,3025,407.5,1350,352.5,1855,455,950,333.197,333.324,1939.1],'cm^-1')),
        HinderedRotor(inertia=(0.111679,'amu*angstrom^2'), symmetry=1, barrier=(8.81339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112277,'amu*angstrom^2'), symmetry=1, barrier=(8.82295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112967,'amu*angstrom^2'), symmetry=1, barrier=(8.81427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112356,'amu*angstrom^2'), symmetry=1, barrier=(8.83526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.291205,0.102743,-0.000162989,1.38785e-07,-4.68149e-11,-85696.1,34.8424], Tmin=(100,'K'), Tmax=(806.941,'K')), NASAPolynomial(coeffs=[11.6457,0.0338861,-1.69884e-05,3.2902e-09,-2.28472e-13,-87307.2,-18.2318], Tmin=(806.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-713.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C([C](F)F)C(F)F(7423)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {8,D} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-743.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,374.813,374.914,374.932,374.943],'cm^-1')),
        HinderedRotor(inertia=(0.118405,'amu*angstrom^2'), symmetry=1, barrier=(11.8111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118439,'amu*angstrom^2'), symmetry=1, barrier=(11.8119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118425,'amu*angstrom^2'), symmetry=1, barrier=(11.8119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232459,0.0970303,-0.000127125,8.50037e-08,-2.24778e-11,-89256.6,33.2258], Tmin=(100,'K'), Tmax=(925.568,'K')), NASAPolynomial(coeffs=[16.1358,0.0262894,-1.24765e-05,2.42139e-09,-1.71059e-13,-92286.4,-44.4749], Tmin=(925.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-743.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F[C]F(156)',
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
    label = 'O=[C]C[CH]C(F)F(6299)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-277.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.203319,'amu*angstrom^2'), symmetry=1, barrier=(4.6747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204828,'amu*angstrom^2'), symmetry=1, barrier=(4.7094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204095,'amu*angstrom^2'), symmetry=1, barrier=(4.69254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49715,0.0592632,-7.66315e-05,5.68602e-08,-1.74325e-11,-33299.1,25.6696], Tmin=(100,'K'), Tmax=(789.372,'K')), NASAPolynomial(coeffs=[8.20032,0.0252943,-1.20791e-05,2.33942e-09,-1.64451e-13,-34357.3,-5.08352], Tmin=(789.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=O(514)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.49],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.0892,0.0020037,-1.61643e-05,2.55031e-08,-1.16411e-11,52802.7,4.52493], Tmin=(100,'K'), Tmax=(856.126,'K')), NASAPolynomial(coeffs=[0.961549,0.00569059,-3.48052e-06,7.19222e-10,-5.08058e-14,53738.7,21.4667], Tmin=(856.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([C](F)F)C(F)F(7216)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-577.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.173692,'amu*angstrom^2'), symmetry=1, barrier=(3.99353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174177,'amu*angstrom^2'), symmetry=1, barrier=(4.00467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281333,'amu*angstrom^2'), symmetry=1, barrier=(6.46839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.489151,0.0864055,-0.000145137,1.29895e-07,-4.49194e-11,-69382.3,27.5062], Tmin=(100,'K'), Tmax=(847.692,'K')), NASAPolynomial(coeffs=[8.60787,0.0304288,-1.4824e-05,2.82461e-09,-1.93248e-13,-70124,-6.57641], Tmin=(847.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-577.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Isobutyl)"""),
)

species(
    label = 'O=C1CC(C(F)F)C1(F)F(8763)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {14,S}
10 C u0 p0 c0 {5,D} {7,S} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-953.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587388,0.0698829,-5.75393e-05,2.33892e-08,-3.81526e-12,-114507,28.2193], Tmin=(100,'K'), Tmax=(1445.97,'K')), NASAPolynomial(coeffs=[16.5882,0.0256194,-1.16215e-05,2.21862e-09,-1.54955e-13,-119134,-54.8762], Tmin=(1445.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-953.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-OdCsCs) + longDistanceInteraction_cyclic(Cs(F)2-CO) + ring(Cyclobutanone)"""),
)

species(
    label = 'O=CCC(=C(F)F)C(F)F(8764)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {10,D}
9  C u0 p0 c0 {5,D} {6,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-918.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18035,0.0707258,-6.16303e-05,2.71023e-08,-5.08041e-12,-110434,27.1729], Tmin=(100,'K'), Tmax=(1199.19,'K')), NASAPolynomial(coeffs=[11.0546,0.0377894,-2.04322e-05,4.19904e-09,-3.05704e-13,-112802,-22.2585], Tmin=(1199.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-918.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF)"""),
)

species(
    label = 'O=C=CC(C(F)F)C(F)F(8765)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {13,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
9  C u0 p0 c0 {6,S} {10,D} {14,S}
10 C u0 p0 c0 {5,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-955.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0500979,0.0934727,-0.000123422,8.70263e-08,-2.48345e-11,-114791,29.6046], Tmin=(100,'K'), Tmax=(851.166,'K')), NASAPolynomial(coeffs=[12.9604,0.0327987,-1.64921e-05,3.27123e-09,-2.33386e-13,-116988,-30.5994], Tmin=(851.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-955.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
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
    label = '[CH2][C]=O(4373)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'CHF2(82)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]CC=C(F)F(6370)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,D}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,1855,455,950,365.082],'cm^-1')),
        HinderedRotor(inertia=(0.00126481,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18373,'amu*angstrom^2'), symmetry=1, barrier=(17.3831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83747,0.0433862,-3.22194e-05,1.08234e-08,-1.422e-12,-39127.7,22.1243], Tmin=(100,'K'), Tmax=(1780.4,'K')), NASAPolynomial(coeffs=[15.1504,0.0134763,-7.0201e-06,1.38754e-09,-9.70449e-14,-43868.2,-49.7823], Tmin=(1780.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(CCCJ=O)"""),
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
    label = 'O=[C]CC(=C(F)F)C(F)F(8766)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 C u1 p0 c0 {5,D} {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-759.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,182,240,577,636,1210,1413,1855,455,950,495.342,495.349],'cm^-1')),
        HinderedRotor(inertia=(0.000687057,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737058,'amu*angstrom^2'), symmetry=1, barrier=(12.8344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737081,'amu*angstrom^2'), symmetry=1, barrier=(12.8343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943524,0.0753542,-8.14252e-05,4.64122e-08,-1.11292e-11,-91185.8,28.8616], Tmin=(100,'K'), Tmax=(979.283,'K')), NASAPolynomial(coeffs=[10.9116,0.034638,-1.90586e-05,3.95463e-09,-2.90186e-13,-93138.1,-19.0201], Tmin=(979.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-759.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(CCCJ=O)"""),
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
    label = 'O=C=CC([C](F)F)C(F)F(8767)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {13,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {5,D} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-753.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,2120,512.5,787.5,180.164,182.934,192.573],'cm^-1')),
        HinderedRotor(inertia=(0.433706,'amu*angstrom^2'), symmetry=1, barrier=(11.1198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815673,'amu*angstrom^2'), symmetry=1, barrier=(22.1858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42316,'amu*angstrom^2'), symmetry=1, barrier=(11.1234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0630951,0.0965849,-0.000138299,9.4387e-08,-2.02648e-11,-90458,31.2174], Tmin=(100,'K'), Tmax=(609.903,'K')), NASAPolynomial(coeffs=[12.4342,0.0311872,-1.61635e-05,3.1908e-09,-2.24846e-13,-92259.7,-24.7486], Tmin=(609.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-753.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CsCsF1sF1s)"""),
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
    label = 'O=[C]CC([CH]F)=C(F)F(8768)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {1,S} {6,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u1 p0 c0 {4,D} {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-394.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.00578554,'amu*angstrom^2'), symmetry=1, barrier=(4.1537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0284966,'amu*angstrom^2'), symmetry=1, barrier=(20.4584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120313,'amu*angstrom^2'), symmetry=1, barrier=(20.4635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28933,0.0635928,-5.83903e-05,2.63602e-08,-4.88688e-12,-47387.8,27.6334], Tmin=(100,'K'), Tmax=(1252.16,'K')), NASAPolynomial(coeffs=[12.7728,0.0269091,-1.44457e-05,2.96334e-09,-2.15575e-13,-50263.6,-30.3498], Tmin=(1252.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[C](C(F)F)C(F)F(8769)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u1 p0 c0 {6,S} {7,S} {8,S}
10 C u1 p0 c0 {5,D} {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-735.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.329806,0.110809,-0.000202014,1.91518e-07,-6.88055e-11,-88286.1,33.8492], Tmin=(100,'K'), Tmax=(853.336,'K')), NASAPolynomial(coeffs=[7.05892,0.0422588,-2.1898e-05,4.22692e-09,-2.89884e-13,-88312.2,6.61016], Tmin=(853.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-735.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CC([C](F)F)C(F)F(8770)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {13,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {5,S} {8,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-733.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383757,0.0930966,-0.000105817,5.87663e-08,-1.27215e-11,-88111.5,33.4693], Tmin=(100,'K'), Tmax=(1131.65,'K')), NASAPolynomial(coeffs=[20.2766,0.0200696,-9.01997e-06,1.74273e-09,-1.24097e-13,-92787.5,-68.7602], Tmin=(1131.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'O=[C][CH]C(C(F)F)C(F)F(8771)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {13,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
9  C u1 p0 c0 {6,S} {10,S} {14,S}
10 C u1 p0 c0 {5,D} {9,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-753.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0265525,0.0940952,-0.000126018,8.96581e-08,-2.55939e-11,-90439,31.405], Tmin=(100,'K'), Tmax=(854.568,'K')), NASAPolynomial(coeffs=[13.514,0.0307135,-1.47622e-05,2.86223e-09,-2.01364e-13,-92753.2,-31.7919], Tmin=(854.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-753.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC[C]([C](F)F)C(F)F(8772)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {5,D} {6,S} {14,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-696.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,522,611,926,1093,1137,1374,1416,3112,360,370,350,2782.5,750,1395,475,1775,1000,190,488,555,1236,1407,180,180,2098.14],'cm^-1')),
        HinderedRotor(inertia=(0.29985,'amu*angstrom^2'), symmetry=1, barrier=(6.89414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299663,'amu*angstrom^2'), symmetry=1, barrier=(6.88984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297916,'amu*angstrom^2'), symmetry=1, barrier=(6.84967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299682,'amu*angstrom^2'), symmetry=1, barrier=(6.89028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1802,0.108567,-0.000200241,1.94263e-07,-7.12002e-11,-83598.5,34.6458], Tmin=(100,'K'), Tmax=(849.646,'K')), NASAPolynomial(coeffs=[5.3098,0.0457076,-2.392e-05,4.64252e-09,-3.1966e-13,-83195.4,16.9164], Tmin=(849.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-696.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'O=CCC([C](F)F)[C](F)F(8773)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {5,D} {7,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-682.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,146,234,414,562,504,606,1176,1296,1354,1460,2782.5,750,1395,475,1775,1000,291.369,291.369,2798.62],'cm^-1')),
        HinderedRotor(inertia=(0.129425,'amu*angstrom^2'), symmetry=1, barrier=(7.79707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129425,'amu*angstrom^2'), symmetry=1, barrier=(7.79707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129425,'amu*angstrom^2'), symmetry=1, barrier=(7.79707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.596884,'amu*angstrom^2'), symmetry=1, barrier=(35.9587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.16232,0.101497,-0.000165048,1.46749e-07,-5.14882e-11,-81961,34.2159], Tmin=(100,'K'), Tmax=(812.915,'K')), NASAPolynomial(coeffs=[9.51567,0.0382641,-1.95642e-05,3.81447e-09,-2.65612e-13,-83018.7,-7.29746], Tmin=(812.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-682.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'O=[C]CC([CH]F)C(F)(F)F(8774)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
9  C u1 p0 c0 {4,S} {6,S} {14,S}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-753.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.188851,0.101173,-0.000161858,1.40791e-07,-4.84737e-11,-90446.5,33.4515], Tmin=(100,'K'), Tmax=(810.778,'K')), NASAPolynomial(coeffs=[10.4414,0.0361351,-1.82358e-05,3.53975e-09,-2.46022e-13,-91756.4,-13.0508], Tmin=(810.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-753.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)CC([CH]F)[C](F)F(8775)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {14,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {1,S} {5,D} {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-721.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,190,488,555,1236,1407,486,617,768,1157,1926,188.446,188.463,188.475,2085.52],'cm^-1')),
        HinderedRotor(inertia=(0.223742,'amu*angstrom^2'), symmetry=1, barrier=(5.63946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223736,'amu*angstrom^2'), symmetry=1, barrier=(5.63948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00474603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63662,'amu*angstrom^2'), symmetry=1, barrier=(41.254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.188512,0.103386,-0.000172689,1.57233e-07,-5.60621e-11,-86601,34.4027], Tmin=(100,'K'), Tmax=(819.203,'K')), NASAPolynomial(coeffs=[8.64352,0.0403777,-2.09112e-05,4.08845e-09,-2.84713e-13,-87380.9,-2.37336], Tmin=(819.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-721.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(CsCsFFH) + group(COCsFO) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
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
    E0 = (-301.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (64.7428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (229.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-141.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-136.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-56.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (175.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (281.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-293.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-238.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-238.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-243.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-190.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-149.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-124.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-185.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-120.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (29.7303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-50.7652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-61.4804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-127.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-142.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-179.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-171.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-189.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-89.8396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-59.2008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['CH2CO(28)', 'FC(F)=CC(F)F(2834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'O=[C]CC[C](F)F(6187)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33604e-10,'m^3/(mol*s)'), n=4.72997, Ea=(125.584,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'O=[C]CC(F)[C](F)F(7863)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.96165e-06,'m^3/(mol*s)'), n=3.30609, Ea=(138.972,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s',), comment="""Estimated from node CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=[C]CC(F)(F)[CH]C(F)F(8660)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C]C[CH]C(F)(F)C(F)F(8762)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['C=C([O])C([C](F)F)C(F)F(7423)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]F(156)', 'O=[C]C[CH]C(F)F(6299)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]=O(514)', '[CH2]C([C](F)F)C(F)F(7216)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=C1CC(C(F)F)C1(F)F(8763)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=CCC(=C(F)F)C(F)F(8764)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=C=CC(C(F)F)C(F)F(8765)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CO(13)', '[CH2]C([C](F)F)C(F)F(7216)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(33.2011,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=O(4373)', 'FC(F)=CC(F)F(2834)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(22.7438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHF2(82)', 'O=[C]CC=C(F)F(6370)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(13.0602,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'O=[C]CC(=C(F)F)C(F)F(8766)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(168,'m^3/(mol*s)'), n=1.64, Ea=(3.34727,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C_Ext-7R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2CO(28)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(9.69925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'O=C=CC([C](F)F)C(F)F(8767)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(0.741057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=O(4373)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(4.12539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'O=[C]CC([CH]F)=C(F)F(8768)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(205.321,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CF2(43)', 'O=[C]C[CH]C(F)F(6299)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=[C]C[C](C(F)F)C(F)F(8769)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['[O]C=CC([C](F)F)C(F)F(8770)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=[C][CH]C(C(F)F)C(F)F(8771)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.01948e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=CC[C]([C](F)F)C(F)F(8772)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.98524e+06,'s^-1'), n=1.75122, Ea=(105.056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=CCC([C](F)F)[C](F)F(8773)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(893.101,'s^-1'), n=2.5049, Ea=(73.5615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;XH_out] for rate rule [R4H_SSS;C_rad_out_noH;CO_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    products = ['O=[C]CC([CH]F)C(F)(F)F(8774)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(212.004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C(F)CC([CH]F)[C](F)F(8775)'],
    products = ['O=[C]CC([C](F)F)C(F)F(8661)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(242.191,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

network(
    label = 'PDepNetwork #2319',
    isomers = [
        'O=[C]CC([C](F)F)C(F)F(8661)',
    ],
    reactants = [
        ('CH2CO(28)', 'FC(F)=CC(F)F(2834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2319',
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

