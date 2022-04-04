species(
    label = 'F[CH]C(F)C1(F)C=C[C]1F(2530)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-373.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,346,659,817,1284,2750,3150,900,1100,334,575,1197,1424,3202,438.055,439.811,439.849,439.877,439.977,440.22,440.918,442.832],'cm^-1')),
        HinderedRotor(inertia=(0.0008665,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344607,'amu*angstrom^2'), symmetry=1, barrier=(47.558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394418,0.0846886,-9.34671e-05,5.44888e-08,-1.30483e-11,-44841,28.5446], Tmin=(100,'K'), Tmax=(997.307,'K')), NASAPolynomial(coeffs=[13.0584,0.0338954,-1.70709e-05,3.41999e-09,-2.46483e-13,-47367,-32.5175], Tmin=(997.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-373.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC1=CC=C1F(304)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {3,S} {6,D} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (40.9583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,247.893,1569.59,1569.59,2285.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89822,0.00612407,0.000102608,-2.16319e-07,1.34206e-10,4932.65,9.22518], Tmin=(10,'K'), Tmax=(550.566,'K')), NASAPolynomial(coeffs=[3.92734,0.0284965,-1.98739e-05,6.49561e-09,-8.00611e-13,4587.15,5.99356], Tmin=(550.566,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(40.9583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCCDC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=C(F)[CH][CH]1(923)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (89.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,206,336,431,607,515,611,528,696,1312,1446,180,692.409,1476.48,1721.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3235.35,'J/mol'), sigma=(5.2451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=505.35 K, Pc=50.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19216,0.0393646,-3.74e-05,1.78194e-08,-3.40083e-12,10888.1,15.3831], Tmin=(100,'K'), Tmax=(1250.71,'K')), NASAPolynomial(coeffs=[10.2921,0.0134598,-6.33194e-06,1.25922e-09,-9.07053e-14,8861.91,-25.5064], Tmin=(1250.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[CH]C(F)[C](F)C1(F)C=C1(3200)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {5,S} {8,D} {13,S}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-208.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,212,367,445,1450,2750,3150,900,1100,334,575,1197,1424,3202,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0779884,0.0954394,-0.000131312,9.5026e-08,-2.56154e-11,-24886,31.5012], Tmin=(100,'K'), Tmax=(630.325,'K')), NASAPolynomial(coeffs=[10.9165,0.0375071,-1.92646e-05,3.82223e-09,-2.71422e-13,-26467.8,-17.4956], Tmin=(630.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs(F)(C)-Cd) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)C1=CC(F)[C]1F(3201)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {3,S} {5,S} {7,S}
9  C u0 p0 c0 {5,S} {7,D} {13,S}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-343.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([96,252,209,325,527,655,638,804,1079,1135,1231,1325,1303,1393,3218,3328,346,659,817,1284,2950,1000,334,575,1197,1424,3202,503.939,505.443,505.831,506.025,506.84,507.451,509.242],'cm^-1')),
        HinderedRotor(inertia=(0.0115955,'amu*angstrom^2'), symmetry=1, barrier=(2.10596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011587,'amu*angstrom^2'), symmetry=1, barrier=(2.12297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.548825,0.0805615,-8.28058e-05,4.45375e-08,-9.87939e-12,-41139.5,28.8315], Tmin=(100,'K'), Tmax=(1067.75,'K')), NASAPolynomial(coeffs=[13.073,0.0336434,-1.68943e-05,3.38459e-09,-2.43977e-13,-43814.1,-32.4116], Tmin=(1067.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)C1C=C(F)[C]1F(2531)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {9,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-364.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,346,659,817,1284,271,519,563,612,1379,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.21,'J/mol'), sigma=(5.49523,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.17 K, Pc=46.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.427336,0.0768709,-7.02378e-05,3.17279e-08,-5.7617e-12,-43691.7,28.6364], Tmin=(100,'K'), Tmax=(1305.72,'K')), NASAPolynomial(coeffs=[16.6331,0.0272252,-1.32049e-05,2.6082e-09,-1.86267e-13,-47923.7,-53.87], Tmin=(1305.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd(F)-Cd-Cs-Cs) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = '[CH]F(137)',
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
    label = 'F[CH]C1(F)C=C[C]1F(3202)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u1 p0 c0 {2,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,D} {9,S}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-149.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,334,575,1197,1424,3202,470.04,470.041,470.042,470.042,470.043,470.043,470.048],'cm^-1')),
        HinderedRotor(inertia=(0.000763014,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47491,0.058568,-5.7649e-05,2.97324e-08,-6.34609e-12,-17891.2,22.2959], Tmin=(100,'K'), Tmax=(1104.16,'K')), NASAPolynomial(coeffs=[10.6414,0.0253609,-1.25376e-05,2.49536e-09,-1.79211e-13,-19915.4,-22.8357], Tmin=(1104.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(CsCsF1sH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1C(F)C2(F)C=CC12F(2532)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {4,S} {5,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {6,S} {9,D} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-515.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96607,0.0734498,-7.4211e-05,4.2387e-08,-1.0431e-11,-61847.5,21.6717], Tmin=(100,'K'), Tmax=(948.047,'K')), NASAPolynomial(coeffs=[9.16491,0.0388562,-1.94751e-05,3.8955e-09,-2.80444e-13,-63402.1,-17.4453], Tmin=(948.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC=C(F)C1(F)C=CC1F(3203)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-573.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709599,0.0793029,-8.96049e-05,5.72231e-08,-1.54558e-11,-68886,25.305], Tmin=(100,'K'), Tmax=(879.104,'K')), NASAPolynomial(coeffs=[9.77568,0.0380517,-1.92193e-05,3.84671e-09,-2.76727e-13,-70480,-17.2655], Tmin=(879.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs(F)-Cs-Cd-Cd)"""),
)

species(
    label = 'FCC(F)C1(F)C=C=C1F(2534)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {10,D} {14,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u0 p0 c0 {8,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-336.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0364537,0.0971608,-0.000154504,1.40671e-07,-5.1481e-11,-40312.8,28.2821], Tmin=(100,'K'), Tmax=(779.368,'K')), NASAPolynomial(coeffs=[7.87468,0.0427684,-2.25587e-05,4.4874e-09,-3.17266e-13,-41104.4,-4.81926], Tmin=(779.368,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'F[CH]C(F)C1(F)C2[CH]C21F(3204)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {12,S}
9  C u1 p0 c0 {6,S} {7,S} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-274.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.309167,0.0916071,-0.000120993,8.06822e-08,-1.70036e-11,-32850.7,27.7925], Tmin=(100,'K'), Tmax=(584.092,'K')), NASAPolynomial(coeffs=[9.69864,0.0400545,-2.13407e-05,4.3105e-09,-3.09485e-13,-34165,-14.3192], Tmin=(584.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-274.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_3_ane) + radical(bicyclo[1.1.0]butane-secondary) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C2[CH]C1(F)C(F)C2F(2982)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 C u1 p0 c0 {5,S} {6,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-287.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37482,0.0587623,-3.87903e-05,1.20045e-08,-1.49998e-12,-34504.9,25.7813], Tmin=(100,'K'), Tmax=(1777.49,'K')), NASAPolynomial(coeffs=[14.5505,0.0291124,-1.37692e-05,2.62002e-09,-1.80075e-13,-39188.8,-45.3623], Tmin=(1777.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFH) + polycyclic(s3_4_5_ane) + radical(CsCsCsF1s) + radical(bicyclo[2.1.1]hexane-C5) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C2C(F)C(F)C12F(3000)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
9  C u1 p0 c0 {5,S} {10,S} {14,S}
10 C u1 p0 c0 {4,S} {6,S} {9,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-277.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17319,0.0696913,-7.22511e-05,4.78953e-08,-1.45975e-11,-33256.8,25.6447], Tmin=(100,'K'), Tmax=(757.794,'K')), NASAPolynomial(coeffs=[5.87518,0.0448716,-2.31217e-05,4.67317e-09,-3.38168e-13,-33969.4,4.26451], Tmin=(757.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + polycyclic(s2_4_4_ane) + radical(bicyclo[2.2.0]hexane-secondary) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]=CC(F)=C(F)C(F)[CH]F(3205)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {9,S}
8  C u1 p0 c0 {4,S} {5,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u1 p0 c0 {9,D} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-266.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,280,518,736,852,873,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.718407,'amu*angstrom^2'), symmetry=1, barrier=(16.5176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.71661,'amu*angstrom^2'), symmetry=1, barrier=(16.4763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715972,'amu*angstrom^2'), symmetry=1, barrier=(16.4616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.413924,0.104327,-0.000149783,1.1234e-07,-3.36399e-11,-31877.2,31.8246], Tmin=(100,'K'), Tmax=(816.706,'K')), NASAPolynomial(coeffs=[14.4556,0.0315044,-1.60405e-05,3.17368e-09,-2.25184e-13,-34306.1,-36.9024], Tmin=(816.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Cds_P)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]C(F)C1=C(F)C=C1(3206)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {5,D} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {13,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-25.6549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2750,3150,900,1100,280,518,736,852,873,334,575,1197,1424,3202,180,180,180,575.933,578.97,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166955,'amu*angstrom^2'), symmetry=1, barrier=(3.83862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166955,'amu*angstrom^2'), symmetry=1, barrier=(3.83862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.45118,0.0805445,-9.41256e-05,5.66656e-08,-1.36442e-11,-2959.78,26.9885], Tmin=(100,'K'), Tmax=(1007.46,'K')), NASAPolynomial(coeffs=[14.3346,0.0254225,-1.20558e-05,2.35808e-09,-1.67952e-13,-5757.2,-40.0942], Tmin=(1007.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.6549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]F(141)',
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
    label = 'FC=CC1(F)C=C[C]1F(3207)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u1 p0 c0 {2,S} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {5,S} {6,D} {12,S}
9  C u0 p0 c0 {3,S} {7,D} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-252.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,478.388,478.394,478.432,478.433,478.434,478.434,478.443],'cm^-1')),
        HinderedRotor(inertia=(0.000736565,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10767,0.0635149,-5.34813e-05,2.28724e-08,-4.0154e-12,-30224.1,24.7092], Tmin=(100,'K'), Tmax=(1325.48,'K')), NASAPolynomial(coeffs=[12.9505,0.0277761,-1.30371e-05,2.53062e-09,-1.78736e-13,-33363.6,-35.7627], Tmin=(1325.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-252.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC=C(F)C1(F)C=C[C]1F(3208)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u1 p0 c0 {2,S} {5,S} {9,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u0 p0 c0 {6,S} {7,D} {12,S}
10 C u0 p0 c0 {4,S} {8,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-435.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,323,467,575,827,1418,194,682,905,1196,1383,3221,213.647,213.647,213.65,213.658,945.191,1555.93,1555.94],'cm^-1')),
        HinderedRotor(inertia=(1.22092,'amu*angstrom^2'), symmetry=1, barrier=(39.5465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735206,0.0775438,-8.64882e-05,5.19347e-08,-1.29232e-11,-52294.6,26.0808], Tmin=(100,'K'), Tmax=(957.852,'K')), NASAPolynomial(coeffs=[11.4121,0.0329583,-1.66691e-05,3.34176e-09,-2.4081e-13,-54340,-24.9695], Tmin=(957.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)C1(F)C=C=C1F(3209)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {10,D} {12,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {4,S} {6,S} {13,S}
10 C u0 p0 c0 {7,D} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-144.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,2950,1000,145,326,398,834,1303,334,575,1197,1424,3202,180,180,180,180,890.828,890.937],'cm^-1')),
        HinderedRotor(inertia=(0.123786,'amu*angstrom^2'), symmetry=1, barrier=(2.84608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123799,'amu*angstrom^2'), symmetry=1, barrier=(2.84639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.174141,0.105196,-0.00018699,1.76591e-07,-6.44032e-11,-17185.4,29.5631], Tmin=(100,'K'), Tmax=(821.683,'K')), NASAPolynomial(coeffs=[7.70202,0.0408815,-2.21685e-05,4.39532e-09,-3.07693e-13,-17602.9,-1.5526], Tmin=(821.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
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
    label = 'F[CH]C=C1[CH]C=C1F(3210)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {4,S} {5,D} {9,S}
7  C u0 p0 c0 {3,D} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {7,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (46.0915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,676.247,676.25,676.252,676.254,676.255,676.257,676.263,676.264,676.264],'cm^-1')),
        HinderedRotor(inertia=(0.000368627,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.093,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46828,0.0360139,3.59352e-05,-7.98633e-08,3.47504e-11,5652.51,20.7651], Tmin=(100,'K'), Tmax=(959.597,'K')), NASAPolynomial(coeffs=[18.7724,0.0126793,-3.86519e-06,7.79009e-10,-6.46787e-14,84.8802,-73.7095], Tmin=(959.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.0915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCCF) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Csj(Cd-CdH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)=C1[CH]C=C1F(3211)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,D}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u0 p0 c0 {2,S} {4,D} {9,S}
9  C u1 p0 c0 {3,S} {8,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-139.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,271,519,563,612,1379,234,589,736,816,1240,3237,713.293,713.663,714.975,715.038,716.169,716.402,717.288,717.407,717.961],'cm^-1')),
        HinderedRotor(inertia=(0.000327829,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.926382,0.0502763,2.44923e-06,-4.93575e-08,2.49253e-11,-16597,23.0841], Tmin=(100,'K'), Tmax=(959.843,'K')), NASAPolynomial(coeffs=[20.2715,0.0119051,-3.60698e-06,7.04392e-10,-5.72395e-14,-22256.7,-79.5895], Tmin=(959.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCCF) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)C1=C(F)C=[C]1(763)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {2,S} {5,D} {7,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u1 p0 c0 {3,S} {4,S} {12,S}
9  C u1 p0 c0 {5,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (193.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,280,518,736,852,873,2950,1000,334,575,1197,1424,3202,324.964,325.08,325.106,325.203,325.496,325.695,325.948,1710.52],'cm^-1')),
        HinderedRotor(inertia=(0.00431945,'amu*angstrom^2'), symmetry=1, barrier=(8.96175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00159374,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.329319,0.0861636,-0.000122386,9.0432e-08,-2.65868e-11,23433.3,27.9016], Tmin=(100,'K'), Tmax=(832.794,'K')), NASAPolynomial(coeffs=[13.034,0.0251407,-1.24723e-05,2.44228e-09,-1.72355e-13,21317.3,-31.0664], Tmin=(832.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(cyclobutadiene-C1)"""),
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
    label = 'FC[C](F)C1(F)C=C[C]1F(2971)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-375.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,346,659,817,1284,2750,3150,900,1100,371.003,371.049,371.072,371.097,371.105,371.154,1204.49,2051.96],'cm^-1')),
        HinderedRotor(inertia=(0.364541,'amu*angstrom^2'), symmetry=1, barrier=(35.6122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364388,'amu*angstrom^2'), symmetry=1, barrier=(35.6116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635545,0.0807837,-8.91578e-05,5.43765e-08,-1.39462e-11,-45065.4,29.1694], Tmin=(100,'K'), Tmax=(924.534,'K')), NASAPolynomial(coeffs=[10.6239,0.0375689,-1.90446e-05,3.81909e-09,-2.75154e-13,-46912.3,-18.2351], Tmin=(924.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCsCsF1s) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)C1(F)C=[C]C1F(691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {10,D} {13,S}
9  C u1 p0 c0 {4,S} {6,S} {14,S}
10 C u1 p0 c0 {7,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-258.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,2950,1000,334,575,1197,1424,3202,312.338,312.385,312.398,312.448,312.522,1271.04],'cm^-1')),
        HinderedRotor(inertia=(0.00761703,'amu*angstrom^2'), symmetry=1, barrier=(28.2461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407061,'amu*angstrom^2'), symmetry=1, barrier=(28.2466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.318923,0.0879857,-0.00010748,7.13154e-08,-1.94987e-11,-30948.8,28.4723], Tmin=(100,'K'), Tmax=(879.224,'K')), NASAPolynomial(coeffs=[11.8158,0.0356839,-1.82553e-05,3.66472e-09,-2.63816e-13,-32970.5,-25.514], Tmin=(879.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH][C](F)C1(F)C=CC1F(3212)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-321.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,2750,3150,900,1100,212,367,445,1450,334,575,1197,1424,3202,253.195,253.264,253.306,253.389,253.477,852.277,1986.48,1986.56],'cm^-1')),
        HinderedRotor(inertia=(0.64625,'amu*angstrom^2'), symmetry=1, barrier=(29.3972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104984,'amu*angstrom^2'), symmetry=1, barrier=(29.3974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675331,0.0853764,-9.31558e-05,2.22415e-08,2.69313e-11,-38541.3,28.7849], Tmin=(100,'K'), Tmax=(519.572,'K')), NASAPolynomial(coeffs=[9.16521,0.0401747,-2.08576e-05,4.14998e-09,-2.94708e-13,-39695.6,-9.23359], Tmin=(519.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)C1(F)[C]=CC1F(3213)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {10,D} {13,S}
9  C u1 p0 c0 {4,S} {7,S} {14,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-259.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,2950,1000,334,575,1197,1424,3202,192.254,194.814,202.503,209.415,211.9,994.772,1343.91],'cm^-1')),
        HinderedRotor(inertia=(0.00462983,'amu*angstrom^2'), symmetry=1, barrier=(0.131872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32573,'amu*angstrom^2'), symmetry=1, barrier=(36.1711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0716936,0.0943323,-0.00013379,1.06844e-07,-3.50902e-11,-31056.5,29.2981], Tmin=(100,'K'), Tmax=(739.169,'K')), NASAPolynomial(coeffs=[10.4624,0.0381018,-1.96781e-05,3.92237e-09,-2.79491e-13,-32592.5,-17.6904], Tmin=(739.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(cyclobutene-vinyl) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)C1(F)[C]=C[C]1F(3214)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
8  C u1 p0 c0 {4,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {14,S}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-313.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,346,659,817,1284,2950,1000,337.214,337.215,337.215,337.215,337.215,337.215,1005.57],'cm^-1')),
        HinderedRotor(inertia=(0.392295,'amu*angstrom^2'), symmetry=1, barrier=(31.6558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461051,'amu*angstrom^2'), symmetry=1, barrier=(37.2039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529182,0.0817731,-8.75871e-05,5.00057e-08,-1.18165e-11,-37602.1,28.008], Tmin=(100,'K'), Tmax=(1005.67,'K')), NASAPolynomial(coeffs=[12.3587,0.0347207,-1.74051e-05,3.48055e-09,-2.50596e-13,-39981.3,-29.1295], Tmin=(1005.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(cyclobutene-vinyl) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)C1(F)C=[C][C]1F(3215)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-312.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,346,659,817,1284,2950,1000,403.503,403.671,403.787,403.799,403.801,403.995,404.168],'cm^-1')),
        HinderedRotor(inertia=(0.309531,'amu*angstrom^2'), symmetry=1, barrier=(35.7996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309341,'amu*angstrom^2'), symmetry=1, barrier=(35.7988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.432574,0.0798364,-7.8366e-05,3.91136e-08,-7.92306e-12,-37480,28.3915], Tmin=(100,'K'), Tmax=(1173.54,'K')), NASAPolynomial(coeffs=[15.0662,0.0299578,-1.46119e-05,2.89606e-09,-2.07613e-13,-40914.7,-44.549], Tmin=(1173.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)[C]1C=CC1(F)F(3216)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-409.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.436657,0.0750309,-6.75271e-05,3.05296e-08,-5.54932e-12,-49094.4,26.9272], Tmin=(100,'K'), Tmax=(1308.67,'K')), NASAPolynomial(coeffs=[16.2279,0.0267637,-1.22025e-05,2.34561e-09,-1.65174e-13,-53227.5,-53.5045], Tmin=(1308.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-409.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(Allyl_T) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[C]1C=C[C]1C(F)C(F)F(3217)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {12,S}
7  C u1 p0 c0 {5,S} {8,S} {9,S}
8  C u1 p0 c0 {4,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {8,S} {9,D} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-452.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338745,0.0694371,-4.86079e-05,1.11712e-08,8.86261e-13,-54251.6,28.2788], Tmin=(100,'K'), Tmax=(1112.07,'K')), NASAPolynomial(coeffs=[17.7094,0.0239877,-1.02762e-05,1.96335e-09,-1.39652e-13,-59168.2,-62.1046], Tmin=(1112.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Allyl_T) + radical(CsCdCsF1s)"""),
)

species(
    label = 'F[CH][CH]C1(F)C=CC1(F)F(3218)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u1 p0 c0 {5,S} {10,S} {11,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-356.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,274,345,380,539,705,1166,1213,2750,3150,900,1100,3025,407.5,1350,352.5,334,575,1197,1424,3202,305.617,305.679,305.829,306.15,1516.22,1516.24,1516.27,1516.53],'cm^-1')),
        HinderedRotor(inertia=(0.00684113,'amu*angstrom^2'), symmetry=1, barrier=(11.1659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.435709,'amu*angstrom^2'), symmetry=1, barrier=(28.929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.488899,0.0848266,-0.000108788,8.04811e-08,-2.50047e-11,-42724.6,29.9968], Tmin=(100,'K'), Tmax=(773.398,'K')), NASAPolynomial(coeffs=[9.32787,0.0391068,-2.01053e-05,4.0289e-09,-2.8907e-13,-44091.7,-10.374], Tmin=(773.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C=CC1(F)[CH]C(F)F(2916)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
7  C u1 p0 c0 {5,S} {6,S} {12,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-402.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86595,0.0755531,-7.53752e-05,4.07764e-08,-9.36386e-12,-48340.3,29.2178], Tmin=(100,'K'), Tmax=(1016.21,'K')), NASAPolynomial(coeffs=[10.5558,0.0374119,-1.90761e-05,3.84238e-09,-2.77636e-13,-50309.7,-17.686], Tmin=(1016.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cs_S) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)[C](F)C1C=C1F(669)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-186.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,212,367,445,1450,323,467,575,827,1418,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314935,0.0910294,-0.00011782,7.40474e-08,-1.2683e-11,-22358.6,32.2717], Tmin=(100,'K'), Tmax=(590.184,'K')), NASAPolynomial(coeffs=[10.227,0.0380938,-1.94822e-05,3.85908e-09,-2.73784e-13,-23776.7,-12.4231], Tmin=(590.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd-Cd-Cs(C-F)) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'FC1=CC2C(F)C(F)C12F(2978)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-553.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.779855,0.0685549,-5.50752e-05,2.1699e-08,-3.45953e-12,-66396.7,22.2956], Tmin=(100,'K'), Tmax=(1458.7,'K')), NASAPolynomial(coeffs=[15.624,0.0278496,-1.32174e-05,2.56886e-09,-1.80893e-13,-70727.4,-54.9232], Tmin=(1458.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-553.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC=C(F)C1(F)CC=C1F(3219)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-587.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847589,0.0757729,-7.83651e-05,4.45683e-08,-1.07359e-11,-70524.3,24.493], Tmin=(100,'K'), Tmax=(975.281,'K')), NASAPolynomial(coeffs=[10.3164,0.0369372,-1.86347e-05,3.73843e-09,-2.69588e-13,-72371.3,-20.9516], Tmin=(975.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cd(F)-Cd-Cs-Cs)"""),
)

species(
    label = 'FC1C(F)C2(F)[CH][CH]C12F(3220)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {4,S} {6,S} {7,S} {12,S}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u1 p0 c0 {6,S} {9,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-238.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73305,0.0605122,-2.25556e-05,-7.06695e-08,8.23632e-11,-28652.1,25.1201], Tmin=(100,'K'), Tmax=(452.425,'K')), NASAPolynomial(coeffs=[5.15179,0.0455011,-2.3231e-05,4.65736e-09,-3.34742e-13,-29117.2,9.61725], Tmin=(452.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s2_4_4_ane) + radical(bicyclo[2.2.0]hexane-secondary) + radical(bicyclo[2.2.0]hexane-secondary) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=CC=C(F)C(F)[CH]F(3221)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {12,S}
8  C u1 p0 c0 {3,S} {5,S} {14,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-279.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,2995,3025,975,1000,1300,1375,400,500,1630,1680,334,575,1197,1424,3202,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.657777,'amu*angstrom^2'), symmetry=1, barrier=(15.1236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658219,'amu*angstrom^2'), symmetry=1, barrier=(15.1338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.659453,'amu*angstrom^2'), symmetry=1, barrier=(15.1621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.585562,0.107407,-0.000162187,1.27309e-07,-3.95137e-11,-33434.8,32.3243], Tmin=(100,'K'), Tmax=(792.019,'K')), NASAPolynomial(coeffs=[14.8949,0.0292323,-1.41454e-05,2.71002e-09,-1.87673e-13,-35887.2,-38.7519], Tmin=(792.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[CH]C(F)C1=C=C[CH]1(3222)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {8,D}
5  C u1 p0 c0 {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u1 p0 c0 {2,S} {3,S} {12,S}
8  C u0 p0 c0 {4,D} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (366.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2750,3150,900,1100,334,575,1197,1424,3202,180,1042.13,1042.13,1042.13,1042.13,1042.13,1042.13,1042.13,1042.13,1042.13,2279.93],'cm^-1')),
        HinderedRotor(inertia=(0.0414968,'amu*angstrom^2'), symmetry=1, barrier=(0.954093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0414968,'amu*angstrom^2'), symmetry=1, barrier=(0.954093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.093,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3759,0.0568431,-4.74145e-05,2.02978e-08,-3.56959e-12,44219.8,24.4939], Tmin=(100,'K'), Tmax=(1325.81,'K')), NASAPolynomial(coeffs=[11.9341,0.024989,-1.13754e-05,2.17607e-09,-1.52488e-13,41420.2,-29.4208], Tmin=(1325.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(C=CCJC=C) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)C1(F)[C]=C=C1(684)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {9,D} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {12,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (269.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,2950,1000,334,575,1197,1424,3202,180,180,180,180,180,180,965.742,966.011],'cm^-1')),
        HinderedRotor(inertia=(0.0931499,'amu*angstrom^2'), symmetry=1, barrier=(2.1417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.093313,'amu*angstrom^2'), symmetry=1, barrier=(2.14545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.158854,0.0998496,-0.000188519,1.84231e-07,-6.80799e-11,32553.6,27.1102], Tmin=(100,'K'), Tmax=(838.856,'K')), NASAPolynomial(coeffs=[5.6605,0.0401659,-2.19825e-05,4.34343e-09,-3.02198e-13,32807.5,8.54956], Tmin=(838.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(F)C1(F)C[C]=C1F(3223)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {4,S} {6,S} {14,S}
10 C u1 p0 c0 {7,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-263.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,2750,3150,900,1100,246,474,533,1155,334,575,1197,1424,3202,471.564,471.591,471.654,471.692,471.859,471.954,471.985,472.259],'cm^-1')),
        HinderedRotor(inertia=(0.000757109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000757129,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463479,0.0843458,-9.65662e-05,5.93775e-08,-1.50956e-11,-31602.3,27.687], Tmin=(100,'K'), Tmax=(939.646,'K')), NASAPolynomial(coeffs=[12.0276,0.0351184,-1.79822e-05,3.62322e-09,-2.61747e-13,-33775.5,-27.3832], Tmin=(939.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd(F)-Cd-Cs-Cs) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Cdj(Cs-CsHH)(Cd-CsF1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[CH][C](F)C1(F)CC=C1F(3224)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u1 p0 c0 {2,S} {5,S} {10,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-334.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2950,3150,900,1000,1100,323,467,575,827,1418,212,367,445,1450,334,575,1197,1424,3202,311.049,311.057,311.128,311.36,1475.36,1475.36,1475.43,1475.45],'cm^-1')),
        HinderedRotor(inertia=(0.0898111,'amu*angstrom^2'), symmetry=1, barrier=(6.20031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401019,'amu*angstrom^2'), symmetry=1, barrier=(27.7243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228385,0.0909191,-0.000128377,1.04781e-07,-3.55546e-11,-40154.4,29.9574], Tmin=(100,'K'), Tmax=(713.635,'K')), NASAPolynomial(coeffs=[9.37434,0.0396564,-2.06305e-05,4.12857e-09,-2.94927e-13,-41459.8,-11.0811], Tmin=(713.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd(F)-Cd-Cs-Cs) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[CH]C(F)[C]1C(F)=CC1F(765)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u1 p0 c0 {4,S} {6,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-390.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609651,0.0755535,-7.09664e-05,3.4422e-08,-6.83001e-12,-46802.9,26.1697], Tmin=(100,'K'), Tmax=(1191.29,'K')), NASAPolynomial(coeffs=[13.9593,0.0307292,-1.45262e-05,2.83695e-09,-2.01657e-13,-49983.5,-40.5711], Tmin=(1191.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-390.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Allyl_T) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]C1(F)C(F)=CC1F(3225)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  C u1 p0 c0 {5,S} {10,S} {12,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-320.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,323,467,575,827,1418,2950,1000,3025,407.5,1350,352.5,334,575,1197,1424,3202,180,180,180,1088.68,1088.96],'cm^-1')),
        HinderedRotor(inertia=(0.125731,'amu*angstrom^2'), symmetry=1, barrier=(2.89081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00343504,'amu*angstrom^2'), symmetry=1, barrier=(2.89159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484284,0.0849238,-0.000113049,8.88162e-08,-2.94047e-11,-38430,31.5386], Tmin=(100,'K'), Tmax=(727.28,'K')), NASAPolynomial(coeffs=[8.62866,0.0401309,-2.06654e-05,4.13387e-09,-2.95866e-13,-39614.7,-5.15993], Tmin=(727.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC(F)C(F)C1(F)[C]=C[CH]1(3226)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
8  C u1 p0 c0 {5,S} {9,S} {13,S}
9  C u0 p0 c0 {8,S} {10,D} {14,S}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-328.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,250,417,511,1155,1315,1456,3119,235,523,627,1123,1142,1372,1406,3097,2750,3150,900,1100,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621096,0.0767922,-7.39103e-05,3.67997e-08,-7.52323e-12,-39359.9,26.273], Tmin=(100,'K'), Tmax=(1154.9,'K')), NASAPolynomial(coeffs=[13.6163,0.0317835,-1.54527e-05,3.05523e-09,-2.18638e-13,-42361.5,-38.293], Tmin=(1154.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
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
    E0 = (-81.3606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (137.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (138.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (212.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (253.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-177.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-97.1072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-147.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (59.2451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-6.37299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-89.5077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (59.5424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (249.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (157.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (69.4008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-31.2644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-33.2754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (255.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (206.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (243.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (55.4147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (244.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (177.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-49.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (112.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-1.53687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (61.0961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-81.5385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (11.9707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-21.9169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (2.95299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (50.7299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-12.5177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (160.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-178.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-97.1072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-51.0184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (46.5368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (545.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (385.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (124.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-10.4686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-12.6112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (69.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (79.9673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['CHFCHF[Z](59)', 'FC1=CC=C1F(304)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(104.715,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 104.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(F)[C](F)C1(F)C=C1(3200)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[CH]C(F)C1=CC(F)[C]1F(3201)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(293.549,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(F)C1C=C(F)[C]1F(2531)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(389.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(137)', 'F[CH]C1(F)C=C[C]1F(3202)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC1C(F)C2(F)C=CC12F(2532)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC=C(F)C1(F)C=CC1F(3203)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FCC(F)C1(F)C=C=C1F(2534)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(38.1528,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[CH]C(F)C1(F)C2[CH]C21F(3204)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.00057e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[C]1C2[CH]C1(F)C(F)C2F(2982)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.96041e+10,'s^-1'), n=0.12, Ea=(179.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[C]1[CH]C2C(F)C(F)C12F(3000)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.2379e+10,'s^-1'), n=0.224681, Ea=(96.5681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 94.2 to 96.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(F)=C(F)C(F)[CH]F(3205)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'F[CH]C(F)C1=C(F)C=C1(3206)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(14.8639,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[CH][CH]F(141)', 'FC1=CC=C1F(304)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.61741e-05,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'FC=CC1(F)C=C[C]1F(3207)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(60.8745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHFCHF[Z](59)', 'FC1=C(F)[CH][CH]1(923)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00201083,'m^3/(mol*s)'), n=2.87982, Ea=(1.06964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'FC=C(F)C1(F)C=C[C]1F(3208)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(2.86354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'F[CH]C(F)C1(F)C=C=C1F(3209)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1090.54,'m^3/(mol*s)'), n=1.09894, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-5R!H-R_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-5R!H-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH][CH]F(141)', 'FC1=C(F)[CH][CH]1(923)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', 'F[CH]C=C1[CH]C=C1F(3210)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(18.438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'F[CH]C(F)=C1[CH]C=C1F(3211)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(287.778,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'F[CH]C(F)C1=C(F)C=[C]1(763)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(143.739,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CHF(40)', 'F[CH]C1(F)C=C[C]1F(3202)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC[C](F)C1(F)C=C[C]1F(2971)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.30951e+09,'s^-1'), n=1.14834, Ea=(136.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[CH]C(F)C1(F)C=[C]C1F(691)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH][C](F)C1(F)C=CC1F(3212)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(17992.9,'s^-1'), n=2.5106, Ea=(132.012,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy4;Y_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(F)C1(F)[C]=CC1F(3213)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['FCC(F)C1(F)[C]=C[C]1F(3214)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['FCC(F)C1(F)C=[C][C]1F(3215)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.59726e+10,'s^-1'), n=0.976078, Ea=(136.859,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[CH]C(F)[C]1C=CC1(F)F(3216)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(164.159,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[C]1C=C[C]1C(F)C(F)F(3217)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(189.029,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F[CH][CH]C1(F)C=CC1(F)F(3218)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(219.161,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[C]1C=CC1(F)[CH]C(F)F(2916)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(173.558,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction34',
    reactants = ['F[CH]C(F)[C](F)C1C=C1F(669)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC1=CC2C(F)C(F)C12F(2978)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC=C(F)C1(F)CC=C1F(3219)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['FC1C(F)C2(F)[CH][CH]C12F(3220)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.2379e+10,'s^-1'), n=0.224681, Ea=(135.057,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 132.5 to 135.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['F[C]=CC=C(F)C(F)[CH]F(3221)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingle]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['F2(78)', 'F[CH]C(F)C1=C=C[CH]1(3222)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction40',
    reactants = ['HF(38)', 'F[CH]C(F)C1(F)[C]=C=C1(684)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(209.062,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction41',
    reactants = ['F[CH]C(F)C1(F)C[C]=C1F(3223)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['F[CH][C](F)C1(F)CC=C1F(3224)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3367.88,'s^-1'), n=2.74476, Ea=(136.665,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_23cy4;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    products = ['F[CH]C(F)[C]1C(F)=CC1F(765)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(173.465,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction44',
    reactants = ['F[CH][CH]C1(F)C(F)=CC1F(3225)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(202.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction45',
    reactants = ['FC(F)C(F)C1(F)[C]=C[CH]1(3226)'],
    products = ['F[CH]C(F)C1(F)C=C[C]1F(2530)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(220.422,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #394',
    isomers = [
        'F[CH]C(F)C1(F)C=C[C]1F(2530)',
    ],
    reactants = [
        ('CHFCHF[Z](59)', 'FC1=CC=C1F(304)'),
        ('CHFCHF[Z](59)', 'FC1=C(F)[CH][CH]1(923)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #394',
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

