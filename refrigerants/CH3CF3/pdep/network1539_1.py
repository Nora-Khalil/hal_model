species(
    label = 'O=C[C](F)OC(F)(F)[C]=CF(5260)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
8  C u1 p0 c0 {3,S} {5,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-715.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,280,501,1494,1531,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,180,180,1590.28],'cm^-1')),
        HinderedRotor(inertia=(0.308241,'amu*angstrom^2'), symmetry=1, barrier=(7.08708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308375,'amu*angstrom^2'), symmetry=1, barrier=(7.09014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26114,'amu*angstrom^2'), symmetry=1, barrier=(51.988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26158,'amu*angstrom^2'), symmetry=1, barrier=(51.9983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.388823,0.108847,-0.00018847,1.73094e-07,-6.18421e-11,-85900.5,33.5812], Tmin=(100,'K'), Tmax=(819.238,'K')), NASAPolynomial(coeffs=[9.31038,0.039487,-2.11876e-05,4.18195e-09,-2.92192e-13,-86751.3,-6.77162], Tmin=(819.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(CsCOF1sO2s) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
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
    label = 'O=C[C]F(4375)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u2 p0 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (22.7053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,163,1167],'cm^-1')),
        HinderedRotor(inertia=(0.0337629,'amu*angstrom^2'), symmetry=1, barrier=(13.6228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34132,0.0165481,-1.72264e-05,1.26789e-08,-4.5452e-12,2752.64,10.5202], Tmin=(100,'K'), Tmax=(638.16,'K')), NASAPolynomial(coeffs=[4.07864,0.0119262,-6.36171e-06,1.32811e-09,-9.81565e-14,2658.54,7.2943], Tmin=(638.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C(F)(F)[C]=CF(547)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6 C u0 p0 c0 {3,S} {7,D} {8,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-300.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.813738,'amu*angstrom^2'), symmetry=1, barrier=(18.7094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97056,0.0474072,-5.58215e-05,3.38039e-08,-8.23918e-12,-36115.3,20.5538], Tmin=(100,'K'), Tmax=(990.624,'K')), NASAPolynomial(coeffs=[9.79057,0.0158305,-8.0071e-06,1.62531e-09,-1.18216e-13,-37664.6,-17.0994], Tmin=(990.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = '[C]=CF-2(1206)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.00704726,-1.17425e-06,-1.98484e-09,8.12221e-13,48709.9,8.54799], Tmin=(100,'K'), Tmax=(1284.6,'K')), NASAPolynomial(coeffs=[5.40191,0.00467991,-2.11332e-06,4.24428e-10,-3.06823e-14,47991.2,-1.49789], Tmin=(1284.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O=C[C](F)O[C](F)F(4092)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {4,S} {7,S}
7 C u0 p0 c0 {5,D} {6,S} {9,S}
8 C u1 p0 c0 {2,S} {3,S} {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-626.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,493,600,700,1144,1293,180,180,2179.74],'cm^-1')),
        HinderedRotor(inertia=(2.38341,'amu*angstrom^2'), symmetry=1, barrier=(54.7993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177265,'amu*angstrom^2'), symmetry=1, barrier=(4.07567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.38265,'amu*angstrom^2'), symmetry=1, barrier=(54.7817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45321,0.0645946,-0.000112687,1.07892e-07,-4.0038e-11,-75304.1,24.5605], Tmin=(100,'K'), Tmax=(819.35,'K')), NASAPolynomial(coeffs=[5.28342,0.0286641,-1.53621e-05,3.03579e-09,-2.12434e-13,-75353.4,10.3748], Tmin=(819.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-626.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFHO) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=CC1(F)OC(F)(F)C1=CF(5264)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1017.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.408671,0.0785493,-7.41011e-05,3.36164e-08,-6.09762e-12,-122225,26.2022], Tmin=(100,'K'), Tmax=(1307.65,'K')), NASAPolynomial(coeffs=[17.6482,0.0258147,-1.36092e-05,2.77619e-09,-2.01473e-13,-126733,-61.5924], Tmin=(1307.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1017.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C=C(F)OC(F)(F)C=CF(5274)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {8,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-950.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.300281,0.0996096,-0.000131256,8.60673e-08,-2.22973e-11,-114152,30.3282], Tmin=(100,'K'), Tmax=(942.993,'K')), NASAPolynomial(coeffs=[17.1774,0.0254732,-1.33294e-05,2.69798e-09,-1.9518e-13,-117449,-52.9658], Tmin=(942.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-950.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'FC=[C]C(F)(F)OC1(F)[CH]O1(5277)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {11,S}
9  C u1 p0 c0 {6,S} {7,S} {12,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-569.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.314655,0.103218,-0.000156656,1.25938e-07,-4.06624e-11,-68330.8,32.4676], Tmin=(100,'K'), Tmax=(757.278,'K')), NASAPolynomial(coeffs=[12.9739,0.0330247,-1.7616e-05,3.53158e-09,-2.51625e-13,-70343.3,-27.9474], Tmin=(757.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-569.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = 'FC=C1O[CH][C](F)OC1(F)F(5278)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u1 p0 c0 {6,S} {9,S} {12,S}
11 C u0 p0 c0 {4,S} {8,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-790.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.933422,0.0960609,-9.9437e-05,4.62667e-08,-8.22651e-12,-94925.9,27.5131], Tmin=(100,'K'), Tmax=(1379.38,'K')), NASAPolynomial(coeffs=[28.4342,0.0108991,-6.82815e-06,1.50796e-09,-1.14397e-13,-103028,-123.614], Tmin=(1379.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-790.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCFFO) + group(Cds-CdsCsOs) + group(CdCFH) + ring(Cyclohexanone) + radical(CsCsF1sO2s) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C1[C](F)OC(F)(F)C1=CF(5279)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {3,S} {5,S} {7,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-734.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0477635,0.0793605,-7.5619e-05,3.42945e-08,-6.0527e-12,-88131.2,31.7962], Tmin=(100,'K'), Tmax=(1377.19,'K')), NASAPolynomial(coeffs=[21.3249,0.0172851,-8.00898e-06,1.56644e-09,-1.11694e-13,-94018.2,-78.1552], Tmin=(1377.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-734.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
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
    label = 'O=C[C](F)OC(F)=C=CF(5280)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {1,S} {4,S} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {10,D}
8  C u0 p0 c0 {5,D} {6,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u0 p0 c0 {7,D} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-540.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,197,221,431,657,2782.5,750,1395,475,1775,1000,113,247,382,1207,3490,540,610,2055,620.96,620.964,620.968,620.97,620.971],'cm^-1')),
        HinderedRotor(inertia=(0.0118221,'amu*angstrom^2'), symmetry=1, barrier=(3.23486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18285,'amu*angstrom^2'), symmetry=1, barrier=(50.0331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182853,'amu*angstrom^2'), symmetry=1, barrier=(50.0332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550058,0.0811672,-9.95766e-05,6.35024e-08,-1.6383e-11,-64903.3,27.709], Tmin=(100,'K'), Tmax=(936.672,'K')), NASAPolynomial(coeffs=[12.9985,0.0280059,-1.44419e-05,2.90752e-09,-2.09775e-13,-67235.3,-31.5329], Tmin=(936.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-540.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(CdCddFO) + group(CdCddFH) + group(Cdd-CdsCds) + radical(CsCOF1sO2s)"""),
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
    label = 'O=C=C(F)OC(F)(F)[C]=CF(5281)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {11,D}
9  C u0 p0 c0 {4,S} {10,D} {12,S}
10 C u1 p0 c0 {7,S} {9,D}
11 C u0 p0 c0 {6,D} {8,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-687.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,197,221,431,657,615,860,1140,1343,3152,1685,370,2120,512.5,787.5,305.179,305.179,305.179,305.18,3116.93],'cm^-1')),
        HinderedRotor(inertia=(0.484339,'amu*angstrom^2'), symmetry=1, barrier=(32.0101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107896,'amu*angstrom^2'), symmetry=1, barrier=(7.13086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48434,'amu*angstrom^2'), symmetry=1, barrier=(32.0101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.549964,0.107738,-0.000170227,1.34699e-07,-4.19025e-11,-82479.6,32.2309], Tmin=(100,'K'), Tmax=(791.03,'K')), NASAPolynomial(coeffs=[15.778,0.0251696,-1.36491e-05,2.73345e-09,-1.93968e-13,-85062.6,-42.7138], Tmin=(791.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cd(Cdd-Od)FO) + group(CdCFH) + missing(Cdd-CdO2d) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = 'C#CC(F)(F)O[C](F)C=O(5282)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {11,S}
9  C u0 p0 c0 {6,S} {10,T}
10 C u0 p0 c0 {9,T} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-620.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([157,399,436,611,601,744,1127,280,501,1494,1531,2782.5,750,1395,475,1775,1000,2175,525,750,770,3400,2100,180,180,1199.75],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.32642,'amu*angstrom^2'), symmetry=1, barrier=(53.4889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.32036,'amu*angstrom^2'), symmetry=1, barrier=(53.3497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0271355,0.0949509,-0.000150742,1.27715e-07,-4.28206e-11,-74551.3,27.8405], Tmin=(100,'K'), Tmax=(810.56,'K')), NASAPolynomial(coeffs=[11.2728,0.0306348,-1.53975e-05,2.97319e-09,-2.05917e-13,-76084.6,-22.2642], Tmin=(810.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-620.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OC(F)(F)C#CF(5283)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u1 p0 c0 {3,S} {5,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {12,S}
10 C u0 p0 c0 {7,S} {11,T}
11 C u0 p0 c0 {4,S} {10,T}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-720.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([157,399,436,611,601,744,1127,280,501,1494,1531,2782.5,750,1395,475,1775,1000,2175,525,239,401,1367,190.692,190.695,190.696,1556.59],'cm^-1')),
        HinderedRotor(inertia=(0.205613,'amu*angstrom^2'), symmetry=1, barrier=(5.30583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02301,'amu*angstrom^2'), symmetry=1, barrier=(52.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02316,'amu*angstrom^2'), symmetry=1, barrier=(52.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0303624,'amu*angstrom^2'), symmetry=1, barrier=(52.2046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.213012,0.104212,-0.000181261,1.64929e-07,-5.81426e-11,-86545.6,31.0216], Tmin=(100,'K'), Tmax=(828.227,'K')), NASAPolynomial(coeffs=[9.66831,0.0355337,-1.89256e-05,3.71077e-09,-2.57858e-13,-87463.6,-10.4488], Tmin=(828.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + group(Cds-OdCsH) + group(Ct-CtCs) + group(CtCF) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C]F-2(1215)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]C(F)OC(F)(F)[C]=CF(5284)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {10,D} {13,S}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {6,D} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-696.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,355,410,600,1181,1341,1420,3056,615,860,1140,1343,3152,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.286213,'amu*angstrom^2'), symmetry=1, barrier=(6.58061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286794,'amu*angstrom^2'), symmetry=1, barrier=(6.59395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28653,'amu*angstrom^2'), symmetry=1, barrier=(6.5879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874989,'amu*angstrom^2'), symmetry=1, barrier=(20.1177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983075,0.123909,-0.00022673,2.05879e-07,-7.11389e-11,-83557.8,35.1749], Tmin=(100,'K'), Tmax=(851.202,'K')), NASAPolynomial(coeffs=[12.1432,0.0339234,-1.82815e-05,3.55882e-09,-2.44458e-13,-84767.1,-20.0143], Tmin=(851.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-696.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=C[C](F)OC(F)(F)C=[C]F(5285)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {11,D} {12,S}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {13,S}
11 C u1 p0 c0 {4,S} {8,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-728.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.350864,0.105957,-0.000175448,1.56415e-07,-5.52233e-11,-87484.4,32.9226], Tmin=(100,'K'), Tmax=(793.249,'K')), NASAPolynomial(coeffs=[10.4604,0.0376458,-2.01898e-05,4.01007e-09,-2.82271e-13,-88765.6,-13.9956], Tmin=(793.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-728.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(CsCOF1sO2s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'O=[C][C](F)OC(F)(F)C=CF(5286)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u1 p0 c0 {4,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 C u1 p0 c0 {6,D} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-820.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.609613,0.113759,-0.000199341,1.79894e-07,-6.25883e-11,-98491.2,33.0377], Tmin=(100,'K'), Tmax=(837.37,'K')), NASAPolynomial(coeffs=[10.9072,0.0358993,-1.89452e-05,3.69156e-09,-2.55138e-13,-99619,-15.6974], Tmin=(837.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-820.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)OC(F)(F)[C]=[C]F(5287)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-604.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,140,424,499,621,667,843,876,1082,2782.5,750,1395,475,1775,1000,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.215067,'amu*angstrom^2'), symmetry=1, barrier=(4.94482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439422,'amu*angstrom^2'), symmetry=1, barrier=(10.1032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215041,'amu*angstrom^2'), symmetry=1, barrier=(4.94422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89138,'amu*angstrom^2'), symmetry=1, barrier=(43.4865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.737946,0.116271,-0.000203415,1.83103e-07,-6.40199e-11,-72550.4,35.1082], Tmin=(100,'K'), Tmax=(817.522,'K')), NASAPolynomial(coeffs=[11.7706,0.0355378,-1.94474e-05,3.85826e-09,-2.6998e-13,-73942.9,-18.7263], Tmin=(817.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'O=CC(F)(F)OC(F)=[C][CH]F(5288)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {11,D}
10 C u1 p0 c0 {4,S} {11,S} {13,S}
11 C u1 p0 c0 {9,D} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-715.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000,293,496,537,1218,234,589,736,816,1240,3237,1685,370,312.158,312.158,312.159,312.159],'cm^-1')),
        HinderedRotor(inertia=(0.103494,'amu*angstrom^2'), symmetry=1, barrier=(7.15635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256364,'amu*angstrom^2'), symmetry=1, barrier=(17.7269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103494,'amu*angstrom^2'), symmetry=1, barrier=(7.15636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549372,'amu*angstrom^2'), symmetry=1, barrier=(37.9876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702975,0.111397,-0.00017099,1.32555e-07,-4.05576e-11,-85867.7,32.5531], Tmin=(100,'K'), Tmax=(802.84,'K')), NASAPolynomial(coeffs=[16.054,0.0279036,-1.4985e-05,3.00451e-09,-2.13917e-13,-88558.2,-44.6091], Tmin=(802.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFHH) + group(Cds-CdsCsH) + group(CdCFO) + group(Cds-OdCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-F1sO2s))"""),
)

species(
    label = 'O=C[C](F)OC(F)=C(F)[CH]F(5289)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {7,D}
9  C u1 p0 c0 {3,S} {5,S} {11,S}
10 C u1 p0 c0 {4,S} {7,S} {12,S}
11 C u0 p0 c0 {6,D} {9,S} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-767.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.224084,0.100435,-0.00013943,1.01104e-07,-2.94679e-11,-92128.4,31.0418], Tmin=(100,'K'), Tmax=(835.775,'K')), NASAPolynomial(coeffs=[14.0348,0.0321919,-1.69513e-05,3.40704e-09,-2.44359e-13,-94511.9,-35.1908], Tmin=(835.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = '[CH]=[C]C(F)(F)OC(F)(F)C=O(5290)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {10,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-652.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,140,424,499,621,667,843,876,1082,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.404725,'amu*angstrom^2'), symmetry=1, barrier=(9.30543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405217,'amu*angstrom^2'), symmetry=1, barrier=(9.31675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405099,'amu*angstrom^2'), symmetry=1, barrier=(9.31401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65828,'amu*angstrom^2'), symmetry=1, barrier=(38.127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.820666,0.118284,-0.000208074,1.8743e-07,-6.55697e-11,-78321.4,34.8351], Tmin=(100,'K'), Tmax=(814.927,'K')), NASAPolynomial(coeffs=[12.1135,0.0353312,-1.95543e-05,3.89533e-09,-2.7316e-13,-79783.1,-20.9516], Tmin=(814.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-652.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)O[C](F)C=O(5291)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {11,D}
9  C u1 p0 c0 {4,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 C u1 p0 c0 {8,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-720.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.468241,0.108085,-0.000177678,1.55689e-07,-5.39984e-11,-86508.9,32.7534], Tmin=(100,'K'), Tmax=(794.984,'K')), NASAPolynomial(coeffs=[11.4662,0.0360511,-1.91496e-05,3.78586e-09,-2.65735e-13,-88027.8,-19.7035], Tmin=(794.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCsCdF) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCOF1sO2s) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)[C]=CF(5257)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-663.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,486,617,768,1157,1926,615,860,1140,1343,3152,1685,370,238.573,238.573,238.574,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00296182,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15593,'amu*angstrom^2'), symmetry=1, barrier=(6.29798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683307,'amu*angstrom^2'), symmetry=1, barrier=(27.5986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0302886,0.0951656,-0.000136338,1.03265e-07,-3.1555e-11,-79650.6,37.3121], Tmin=(100,'K'), Tmax=(797.186,'K')), NASAPolynomial(coeffs=[12.6759,0.0317158,-1.69526e-05,3.42831e-09,-2.46437e-13,-81666.8,-20.8296], Tmin=(797.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-663.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = '[CH]=C(F)OC(F)(F)[C]=CF(5292)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {10,D}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {6,S} {8,D}
10 C u1 p0 c0 {7,D} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-394.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,293,496,537,1218,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,255.176,255.262,255.305,255.357],'cm^-1')),
        HinderedRotor(inertia=(0.220669,'amu*angstrom^2'), symmetry=1, barrier=(10.2137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220798,'amu*angstrom^2'), symmetry=1, barrier=(10.2126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.650543,'amu*angstrom^2'), symmetry=1, barrier=(30.1011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377959,0.103171,-0.000157501,1.19884e-07,-3.58413e-11,-47272.9,30.2302], Tmin=(100,'K'), Tmax=(822.586,'K')), NASAPolynomial(coeffs=[15.9505,0.0237721,-1.27192e-05,2.54827e-09,-1.8145e-13,-49959.3,-45.3564], Tmin=(822.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'FC=C1OC=C(F)OC1(F)F(5293)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {12,S}
11 C u0 p0 c0 {4,S} {8,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1001.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0848401,0.0725683,-4.84837e-05,8.61846e-09,1.44855e-12,-120305,24.3582], Tmin=(100,'K'), Tmax=(1207.3,'K')), NASAPolynomial(coeffs=[21.2988,0.0227211,-1.19452e-05,2.46437e-09,-1.80776e-13,-126917,-88.1524], Tmin=(1207.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1001.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsOs) + group(CdCFO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclohexane)"""),
)

species(
    label = '[O][CH]C1(F)OC(F)(F)C1=CF(5294)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-690.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.131797,0.100169,-0.000144973,1.13701e-07,-3.64435e-11,-82862.6,28.3647], Tmin=(100,'K'), Tmax=(757.453,'K')), NASAPolynomial(coeffs=[11.7643,0.0373477,-2.05661e-05,4.2047e-09,-3.03959e-13,-84664.8,-25.7225], Tmin=(757.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = 'O=C=[C]OC(F)(F)[C]=CF(5295)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u1 p0 c0 {6,S} {7,D}
9  C u1 p0 c0 {4,S} {10,D}
10 C u0 p0 c0 {5,D} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-259.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,615,860,1140,1343,3152,1670,1700,300,440,2120,512.5,787.5,223.395,223.397,223.409,223.42],'cm^-1')),
        HinderedRotor(inertia=(0.326456,'amu*angstrom^2'), symmetry=1, barrier=(11.559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32642,'amu*angstrom^2'), symmetry=1, barrier=(11.559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524556,'amu*angstrom^2'), symmetry=1, barrier=(18.5777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.300833,0.100504,-0.000162966,1.29334e-07,-3.98492e-11,-31016.4,32.2143], Tmin=(100,'K'), Tmax=(802.544,'K')), NASAPolynomial(coeffs=[16.1623,0.0184309,-9.53342e-06,1.85133e-09,-1.28265e-13,-33658.3,-43.5856], Tmin=(802.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=C(F)OC(F)(F)[C]=CF(5296)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,S} {13,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {11,D}
9  C u0 p0 c0 {4,S} {10,D} {12,S}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {6,S} {8,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-599.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,140,424,499,621,667,843,876,1082,293,496,537,1218,615,860,1140,1343,3152,1670,1700,300,440,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.981065,0.117139,-0.000184982,1.44501e-07,-4.41147e-11,-71872.4,35.8744], Tmin=(100,'K'), Tmax=(807.583,'K')), NASAPolynomial(coeffs=[17.6984,0.0246099,-1.31045e-05,2.60175e-09,-1.83591e-13,-74889.2,-50.2495], Tmin=(807.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(C=CJO)"""),
)

species(
    label = 'OC=C(F)OC(F)(F)[C]=[C]F(5297)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u1 p0 c0 {7,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-588.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,140,424,499,621,667,843,876,1082,326,540,652,719,1357,3010,987.5,1337.5,450,1655,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.844678,'amu*angstrom^2'), symmetry=1, barrier=(19.4208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845541,'amu*angstrom^2'), symmetry=1, barrier=(19.4407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845155,'amu*angstrom^2'), symmetry=1, barrier=(19.4318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.844922,'amu*angstrom^2'), symmetry=1, barrier=(19.4264,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13261,0.117286,-0.000171555,1.19874e-07,-3.23628e-11,-70616,33.7835], Tmin=(100,'K'), Tmax=(915.027,'K')), NASAPolynomial(coeffs=[21.4511,0.0185554,-9.69496e-06,1.93912e-09,-1.38729e-13,-74748.7,-73.1625], Tmin=(915.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-588.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=C=C(F)O[C](F)[CH]OF(5298)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {8,S}
7  C u1 p0 c0 {1,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u0 p0 c0 {3,S} {11,D} {13,S}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-257.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,395,473,707,1436,3025,407.5,1350,352.5,197,221,431,657,113,247,382,1207,3490,540,610,2055,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24998,0.123875,-0.000199254,1.57097e-07,-4.83284e-11,-30773.5,35.6627], Tmin=(100,'K'), Tmax=(801.673,'K')), NASAPolynomial(coeffs=[18.6344,0.0246643,-1.36308e-05,2.73937e-09,-1.94396e-13,-33961.8,-55.8732], Tmin=(801.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCddFO) + group(CdCddFH) + group(Cdd-CdsCds) + radical(CsCsF1sO2s) + radical(CCsJO)"""),
)

species(
    label = 'FC=[C]C(F)(F)O[C]=COF(5299)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {4,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {6,S} {11,D} {12,S}
9  C u0 p0 c0 {3,S} {10,D} {13,S}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {5,S} {8,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-282.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,140,424,499,621,667,843,876,1082,3010,987.5,1337.5,450,1655,615,860,1140,1343,3152,1670,1700,300,440,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01671,0.114112,-0.000166282,1.1653e-07,-3.15472e-11,-33835.1,36.4533], Tmin=(100,'K'), Tmax=(912.955,'K')), NASAPolynomial(coeffs=[20.8113,0.0184754,-9.14813e-06,1.78632e-09,-1.26023e-13,-37820.7,-66.8661], Tmin=(912.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC(F)(F)C(F)=CF(5300)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {12,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-729.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26575,0.109098,-0.000139893,8.32868e-08,-1.88227e-11,-87528,35.3112], Tmin=(100,'K'), Tmax=(1100.34,'K')), NASAPolynomial(coeffs=[26.2582,0.0090409,-3.49134e-06,6.44033e-10,-4.58011e-14,-93585,-100.107], Tmin=(1100.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C(F)(F)OC(F)=COF(5301)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u1 p0 c0 {7,S} {11,D}
11 C u1 p0 c0 {10,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-256.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,140,424,499,621,667,843,876,1082,326,540,652,719,1357,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.498982,'amu*angstrom^2'), symmetry=1, barrier=(11.4726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38233,'amu*angstrom^2'), symmetry=1, barrier=(31.7824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497473,'amu*angstrom^2'), symmetry=1, barrier=(11.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38477,'amu*angstrom^2'), symmetry=1, barrier=(31.8386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.904102,0.118318,-0.000188734,1.49113e-07,-4.51144e-11,-30624.9,34.2619], Tmin=(100,'K'), Tmax=(668.271,'K')), NASAPolynomial(coeffs=[15.7342,0.0302335,-1.68462e-05,3.40185e-09,-2.42164e-13,-33105.5,-41.2242], Tmin=(668.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-HH)) + radical(Cds_P)"""),
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
    E0 = (-285.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (151.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (208.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-277.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-246.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-101.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-191.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-213.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (8.51096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-149.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-45.2232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-49.2253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-78.8953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (14.9066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (164.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-129.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-81.2454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-133.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-141.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-58.1737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-119.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-54.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-134.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-96.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (278.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-277.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-210.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (172.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-6.18575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-3.01085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (215.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (215.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-93.7723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (216.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[C]F(4375)', '[O]C(F)(F)[C]=CF(547)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=CF-2(1206)', 'O=C[C](F)O[C](F)F(4092)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=CC1(F)OC(F)(F)C1=CF(5264)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=C=C(F)OC(F)(F)C=CF(5274)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['FC=[C]C(F)(F)OC1(F)[CH]O1(5277)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['FC=C1O[CH][C](F)OC1(F)F(5278)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.18632e+11,'s^-1'), n=0.442892, Ea=(58.383,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['[O]C1[C](F)OC(F)(F)C1=CF(5279)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.83275e+09,'s^-1'), n=0.690186, Ea=(94.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_cddouble] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(=O)F(4234)', 'F[CH][C]=C(F)F(2842)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.000388714,'m^3/(mol*s)'), n=2.45366, Ea=(40.1583,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.450791267159554, var=0.5685704571032071, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O=C[C](F)OC(F)=C=CF(5280)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(46.2032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]C(=O)F(398)', 'FC=C=C(F)F(1325)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C=C(F)OC(F)(F)[C]=CF(5281)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CC(F)(F)O[C](F)C=O(5282)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(68.8176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=C[C](F)OC(F)(F)C#CF(5283)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(398)', 'F[CH][C]=C(F)F(2842)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[C]F-2(1215)', '[O]C(F)(F)[C]=CF(547)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]C(F)OC(F)(F)[C]=CF(5284)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_noH] for rate rule [R2H_S;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=C[C](F)OC(F)(F)C=[C]F(5285)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=[C][C](F)OC(F)(F)C=CF(5286)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=CC(F)OC(F)(F)[C]=[C]F(5287)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=CC(F)(F)OC(F)=[C][CH]F(5288)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(227.069,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['O=C[C](F)OC(F)=C(F)[CH]F(5289)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.58534e+16,'s^-1'), n=-0.733083, Ea=(166.015,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0717849990446407, var=47.27832310158784, Tref=1000.0, N=3, data_mean=0.0, correlation='R2F_Ext-1R!H-R_N-4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_N-4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]C(F)(F)OC(F)(F)C=O(5290)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(167.639,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['[CH]=C(F)C(F)(F)O[C](F)C=O(5291)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(151.145,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(136.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(6)', '[CH]=C(F)OC(F)(F)[C]=CF(5292)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['FC=C1OC=C(F)OC1(F)F(5293)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['[O][CH]C1(F)OC(F)(F)C1=CF(5294)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(74.6464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'O=C=[C]OC(F)(F)[C]=CF(5295)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(282.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O[C]=C(F)OC(F)(F)[C]=CF(5296)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['OC=C(F)OC(F)(F)[C]=[C]F(5297)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.34748e+07,'s^-1'), n=1.47154, Ea=(155.58,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_single;XH_out] for rate rule [R7HJ_1;Cd_rad_out_single;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FC=C=C(F)O[C](F)[CH]OF(5298)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(42.7702,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FC=[C]C(F)(F)O[C]=COF(5299)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(68.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    products = ['[O]C=[C]OC(F)(F)C(F)=CF(5300)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(191.61,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]C(F)(F)OC(F)=COF(5301)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(42.6385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1539',
    isomers = [
        'O=C[C](F)OC(F)(F)[C]=CF(5260)',
    ],
    reactants = [
        ('O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1539',
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

