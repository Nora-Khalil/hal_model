species(
    label = 'O=[C]CC1[CH]C1F(8609)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {12,S}
7  C u1 p0 c0 {2,D} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (129.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,1855,455,950,327.854,327.857,327.859,327.859,933.896,1930.08,4000],'cm^-1')),
        HinderedRotor(inertia=(0.394471,'amu*angstrom^2'), symmetry=1, barrier=(30.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0113822,'amu*angstrom^2'), symmetry=1, barrier=(30.0887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51485,0.0498018,-3.75289e-05,1.44163e-08,-2.24842e-12,15608.9,25.2761], Tmin=(100,'K'), Tmax=(1499.06,'K')), NASAPolynomial(coeffs=[12.3287,0.0209465,-8.65547e-06,1.57549e-09,-1.06937e-13,12366.8,-31.2724], Tmin=(1499.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(CCCJ=O)"""),
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
    label = 'FC1C=C1(6880)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {2,S} {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (53.2336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,725.22,725.471,725.523,725.698],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2885.09,'J/mol'), sigma=(4.85707,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=450.64 K, Pc=57.13 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.11981,-0.0104554,0.000119533,-1.9991e-07,1.08069e-10,6403.99,8.0387], Tmin=(10,'K'), Tmax=(590.766,'K')), NASAPolynomial(coeffs=[1.59625,0.0247584,-1.59029e-05,4.86605e-09,-5.67736e-13,6385.84,16.208], Tmin=(590.766,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(53.2336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1CDC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(=O)C1[CH]C1F(8607)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {2,D} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (82.8005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58893,0.0447985,-2.02566e-05,-3.37562e-09,3.56589e-12,10052.3,25.3298], Tmin=(100,'K'), Tmax=(1130.43,'K')), NASAPolynomial(coeffs=[12.1615,0.0220799,-9.60659e-06,1.84131e-09,-1.3065e-13,6723.22,-31.1253], Tmin=(1130.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.8005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(Cs(F)-Cs-Cs) + radical(CCJCC=O) + radical(C2JC=O)"""),
)

species(
    label = '[C]=O(128)',
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
    label = '[CH2]C1[CH]C1F(8659)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4  C u1 p0 c0 {2,S} {3,S} {8,S}
5  C u1 p0 c0 {2,S} {9,S} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (264.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,3000,3100,440,815,1455,1000,180,180,663.738,1463.19,1463.28,1463.42],'cm^-1')),
        HinderedRotor(inertia=(0.0461459,'amu*angstrom^2'), symmetry=1, barrier=(70.1556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41676,0.027992,1.78753e-06,-2.02849e-08,9.33456e-12,31822.7,18.31], Tmin=(100,'K'), Tmax=(994.052,'K')), NASAPolynomial(coeffs=[8.95231,0.0175157,-6.27962e-06,1.13771e-09,-8.01616e-14,29741.6,-17.1133], Tmin=(994.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Isobutyl)"""),
)

species(
    label = 'O=C1CC2C(F)C12(8706)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
6  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {2,D} {4,S} {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-181.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96582,0.0375768,-8.28867e-06,-9.66214e-09,4.40371e-12,-21726.2,15.7401], Tmin=(100,'K'), Tmax=(1205.06,'K')), NASAPolynomial(coeffs=[10.038,0.0254597,-1.14754e-05,2.20793e-09,-1.55649e-13,-24737.4,-29.1308], Tmin=(1205.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'O=CCC1=CC1F(8707)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u0 p0 c0 {2,D} {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-86.3416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88664,0.0345688,1.16541e-05,-3.54989e-08,1.40338e-11,-10298.3,22.427], Tmin=(100,'K'), Tmax=(1082.89,'K')), NASAPolynomial(coeffs=[12.5118,0.0231675,-1.11247e-05,2.27063e-09,-1.67888e-13,-14232.2,-37.2177], Tmin=(1082.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.3416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cd(C)-Cd-Cs(F))"""),
)

species(
    label = 'O=C=CC1CC1F(8708)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {7,D} {12,S}
7  C u0 p0 c0 {2,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-140.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70463,0.0424412,-1.56498e-05,-6.99285e-09,4.67344e-12,-16861.2,21.7252], Tmin=(100,'K'), Tmax=(1106.39,'K')), NASAPolynomial(coeffs=[11.338,0.0228254,-9.67963e-06,1.83715e-09,-1.29895e-13,-19923.9,-29.9316], Tmin=(1106.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(Cs(F)-Cs-Cs)"""),
)

species(
    label = 'O=CCC1C=C1F(8613)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u0 p0 c0 {2,D} {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-36.9921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18817,0.0453971,-3.13219e-05,1.16288e-08,-1.9658e-12,-4388.75,20.3648], Tmin=(100,'K'), Tmax=(1230.05,'K')), NASAPolynomial(coeffs=[6.20267,0.0323424,-1.54022e-05,3.00058e-09,-2.12169e-13,-5376.36,0.165909], Tmin=(1230.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.9921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cd-Cd(F)-Cs(C))"""),
)

species(
    label = 'O=[C]CC=C[CH]F(7632)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u1 p0 c0 {1,S} {5,S} {12,S}
7  C u1 p0 c0 {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-13.4448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,234,589,736,816,1240,3237,1855,455,950,1466.91],'cm^-1')),
        HinderedRotor(inertia=(0.828331,'amu*angstrom^2'), symmetry=1, barrier=(19.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183933,'amu*angstrom^2'), symmetry=1, barrier=(19.0239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828911,'amu*angstrom^2'), symmetry=1, barrier=(19.0583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37402,0.0503284,-3.47863e-05,1.07436e-08,-1.28402e-12,-1516.1,24.8771], Tmin=(100,'K'), Tmax=(1969.82,'K')), NASAPolynomial(coeffs=[19.3267,0.0138736,-7.02672e-06,1.34874e-09,-9.16947e-14,-8588.89,-73.9051], Tmin=(1969.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.4448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]=O(136)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = 'O=[C]CC1=CC1F(8709)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u1 p0 c0 {2,D} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (73.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,1000,1855,455,950,180,295.948,737.169,737.18,737.184,737.208],'cm^-1')),
        HinderedRotor(inertia=(0.0555469,'amu*angstrom^2'), symmetry=1, barrier=(21.42,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555442,'amu*angstrom^2'), symmetry=1, barrier=(21.42,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0829,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81873,0.0375019,-3.40929e-06,-2.10027e-08,9.6124e-12,8941.74,23.4878], Tmin=(100,'K'), Tmax=(1083.83,'K')), NASAPolynomial(coeffs=[13.1074,0.0186657,-8.93092e-06,1.82507e-09,-1.35211e-13,5154.07,-38.0671], Tmin=(1083.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cd(C)-Cd-Cs(F)) + radical(CCCJ=O)"""),
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
    label = 'O=[C]CC1C=C1(8710)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,D} {10,S}
5  C u0 p0 c0 {2,S} {4,D} {11,S}
6  C u1 p0 c0 {1,D} {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (294.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,180,180,180,191.53,953.554,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.193807,'amu*angstrom^2'), symmetry=1, barrier=(4.456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193807,'amu*angstrom^2'), symmetry=1, barrier=(4.456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91451,0.0482027,-4.89674e-05,2.91456e-08,-7.34594e-12,35488.8,18.237], Tmin=(100,'K'), Tmax=(941.535,'K')), NASAPolynomial(coeffs=[7.57689,0.0241461,-1.06406e-05,2.00702e-09,-1.39815e-13,34422.5,-8.73943], Tmin=(941.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=C1F(8711)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {3,S} {5,D}
7  C u1 p0 c0 {2,D} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (122.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,1855,455,950,242.196,242.268,242.507,1608.26,1608.27,1608.28],'cm^-1')),
        HinderedRotor(inertia=(0.158673,'amu*angstrom^2'), symmetry=1, barrier=(6.6093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00286661,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0829,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86255,0.0511447,-5.53223e-05,3.67515e-08,-1.06384e-11,14863,22.3665], Tmin=(100,'K'), Tmax=(815.639,'K')), NASAPolynomial(coeffs=[6.48824,0.0284606,-1.36068e-05,2.6565e-09,-1.88381e-13,14108.4,0.992688], Tmin=(815.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsH) + ring(Cd-Cd(F)-Cs(C)) + radical(CCCJ=O)"""),
)

species(
    label = 'FC1[CH][CH]1(8666)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {4,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (308.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,3150,900,1100,1800.48,1800.48,1800.48,1800.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60776,0.0050847,2.30871e-05,-2.54932e-08,7.79552e-12,37164.6,14.3884], Tmin=(100,'K'), Tmax=(1144.07,'K')), NASAPolynomial(coeffs=[3.87852,0.0162024,-7.30699e-06,1.43499e-09,-1.02849e-13,36313.1,9.595], Tmin=(1144.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'O=C=CC1[CH]C1F(8712)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {2,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (84.9779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180,1181.34,1182.77,1183.15,1183.33],'cm^-1')),
        HinderedRotor(inertia=(0.00737584,'amu*angstrom^2'), symmetry=1, barrier=(7.33565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0829,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07335,0.0405819,-2.50393e-05,7.25516e-09,-8.40334e-13,10290.8,22.7503], Tmin=(100,'K'), Tmax=(1937.11,'K')), NASAPolynomial(coeffs=[12.4581,0.019138,-8.43403e-06,1.54029e-09,-1.02774e-13,6267.63,-34.2161], Tmin=(1937.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.9779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane)"""),
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
    label = 'O=[C]C[C]1C=C1(8713)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {3,S} {4,D} {10,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (426.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,1855,455,950,215.805,215.844,215.845,215.856,1357.54,1357.54,1357.55,1357.55],'cm^-1')),
        HinderedRotor(inertia=(0.00361805,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155418,'amu*angstrom^2'), symmetry=1, barrier=(5.13602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0845,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10758,0.0438249,-4.71557e-05,3.12399e-08,-8.81774e-12,51355.5,15.7956], Tmin=(100,'K'), Tmax=(846.181,'K')), NASAPolynomial(coeffs=[6.56403,0.0227584,-9.81116e-06,1.81742e-09,-1.24866e-13,50601.3,-4.95979], Tmin=(846.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(Allyl_T) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[C]1CC1F(7750)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {2,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (75.9995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92375,0.0509459,-6.31665e-05,5.77208e-08,-2.24881e-11,9210.27,24.2126], Tmin=(100,'K'), Tmax=(784.173,'K')), NASAPolynomial(coeffs=[2.9022,0.0364709,-1.73367e-05,3.33565e-09,-2.32773e-13,9348.41,21.5893], Tmin=(784.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.9995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C[C]1F(8714)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {3,S} {4,S}
7  C u1 p0 c0 {2,D} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (123.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53319,0.0481406,-3.44718e-05,1.2441e-08,-1.81502e-12,14969.6,25.2142], Tmin=(100,'K'), Tmax=(1598.94,'K')), NASAPolynomial(coeffs=[12.7585,0.0200587,-8.12753e-06,1.45695e-09,-9.76262e-14,11379.8,-34.2101], Tmin=(1598.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CC1[CH]C1F(8715)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {6,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (104.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77235,0.0346628,2.01182e-05,-5.3429e-08,2.33356e-11,12631.6,24.5211], Tmin=(100,'K'), Tmax=(974.275,'K')), NASAPolynomial(coeffs=[14.601,0.0168222,-6.03747e-06,1.16115e-09,-8.76125e-14,8478.85,-45.5188], Tmin=(974.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cs(F)-Cs-Cs) + radical(C=COJ) + radical(cyclopropane)"""),
)

species(
    label = 'O=[C][CH]C1CC1F(8716)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {7,S} {12,S}
7  C u1 p0 c0 {2,D} {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (61.5072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74297,0.0417049,-1.35229e-05,-1.03536e-08,6.37711e-12,7485.71,23.1137], Tmin=(100,'K'), Tmax=(1039.89,'K')), NASAPolynomial(coeffs=[11.0515,0.0221495,-8.7552e-06,1.61723e-09,-1.13486e-13,4671.11,-26.3833], Tmin=(1039.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.5072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cs) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC[C]1[CH]C1F(8717)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {6,S}
6  C u1 p0 c0 {3,S} {5,S} {11,S}
7  C u0 p0 c0 {2,D} {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (154.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,2950,1000,2782.5,750,1395,475,1775,1000,180,180,180,1465.91,1467.22,1470.57],'cm^-1')),
        HinderedRotor(inertia=(0.00434181,'amu*angstrom^2'), symmetry=1, barrier=(6.59924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286833,'amu*angstrom^2'), symmetry=1, barrier=(6.59485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85322,0.0523037,-6.36555e-05,5.71794e-08,-2.24004e-11,18650.4,25.064], Tmin=(100,'K'), Tmax=(756.594,'K')), NASAPolynomial(coeffs=[3.33227,0.0368433,-1.78558e-05,3.47549e-09,-2.44584e-13,18645.2,19.7863], Tmin=(756.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'O=CCC1[CH][C]1F(8718)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {6,S}
6  C u1 p0 c0 {3,S} {5,S} {11,S}
7  C u0 p0 c0 {2,D} {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (198.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,212,367,445,1450,2782.5,750,1395,475,1775,1000,385.22,387.289,597.222,1699.32,1702.06,1703.67,1705.78],'cm^-1')),
        HinderedRotor(inertia=(0.22687,'amu*angstrom^2'), symmetry=1, barrier=(24.0156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222248,'amu*angstrom^2'), symmetry=1, barrier=(24.0161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50464,0.0475067,-2.8899e-05,5.64863e-09,5.56121e-13,24000,25.5938], Tmin=(100,'K'), Tmax=(1187.2,'K')), NASAPolynomial(coeffs=[12.1222,0.0218058,-9.1528e-06,1.70663e-09,-1.18639e-13,20769.1,-30.4412], Tmin=(1187.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'O=C(F)CC1[CH][CH]1(8719)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {6,S} {11,S}
6  C u1 p0 c0 {3,S} {5,S} {12,S}
7  C u0 p0 c0 {1,S} {2,D} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (100.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26533,0.0391942,-2.05707e-05,4.68438e-09,-4.01591e-13,12149.2,25.0145], Tmin=(100,'K'), Tmax=(2344.54,'K')), NASAPolynomial(coeffs=[18.5799,0.0139864,-6.12349e-06,1.05412e-09,-6.54406e-14,3777.31,-69.1348], Tmin=(2344.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(COCsFO) + ring(Cyclopropane) + radical(cyclopropane) + radical(cyclopropane)"""),
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
    E0 = (4.13984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (249.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (578.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (12.4242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (67.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (67.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (29.1131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (87.6271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (57.5402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (112.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (166.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (290.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (210.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (127.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (171.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (344.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (274.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (158.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (165.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (163.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (161.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (134.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (166.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (231.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['CH2CO(28)', 'FC1C=C1(6880)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['[CH2]C(=O)C1[CH]C1F(8607)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(128)', '[CH2]C1[CH]C1F(8659)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=C1CC2C(F)C12(8706)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=CCC1=CC1F(8707)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=C=CC1CC1F(8708)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=CCC1C=C1F(8613)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]CC=C[CH]F(7632)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CO(13)', '[CH2]C1[CH]C1F(8659)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(37.08,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=O(136)', 'FC1C=C1(6880)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(24.0539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'O=[C]CC1=CC1F(8709)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(5.79641,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'O=[C]CC1C=C1(8710)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.8449,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=[C]CC1C=C1F(8711)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.596503,'m^3/(mol*s)'), n=2.29663, Ea=(0.989468,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08546600576587018, var=0.6163705976545122, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2CO(28)', 'FC1[CH][CH]1(8666)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00100541,'m^3/(mol*s)'), n=2.87982, Ea=(4.33339,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'O=C=CC1[CH]C1F(8712)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.42983,'m^3/(mol*s)'), n=1.92039, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.021512250841106063, var=1.4639131859069563, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=O(136)', 'FC1[CH][CH]1(8666)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'O=[C]C[C]1C=C1(8713)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(253.861,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=[C]C[C]1CC1F(7750)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2] for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=[C]CC1C[C]1F(8714)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['[O]C=CC1[CH]C1F(8715)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=[C][CH]C1CC1F(8716)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.5088e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=CC[C]1[CH]C1F(8717)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.98524e+06,'s^-1'), n=1.75122, Ea=(105.056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=CCC1[CH][C]1F(8718)'],
    products = ['O=[C]CC1[CH]C1F(8609)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4559.43,'s^-1'), n=2.33789, Ea=(92.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;XH_out] for rate rule [R4H_SSS;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC1[CH]C1F(8609)'],
    products = ['O=C(F)CC1[CH][CH]1(8719)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(227.433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #3078',
    isomers = [
        'O=[C]CC1[CH]C1F(8609)',
    ],
    reactants = [
        ('CH2CO(28)', 'FC1C=C1(6880)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3078',
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

