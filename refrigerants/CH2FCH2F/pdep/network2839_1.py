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
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25068,'amu*angstrom^2'), symmetry=1, barrier=(28.7557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24767,'amu*angstrom^2'), symmetry=1, barrier=(28.6863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2493,'amu*angstrom^2'), symmetry=1, barrier=(28.724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02382,0.122557,-0.000209261,1.86533e-07,-6.50473e-11,-90851.9,36.6752], Tmin=(100,'K'), Tmax=(809.493,'K')), NASAPolynomial(coeffs=[12.3022,0.0388206,-2.09506e-05,4.15048e-09,-2.90819e-13,-92423.3,-21.179], Tmin=(809.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'F[CH]C(=C(F)F)C([CH]F)=C(F)F(7580)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {11,D}
8  C u0 p0 c0 {7,S} {10,S} {12,D}
9  C u1 p0 c0 {1,S} {7,S} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {7,D}
12 C u0 p0 c0 {5,S} {6,S} {8,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-836.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.904228,0.114211,-0.00015712,1.08447e-07,-2.95611e-11,-100494,32.7531], Tmin=(100,'K'), Tmax=(898.202,'K')), NASAPolynomial(coeffs=[18.1952,0.0291546,-1.50747e-05,3.01702e-09,-2.16298e-13,-103925,-57.3402], Tmin=(898.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-836.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCFF) + group(CdCFF) + radical(CsCdF1sH) + radical(CsCdF1sH)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3410.9,'J/mol'), sigma=(5.50955,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=532.77 K, Pc=46.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279359,0.0980193,-0.00010752,1.40221e-09,6.1285e-11,-72818.6,35.2415], Tmin=(100,'K'), Tmax=(490.853,'K')), NASAPolynomial(coeffs=[10.6187,0.0421147,-2.33214e-05,4.71771e-09,-3.36438e-13,-74175.2,-10.7615], Tmin=(490.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH][C]=C(F)C(=CF)C(F)(F)F(8688)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {4,S} {8,S} {12,D}
10 C u0 p0 c0 {5,S} {8,D} {13,S}
11 C u1 p0 c0 {6,S} {12,S} {14,S}
12 C u1 p0 c0 {9,D} {11,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-809.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.575003,0.1082,-0.000147169,1.0259e-07,-2.85849e-11,-97197.8,33.1037], Tmin=(100,'K'), Tmax=(874.543,'K')), NASAPolynomial(coeffs=[16.003,0.0323788,-1.71281e-05,3.46389e-09,-2.49753e-13,-100098,-44.6541], Tmin=(874.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-809.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'F[CH]C([C](F)F)=C(F)F(7786)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {9,D}
7  C u1 p0 c0 {1,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-683.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,161,297,490,584,780,1358,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.00041304,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000413002,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2987,0.0651462,-8.30533e-05,5.73355e-08,-1.62573e-11,-82156.9,26.972], Tmin=(100,'K'), Tmax=(850.044,'K')), NASAPolynomial(coeffs=[9.68447,0.0256869,-1.34247e-05,2.72914e-09,-1.97812e-13,-83582.6,-12.1223], Tmin=(850.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = '[CH]F(223)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'FC=[C]C(F)(F)[C]=C(F)F(8689)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {9,D} {11,S}
8  C u0 p0 c0 {4,S} {5,S} {10,D}
9  C u1 p0 c0 {6,S} {7,D}
10 C u1 p0 c0 {6,S} {8,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-450.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,615,860,1140,1343,3152,562,600,623,1070,1265,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16112,'amu*angstrom^2'), symmetry=1, barrier=(26.6963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14769,'amu*angstrom^2'), symmetry=1, barrier=(26.3877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00100308,0.101201,-0.000188768,1.76571e-07,-6.27138e-11,-54011,29.9203], Tmin=(100,'K'), Tmax=(844.435,'K')), NASAPolynomial(coeffs=[8.94986,0.0317825,-1.74608e-05,3.43553e-09,-2.37812e-13,-54559.3,-6.04488], Tmin=(844.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1C(F)C(=C(F)F)C1(F)F(8680)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {7,S} {8,S} {12,D}
11 C u0 p0 c0 {4,S} {9,D} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255262,0.0901062,-0.000101817,6.05483e-08,-1.48894e-11,-122669,26.5778], Tmin=(100,'K'), Tmax=(967.713,'K')), NASAPolynomial(coeffs=[13.0535,0.0372053,-1.9818e-05,4.0586e-09,-2.95798e-13,-125146,-34.7463], Tmin=(967.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(CdCFF) + ring(Cyclobutane)"""),
)

species(
    label = 'FC=[C]C(F)(F)[C]1C(F)C1(F)F(8690)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {12,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-642.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.527178,0.111495,-0.000186683,1.70302e-07,-6.12069e-11,-77153.5,36.2938], Tmin=(100,'K'), Tmax=(803.169,'K')), NASAPolynomial(coeffs=[9.33955,0.0430945,-2.29649e-05,4.54966e-09,-3.19598e-13,-78117.1,-5.2772], Tmin=(803.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-642.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCCFF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C(F)(F)C(=CF)C1(F)F(8691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
9  C u1 p0 c0 {7,S} {8,S} {11,S}
10 C u0 p0 c0 {7,S} {8,S} {12,D}
11 C u1 p0 c0 {5,S} {9,S} {13,S}
12 C u0 p0 c0 {6,S} {10,D} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-774.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67026,0.119728,-0.000219339,2.13306e-07,-7.96009e-11,-92972,32.177], Tmin=(100,'K'), Tmax=(819.311,'K')), NASAPolynomial(coeffs=[6.51846,0.0498414,-2.76941e-05,5.53862e-09,-3.89462e-13,-92982.2,6.05434], Tmin=(819.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-774.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFF) + group(CsCCFF) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C(=CF)C1(F)F(8692)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 C u1 p0 c0 {4,S} {5,S} {7,S}
12 C u0 p0 c0 {6,S} {9,D} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-658.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.747708,0.110775,-0.000146812,9.70296e-08,-2.54002e-11,-79020.2,34.8997], Tmin=(100,'K'), Tmax=(932.216,'K')), NASAPolynomial(coeffs=[18.2623,0.0292071,-1.55667e-05,3.17244e-09,-2.30247e-13,-82564.5,-55.4787], Tmin=(932.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(CsCCFF) + group(CsCsFHH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(C(F)=C=CF)=C(F)F(8693)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,D}
7  C u0 p0 c0 {1,S} {6,S} {11,D}
8  C u1 p0 c0 {2,S} {6,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 C u0 p0 c0 {5,S} {11,D} {13,S}
11 C u0 p0 c0 {7,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-588.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,113,247,382,1207,3490,540,610,2055,180,550.77],'cm^-1')),
        HinderedRotor(inertia=(1.25609,'amu*angstrom^2'), symmetry=1, barrier=(28.88,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93714,'amu*angstrom^2'), symmetry=1, barrier=(44.5386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.340951,0.104167,-0.00015686,1.21213e-07,-3.62483e-11,-70587.5,30.8456], Tmin=(100,'K'), Tmax=(676.138,'K')), NASAPolynomial(coeffs=[13.5417,0.0314369,-1.63589e-05,3.23776e-09,-2.28701e-13,-72679.6,-32.2855], Tmin=(676.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-588.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCddCF) + group(CdCFF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(CsCdF1sH)"""),
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
    label = 'C#CC(F)(F)C([CH]F)=C(F)F(8694)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {3,S} {7,S} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {7,D}
10 C u0 p0 c0 {6,S} {11,T}
11 C u0 p0 c0 {10,T} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-636.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.68268,'amu*angstrom^2'), symmetry=1, barrier=(38.6881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68208,'amu*angstrom^2'), symmetry=1, barrier=(38.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6818,'amu*angstrom^2'), symmetry=1, barrier=(38.6679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.616888,0.109625,-0.000176116,1.4693e-07,-4.82358e-11,-76442.1,31.5074], Tmin=(100,'K'), Tmax=(810.658,'K')), NASAPolynomial(coeffs=[13.8458,0.0305155,-1.54022e-05,2.97476e-09,-2.05826e-13,-78532.4,-33.6608], Tmin=(810.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-636.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'FC#CC(F)(F)C([CH]F)=C(F)F(8695)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {3,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u0 p0 c0 {7,S} {12,T}
12 C u0 p0 c0 {6,S} {11,T}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-736.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,2175,525,239,401,1367,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.59146,'amu*angstrom^2'), symmetry=1, barrier=(36.5909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59225,'amu*angstrom^2'), symmetry=1, barrier=(36.6089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59266,'amu*angstrom^2'), symmetry=1, barrier=(36.6185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.859289,0.118916,-0.000206756,1.84329e-07,-6.36523e-11,-88436.2,34.6963], Tmin=(100,'K'), Tmax=(828.569,'K')), NASAPolynomial(coeffs=[12.2451,0.0354074,-1.8926e-05,3.7113e-09,-2.5768e-13,-89912.9,-21.8667], Tmin=(828.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[C]=CC(F)(F)C([CH]F)=C(F)F(8696)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {7,S} {12,D} {13,S}
10 C u1 p0 c0 {3,S} {8,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 C u1 p0 c0 {6,S} {9,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-744.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.26999,'amu*angstrom^2'), symmetry=1, barrier=(29.1996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27318,'amu*angstrom^2'), symmetry=1, barrier=(29.273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27387,'amu*angstrom^2'), symmetry=1, barrier=(29.2888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.992038,0.120593,-0.000200667,1.7539e-07,-6.05154e-11,-89375.3,36.5795], Tmin=(100,'K'), Tmax=(792.563,'K')), NASAPolynomial(coeffs=[13.0284,0.0375356,-2.02e-05,4.01298e-09,-2.82295e-13,-91211.5,-25.3648], Tmin=(792.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-744.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(CF)=C(F)F(8697)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 C u1 p0 c0 {7,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-652.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02152,'amu*angstrom^2'), symmetry=1, barrier=(23.4868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01984,'amu*angstrom^2'), symmetry=1, barrier=(23.448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02175,'amu*angstrom^2'), symmetry=1, barrier=(23.492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18674,0.129948,-0.000237247,2.19662e-07,-7.77871e-11,-78299.5,36.5979], Tmin=(100,'K'), Tmax=(836.833,'K')), NASAPolynomial(coeffs=[10.9563,0.0410105,-2.24528e-05,4.42873e-09,-3.07674e-13,-79250.1,-13.359], Tmin=(836.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-652.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[CH][C]=C(F)C(=C(F)F)C(F)F(8698)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {3,S} {8,S} {12,D}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u1 p0 c0 {6,S} {12,S} {14,S}
12 C u1 p0 c0 {9,D} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-767.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.755774,0.115645,-0.000175214,1.35494e-07,-4.00159e-11,-92204.9,34.2912], Tmin=(100,'K'), Tmax=(643.423,'K')), NASAPolynomial(coeffs=[13.8238,0.0370856,-2.02268e-05,4.08213e-09,-2.91595e-13,-94331.1,-31.5605], Tmin=(643.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(F)=C(F)C([CH]F)=C(F)F(8699)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {10,S} {11,D}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {1,S} {8,D} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {7,D}
12 C u1 p0 c0 {6,S} {9,S} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-807.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.96833,0.118071,-0.000173225,1.30251e-07,-3.89827e-11,-96901,33.5487], Tmin=(100,'K'), Tmax=(817.534,'K')), NASAPolynomial(coeffs=[16.3328,0.0334247,-1.79231e-05,3.61317e-09,-2.58908e-13,-99730,-46.4344], Tmin=(817.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-807.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + radical(CsCdF1sH) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = 'F[C]=C(C(F)F)C(F)(F)[C]=CF(8700)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u1 p0 c0 {7,S} {10,D}
12 C u1 p0 c0 {6,S} {9,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-669.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,615,860,1140,1343,3152,1685,370,167,640,1190,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.976013,'amu*angstrom^2'), symmetry=1, barrier=(22.4404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975414,'amu*angstrom^2'), symmetry=1, barrier=(22.4267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968262,'amu*angstrom^2'), symmetry=1, barrier=(22.2623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25036,0.131366,-0.000240193,2.21862e-07,-7.83605e-11,-80354.4,36.9322], Tmin=(100,'K'), Tmax=(837.095,'K')), NASAPolynomial(coeffs=[11.3347,0.0404779,-2.2226e-05,4.38629e-09,-3.04716e-13,-81383.9,-15.1098], Tmin=(837.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'F[C]C(=CF)C(F)(F)C(F)=CF(8701)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u0 p0 c0 {4,S} {8,D} {13,S}
11 C u0 p0 c0 {5,S} {9,D} {14,S}
12 C u2 p0 c0 {6,S} {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-735.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,323,467,575,827,1418,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.833563,'amu*angstrom^2'), symmetry=1, barrier=(19.1652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43945,'amu*angstrom^2'), symmetry=1, barrier=(33.0958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833811,'amu*angstrom^2'), symmetry=1, barrier=(19.171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28504,0.130837,-0.000234559,2.13396e-07,-7.46898e-11,-88312.9,34.642], Tmin=(100,'K'), Tmax=(831.433,'K')), NASAPolynomial(coeffs=[12.2835,0.039119,-2.13887e-05,4.22302e-09,-2.93862e-13,-89655.3,-22.8178], Tmin=(831.433,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-735.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(=C(F)F)C(F)F(8702)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {6,S} {9,D}
11 C u1 p0 c0 {7,S} {12,D}
12 C u1 p0 c0 {11,D} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-682.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,182,240,577,636,1210,1413,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10206,'amu*angstrom^2'), symmetry=1, barrier=(25.3385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10008,'amu*angstrom^2'), symmetry=1, barrier=(25.293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09066,'amu*angstrom^2'), symmetry=1, barrier=(25.0765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18979,0.128759,-0.000230921,2.11326e-07,-7.45054e-11,-81869.8,36.781], Tmin=(100,'K'), Tmax=(827.237,'K')), NASAPolynomial(coeffs=[11.7356,0.0398793,-2.19223e-05,4.3445e-09,-3.03203e-13,-83105.6,-17.6695], Tmin=(827.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-682.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C([CH]F)=C(F)F(8703)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {3,S} {7,S} {12,D}
10 C u1 p0 c0 {4,S} {8,S} {13,S}
11 C u0 p0 c0 {5,S} {6,S} {8,D}
12 C u1 p0 c0 {9,D} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-736.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,246,474,533,1155,234,589,736,816,1240,3237,182,240,577,636,1210,1413,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37501,'amu*angstrom^2'), symmetry=1, barrier=(31.6142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37509,'amu*angstrom^2'), symmetry=1, barrier=(31.6159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37458,'amu*angstrom^2'), symmetry=1, barrier=(31.6042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10977,0.122725,-0.000202917,1.74695e-07,-5.93061e-11,-88399.8,36.4115], Tmin=(100,'K'), Tmax=(794.267,'K')), NASAPolynomial(coeffs=[14.0349,0.0359397,-1.91591e-05,3.7886e-09,-2.65745e-13,-90473.9,-31.0761], Tmin=(794.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCFF) + group(Cds-CdsHH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)(F)[C]=CF(7585)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {11,D} {13,S}
10 C u0 p0 c0 {6,S} {12,D} {14,S}
11 C u1 p0 c0 {7,S} {9,D}
12 C u1 p0 c0 {8,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-624.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([106,166,285,329,412,480,417,605,567,797,680,834,1102,1258,1148,1222,535,695,819,901,1081,1199,1292,1394,3088,3216,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.237073,'amu*angstrom^2'), symmetry=1, barrier=(5.45079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239935,'amu*angstrom^2'), symmetry=1, barrier=(5.51658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239283,'amu*angstrom^2'), symmetry=1, barrier=(5.50158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3401.78,'J/mol'), sigma=(5.80027,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=531.35 K, Pc=39.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979987,0.121352,-0.000206413,1.84139e-07,-6.45574e-11,-74967.1,37.4247], Tmin=(100,'K'), Tmax=(797.988,'K')), NASAPolynomial(coeffs=[12.2639,0.0389097,-2.12619e-05,4.24296e-09,-2.98902e-13,-76569.5,-20.2779], Tmin=(797.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-624.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = 'FC=[C]C(F)(F)[C]=CF(8080)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {9,D} {11,S}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  C u1 p0 c0 {5,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-249.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,535,695,819,901,1081,1199,1292,1394,3088,3216,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16788,'amu*angstrom^2'), symmetry=1, barrier=(26.8519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16585,'amu*angstrom^2'), symmetry=1, barrier=(26.8051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353558,0.0917964,-0.000166636,1.55121e-07,-5.51454e-11,-29843.7,26.9481], Tmin=(100,'K'), Tmax=(841.741,'K')), NASAPolynomial(coeffs=[8.21727,0.0310154,-1.66021e-05,3.24929e-09,-2.24798e-13,-30338.2,-4.70825], Tmin=(841.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1C(F)(F)C(=CF)C1(F)F(8704)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {7,S} {8,S} {12,D}
11 C u0 p0 c0 {5,S} {9,D} {13,S}
12 C u0 p0 c0 {6,S} {10,D} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-1047.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0712299,0.0939171,-0.00011012,6.72306e-08,-1.68002e-11,-125840,26.3188], Tmin=(100,'K'), Tmax=(958.567,'K')), NASAPolynomial(coeffs=[14.0139,0.0357348,-1.90721e-05,3.90729e-09,-2.84794e-13,-128513,-40.3561], Tmin=(958.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1047.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = 'F[CH][C]1C(F)C(=C(F)F)C1(F)F(8667)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
9  C u1 p0 c0 {7,S} {8,S} {11,S}
10 C u0 p0 c0 {7,S} {8,S} {12,D}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-747.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.514456,0.116263,-0.000212325,2.0846e-07,-7.85744e-11,-89800,31.8435], Tmin=(100,'K'), Tmax=(817.165,'K')), NASAPolynomial(coeffs=[5.55697,0.0513166,-2.84441e-05,5.69113e-09,-4.0058e-13,-89616.1,10.9749], Tmin=(817.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-747.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCCFF) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = 'F[C]C(=C(F)F)C(F)(F)C=CF(8705)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u0 p0 c0 {5,S} {9,D} {14,S}
12 C u2 p0 c0 {6,S} {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-769.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04034,0.123833,-0.000215511,1.94948e-07,-6.873e-11,-92420.2,34.7244], Tmin=(100,'K'), Tmax=(810.869,'K')), NASAPolynomial(coeffs=[11.8413,0.0396934,-2.1765e-05,4.33385e-09,-3.04289e-13,-93832.2,-20.5464], Tmin=(810.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-769.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(=CF)C(F)F(8706)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {14,S}
11 C u1 p0 c0 {7,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-669.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,194,682,905,1196,1383,3221,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00096,'amu*angstrom^2'), symmetry=1, barrier=(23.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0026,'amu*angstrom^2'), symmetry=1, barrier=(23.0518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00495,'amu*angstrom^2'), symmetry=1, barrier=(23.1058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20737,0.130165,-0.000236973,2.18653e-07,-7.72356e-11,-80395.3,36.7794], Tmin=(100,'K'), Tmax=(836.106,'K')), NASAPolynomial(coeffs=[11.2467,0.0405124,-2.21861e-05,4.37732e-09,-3.04185e-13,-81426.8,-14.7892], Tmin=(836.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=C(C(F)(F)F)C(F)(F)[C]=CF(8707)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {7,S} {10,D}
12 C u1 p0 c0 {9,D} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-729.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,219,296,586,564,718,793,1177,1228,350,440,435,1725,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0235,'amu*angstrom^2'), symmetry=1, barrier=(23.5322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02343,'amu*angstrom^2'), symmetry=1, barrier=(23.5306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02353,'amu*angstrom^2'), symmetry=1, barrier=(23.533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19981,0.12683,-0.000220292,1.95652e-07,-6.74726e-11,-87598.9,35.417], Tmin=(100,'K'), Tmax=(820.851,'K')), NASAPolynomial(coeffs=[13.3263,0.0365739,-1.97782e-05,3.90464e-09,-2.72367e-13,-89327.7,-27.7995], Tmin=(820.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C(F)F)C(F)(F)C(F)=CF(8708)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 C u2 p0 c0 {8,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-761.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26497,0.126973,-0.000205681,1.78184e-07,-6.09299e-11,-91367.7,35.4235], Tmin=(100,'K'), Tmax=(813.251,'K')), NASAPolynomial(coeffs=[12.7103,0.0428066,-2.19844e-05,4.27059e-09,-2.96436e-13,-93130.6,-25.9736], Tmin=(813.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-761.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(=CF)C(F)(F)F(8709)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {13,S}
11 C u1 p0 c0 {7,S} {12,D}
12 C u1 p0 c0 {11,D} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-729.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,219,296,586,564,718,793,1177,1228,350,440,435,1725,194,682,905,1196,1383,3221,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04538,'amu*angstrom^2'), symmetry=1, barrier=(24.0353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04803,'amu*angstrom^2'), symmetry=1, barrier=(24.0963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05651,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19984,0.126831,-0.000220293,1.95655e-07,-6.74737e-11,-87598.9,35.4171], Tmin=(100,'K'), Tmax=(820.844,'K')), NASAPolynomial(coeffs=[13.3263,0.0365738,-1.97781e-05,3.90463e-09,-2.72366e-13,-89327.7,-27.7998], Tmin=(820.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (-237.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-139.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-20.1622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (17.5593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (212.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (256.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-256.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-33.7886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-138.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-166.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (24.0887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-73.5796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-11.3274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-33.0665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (90.4563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (180.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-82.9817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-127.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-43.0167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-106.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (8.20665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-12.8027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-27.7046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-72.5127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-7.0231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (276.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-256.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-138.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (38.9417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-17.2636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-145.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-23.8921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-24.5172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-53.2662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['FC=C=C(F)F(5948)', 'FC=C=C(F)F(5948)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(27.3892,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 27.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH]C(=C(F)F)C([CH]F)=C(F)F(7580)'],
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
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH][C]=C(F)C(=CF)C(F)(F)F(8688)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(282.564,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CF(1054)', 'F[CH]C([C](F)F)=C(F)F(7786)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(223)', 'FC=[C]C(F)(F)[C]=C(F)F(8689)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['FC=C1C(F)C(=C(F)F)C1(F)F(8680)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_1H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['FC=[C]C(F)(F)[C]1C(F)C1(F)F(8690)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH][C]1C(F)(F)C(=CF)C1(F)F(8691)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH]C1([C](F)F)C(=CF)C1(F)F(8692)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(98.4045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 98.0 to 98.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'F[CH]C(C(F)=C=CF)=C(F)F(8693)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.5389,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC=C=C(F)F(5948)', 'F[CH][C]=C(F)F(7343)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CC(F)(F)C([CH]F)=C(F)F(8694)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(60.8867,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'FC#CC(F)(F)C([CH]F)=C(F)F(8695)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C]=C(F)F(7343)', 'F[CH][C]=C(F)F(7343)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -3.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHF(40)', 'FC=[C]C(F)(F)[C]=C(F)F(8689)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C]=CC(F)(F)C([CH]F)=C(F)F(8696)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C]=[C]C(F)(F)C(CF)=C(F)F(8697)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH][C]=C(F)C(=C(F)F)C(F)F(8698)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(221.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH]C(F)=C(F)C([CH]F)=C(F)F(8699)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(158.255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]=C(C(F)F)C(F)(F)[C]=CF(8700)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]C(=CF)C(F)(F)C(F)=CF(8701)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(231.15,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]C(F)(F)C(=C(F)F)C(F)F(8702)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(162.648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(F)C(F)(F)C([CH]F)=C(F)F(8703)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(172.145,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]C(F)(F)C(F)(F)[C]=CF(7585)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.14642e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]F(744)', 'FC=[C]C(F)(F)[C]=CF(8080)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['FC=C1C(F)(F)C(=CF)C1(F)F(8704)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_noH]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[CH][C]1C(F)C(=C(F)F)C1(F)F(8667)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CF2(43)', 'FC=[C]C(F)(F)[C]=CF(8080)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0195237,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['F[C]C(=C(F)F)C(F)(F)C=CF(8705)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.07654e+10,'s^-1'), n=1.20849, Ea=(247.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_Cd;XH_out] + [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]=[C]C(F)(F)C(=CF)C(F)F(8706)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C(F)(F)F)C(F)(F)[C]=CF(8707)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    products = ['[CH]C(=C(F)F)C(F)(F)C(F)=CF(8708)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(240.488,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]C(F)(F)C(=CF)C(F)(F)F(8709)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(184.741,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #2839',
    isomers = [
        'F[CH]C(=C(F)F)C(F)(F)[C]=CF(7582)',
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
    label = 'PDepNetwork #2839',
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

