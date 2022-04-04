species(
    label = 'F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)',
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
        HinderedRotor(inertia=(1.24928,'amu*angstrom^2'), symmetry=1, barrier=(28.7234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24939,'amu*angstrom^2'), symmetry=1, barrier=(28.726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24902,'amu*angstrom^2'), symmetry=1, barrier=(28.7175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02395,0.122559,-0.000209268,1.86544e-07,-6.5053e-11,-90851.9,36.6757], Tmin=(100,'K'), Tmax=(809.46,'K')), NASAPolynomial(coeffs=[12.3025,0.0388201,-2.09504e-05,4.15042e-09,-2.90813e-13,-92423.4,-21.1803], Tmin=(809.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'FC=C=C(F)F(5101)',
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
    label = 'F[CH]C(=C(F)F)C([CH]F)=C(F)F(9827)',
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
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,173,295,515,663,653,819,693,939,1188,1292,3205,3269,141,223,164,316,539,615,578,694,1133,1287,1372,1454,257.048],'cm^-1')),
        HinderedRotor(inertia=(0.00261358,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937037,'amu*angstrom^2'), symmetry=1, barrier=(43.6866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934722,'amu*angstrom^2'), symmetry=1, barrier=(43.6917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.904228,0.114211,-0.00015712,1.08447e-07,-2.95612e-11,-100494,32.7531], Tmin=(100,'K'), Tmax=(898.201,'K')), NASAPolynomial(coeffs=[18.1952,0.0291546,-1.50747e-05,3.01702e-09,-2.16298e-13,-103925,-57.3402], Tmin=(898.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-836.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCFF) + group(CdCFF) + radical(CsCdF1sH) + radical(CsCdF1sH)"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)[C]=C(F)F(9831)',
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
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,164,312,561,654,898,1207,1299,3167,615,860,1140,1343,3152,562,600,623,1070,1265,1670,1700,300,440,180,180,1095.32],'cm^-1')),
        HinderedRotor(inertia=(0.244261,'amu*angstrom^2'), symmetry=1, barrier=(5.61603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244249,'amu*angstrom^2'), symmetry=1, barrier=(5.61577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244401,'amu*angstrom^2'), symmetry=1, barrier=(5.61926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3410.9,'J/mol'), sigma=(5.50955,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=532.77 K, Pc=46.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.278961,0.098026,-0.000107558,1.48751e-09,6.12189e-11,-72818.6,35.2428], Tmin=(100,'K'), Tmax=(490.875,'K')), NASAPolynomial(coeffs=[10.6188,0.0421145,-2.33213e-05,4.71767e-09,-3.36435e-13,-74175.2,-10.7619], Tmin=(490.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH][C]=C(F)C(=CF)C(F)(F)F(10778)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.574973,0.1082,-0.000147168,1.02589e-07,-2.85844e-11,-97197.8,33.1036], Tmin=(100,'K'), Tmax=(874.594,'K')), NASAPolynomial(coeffs=[16.0031,0.0323787,-1.7128e-05,3.46387e-09,-2.49751e-13,-100098,-44.6544], Tmin=(874.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-809.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = '[C]=CF-2(3543)',
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
    label = 'F[CH]C([C](F)F)=C(F)F(9960)',
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
        HinderedRotor(inertia=(0.000412985,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000413036,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29869,0.0651463,-8.30537e-05,5.7336e-08,-1.62574e-11,-82156.9,26.972], Tmin=(100,'K'), Tmax=(850.021,'K')), NASAPolynomial(coeffs=[9.68445,0.0256869,-1.34247e-05,2.72914e-09,-1.97813e-13,-83582.5,-12.1222], Tmin=(850.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
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
    label = 'FC=[C]C(F)(F)[C]=C(F)F(10779)',
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
        HinderedRotor(inertia=(1.15114,'amu*angstrom^2'), symmetry=1, barrier=(26.4669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15734,'amu*angstrom^2'), symmetry=1, barrier=(26.6096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00106365,0.101202,-0.000188771,1.76576e-07,-6.27161e-11,-54011,29.9205], Tmin=(100,'K'), Tmax=(844.425,'K')), NASAPolynomial(coeffs=[8.94999,0.0317823,-1.74607e-05,3.43549e-09,-2.3781e-13,-54559.4,-6.04557], Tmin=(844.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1C(F)C(=C(F)F)C1(F)F(10770)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255265,0.0901062,-0.000101817,6.05482e-08,-1.48894e-11,-122669,26.5778], Tmin=(100,'K'), Tmax=(967.717,'K')), NASAPolynomial(coeffs=[13.0535,0.0372052,-1.9818e-05,4.0586e-09,-2.95797e-13,-125146,-34.7464], Tmin=(967.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(CdCFF) + ring(Cyclobutane)"""),
)

species(
    label = 'FC=[C]C(F)(F)[C]1C(F)C1(F)F(10780)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.527283,0.111496,-0.000186688,1.70311e-07,-6.12115e-11,-77153.5,36.2942], Tmin=(100,'K'), Tmax=(803.146,'K')), NASAPolynomial(coeffs=[9.33972,0.0430942,-2.29647e-05,4.54961e-09,-3.19594e-13,-78117.2,-5.27817], Tmin=(803.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-642.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCCFF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C(F)(F)C(=CF)C1(F)F(10781)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.670294,0.119729,-0.00021934,2.13309e-07,-7.96023e-11,-92972,32.1771], Tmin=(100,'K'), Tmax=(819.307,'K')), NASAPolynomial(coeffs=[6.51852,0.0498413,-2.7694e-05,5.53861e-09,-3.89461e-13,-92982.2,6.05399], Tmin=(819.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-774.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFF) + group(CsCCFF) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C(=CF)C1(F)F(10782)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.747693,0.110775,-0.000146811,9.70291e-08,-2.54e-11,-79020.2,34.8997], Tmin=(100,'K'), Tmax=(932.24,'K')), NASAPolynomial(coeffs=[18.2623,0.029207,-1.55666e-05,3.17243e-09,-2.30246e-13,-82564.5,-55.4789], Tmin=(932.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(CsCCFF) + group(CsCsFHH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(C(F)=C=CF)=C(F)F(10783)',
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
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,113,247,382,1207,3490,540,610,2055,180,550.742],'cm^-1')),
        HinderedRotor(inertia=(1.25604,'amu*angstrom^2'), symmetry=1, barrier=(28.8788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93751,'amu*angstrom^2'), symmetry=1, barrier=(44.5472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.340684,0.104164,-0.000156842,1.21181e-07,-3.62297e-11,-70587.5,30.8447], Tmin=(100,'K'), Tmax=(676.051,'K')), NASAPolynomial(coeffs=[13.5414,0.0314374,-1.63592e-05,3.23783e-09,-2.28706e-13,-72679.5,-32.2843], Tmin=(676.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-588.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCddCF) + group(CdCFF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(CsCdF1sH)"""),
)

species(
    label = 'F[CH][C]=C(F)F(5078)',
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
    label = 'C#CC(F)(F)C([CH]F)=C(F)F(10784)',
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
        HinderedRotor(inertia=(1.681,'amu*angstrom^2'), symmetry=1, barrier=(38.6495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6831,'amu*angstrom^2'), symmetry=1, barrier=(38.6978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68168,'amu*angstrom^2'), symmetry=1, barrier=(38.6651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.616889,0.109625,-0.000176116,1.4693e-07,-4.82359e-11,-76442.1,31.5074], Tmin=(100,'K'), Tmax=(810.657,'K')), NASAPolynomial(coeffs=[13.8458,0.0305155,-1.54022e-05,2.97476e-09,-2.05826e-13,-78532.4,-33.6608], Tmin=(810.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-636.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'FC#CC(F)(F)C([CH]F)=C(F)F(10785)',
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
        HinderedRotor(inertia=(1.59229,'amu*angstrom^2'), symmetry=1, barrier=(36.6098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59129,'amu*angstrom^2'), symmetry=1, barrier=(36.5869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59217,'amu*angstrom^2'), symmetry=1, barrier=(36.6072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.859253,0.118915,-0.000206754,1.84326e-07,-6.36509e-11,-88436.2,34.6962], Tmin=(100,'K'), Tmax=(828.577,'K')), NASAPolynomial(coeffs=[12.245,0.0354076,-1.89261e-05,3.71132e-09,-2.57681e-13,-89912.8,-21.8664], Tmin=(828.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[C]=CC(F)(F)C([CH]F)=C(F)F(10786)',
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
        HinderedRotor(inertia=(1.27162,'amu*angstrom^2'), symmetry=1, barrier=(29.2371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2721,'amu*angstrom^2'), symmetry=1, barrier=(29.2482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27215,'amu*angstrom^2'), symmetry=1, barrier=(29.2492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.991845,0.12059,-0.000200656,1.75373e-07,-6.05067e-11,-89375.3,36.5788], Tmin=(100,'K'), Tmax=(792.629,'K')), NASAPolynomial(coeffs=[13.0281,0.0375361,-2.02004e-05,4.01306e-09,-2.82302e-13,-91211.4,-25.3631], Tmin=(792.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-744.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(CF)=C(F)F(10787)',
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
        HinderedRotor(inertia=(1.02112,'amu*angstrom^2'), symmetry=1, barrier=(23.4775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02133,'amu*angstrom^2'), symmetry=1, barrier=(23.4825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02086,'amu*angstrom^2'), symmetry=1, barrier=(23.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18695,0.12995,-0.000237258,2.19678e-07,-7.77954e-11,-78299.4,36.5986], Tmin=(100,'K'), Tmax=(836.8,'K')), NASAPolynomial(coeffs=[10.9567,0.0410097,-2.24523e-05,4.42862e-09,-3.07665e-13,-79250.2,-13.3613], Tmin=(836.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-652.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[CH][C]=C(F)C(=C(F)F)C(F)F(10545)',
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
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,250,446,589,854,899,182,240,577,636,1210,1413,234,589,736,816,1240,3237,1685,370,211.526,2432.42],'cm^-1')),
        HinderedRotor(inertia=(0.229147,'amu*angstrom^2'), symmetry=1, barrier=(7.2875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05396,'amu*angstrom^2'), symmetry=1, barrier=(33.5215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05848,'amu*angstrom^2'), symmetry=1, barrier=(33.5225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.755901,0.115647,-0.000175223,1.35511e-07,-4.0026e-11,-92204.9,34.2917], Tmin=(100,'K'), Tmax=(643.452,'K')), NASAPolynomial(coeffs=[13.8238,0.0370854,-2.02267e-05,4.08211e-09,-2.91592e-13,-94331.1,-31.561], Tmin=(643.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(F)=C(F)C([CH]F)=C(F)F(10788)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.968314,0.118071,-0.000173225,1.3025e-07,-3.89824e-11,-96901,33.5486], Tmin=(100,'K'), Tmax=(817.571,'K')), NASAPolynomial(coeffs=[16.3328,0.0334246,-1.79231e-05,3.61316e-09,-2.58907e-13,-99730,-46.4346], Tmin=(817.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-807.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + radical(CsCdF1sH) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = 'F[C]=C(C(F)F)C(F)(F)[C]=CF(10789)',
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
        HinderedRotor(inertia=(0.973227,'amu*angstrom^2'), symmetry=1, barrier=(22.3764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973184,'amu*angstrom^2'), symmetry=1, barrier=(22.3754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973179,'amu*angstrom^2'), symmetry=1, barrier=(22.3753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25029,0.131365,-0.00024019,2.21856e-07,-7.83579e-11,-80354.4,36.932], Tmin=(100,'K'), Tmax=(837.105,'K')), NASAPolynomial(coeffs=[11.3345,0.0404782,-2.22261e-05,4.38633e-09,-3.04719e-13,-81383.9,-15.1091], Tmin=(837.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'F[C]C(=CF)C(F)(F)C(F)=CF(10790)',
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
        HinderedRotor(inertia=(0.83361,'amu*angstrom^2'), symmetry=1, barrier=(19.1663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43943,'amu*angstrom^2'), symmetry=1, barrier=(33.0953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833783,'amu*angstrom^2'), symmetry=1, barrier=(19.1703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28511,0.130838,-0.000234562,2.13401e-07,-7.46923e-11,-88312.9,34.6422], Tmin=(100,'K'), Tmax=(831.422,'K')), NASAPolynomial(coeffs=[12.2836,0.0391187,-2.13885e-05,4.22299e-09,-2.93859e-13,-89655.4,-22.8185], Tmin=(831.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-735.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(=C(F)F)C(F)F(10791)',
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
        HinderedRotor(inertia=(1.10123,'amu*angstrom^2'), symmetry=1, barrier=(25.3195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09596,'amu*angstrom^2'), symmetry=1, barrier=(25.1984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09487,'amu*angstrom^2'), symmetry=1, barrier=(25.1733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18973,0.128759,-0.000230917,2.11321e-07,-7.4503e-11,-81869.8,36.7808], Tmin=(100,'K'), Tmax=(827.248,'K')), NASAPolynomial(coeffs=[11.7355,0.0398796,-2.19224e-05,4.34453e-09,-3.03205e-13,-83105.6,-17.6688], Tmin=(827.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-682.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C([CH]F)=C(F)F(10792)',
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
        HinderedRotor(inertia=(1.37791,'amu*angstrom^2'), symmetry=1, barrier=(31.6808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3781,'amu*angstrom^2'), symmetry=1, barrier=(31.6852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36882,'amu*angstrom^2'), symmetry=1, barrier=(31.4719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1097,0.122725,-0.000202913,1.74689e-07,-5.9303e-11,-88399.8,36.4113], Tmin=(100,'K'), Tmax=(794.295,'K')), NASAPolynomial(coeffs=[14.0348,0.0359399,-1.91592e-05,3.78864e-09,-2.65747e-13,-90473.9,-31.0755], Tmin=(794.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCFF) + group(Cds-CdsHH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)(F)[C]=CF(9832)',
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
        HinderedRotor(inertia=(0.23723,'amu*angstrom^2'), symmetry=1, barrier=(5.45438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239513,'amu*angstrom^2'), symmetry=1, barrier=(5.50687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239492,'amu*angstrom^2'), symmetry=1, barrier=(5.50639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3401.78,'J/mol'), sigma=(5.80027,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=531.35 K, Pc=39.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.98008,0.121354,-0.000206418,1.84147e-07,-6.45615e-11,-74967.1,37.425], Tmin=(100,'K'), Tmax=(797.963,'K')), NASAPolynomial(coeffs=[12.2641,0.0389094,-2.12617e-05,4.24291e-09,-2.98898e-13,-76569.6,-20.2787], Tmin=(797.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-624.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = 'FC=[C]C(F)(F)[C]=CF(10063)',
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
        HinderedRotor(inertia=(1.16778,'amu*angstrom^2'), symmetry=1, barrier=(26.8495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1603,'amu*angstrom^2'), symmetry=1, barrier=(26.6776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353529,0.0917968,-0.000166638,1.55124e-07,-5.51465e-11,-29843.7,26.9482], Tmin=(100,'K'), Tmax=(841.735,'K')), NASAPolynomial(coeffs=[8.21733,0.0310153,-1.6602e-05,3.24928e-09,-2.24797e-13,-30338.2,-4.70857], Tmin=(841.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1C(F)(F)C(=CF)C1(F)F(10793)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0711915,0.0939175,-0.000110121,6.72319e-08,-1.68007e-11,-125840,26.3189], Tmin=(100,'K'), Tmax=(958.517,'K')), NASAPolynomial(coeffs=[14.0137,0.035735,-1.90722e-05,3.90732e-09,-2.84797e-13,-128513,-40.3555], Tmin=(958.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1047.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = 'F[CH][C]1C(F)C(=C(F)F)C1(F)F(10758)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.514304,0.116261,-0.000212317,2.08448e-07,-7.85681e-11,-89800,31.8429], Tmin=(100,'K'), Tmax=(817.185,'K')), NASAPolynomial(coeffs=[5.5567,0.0513171,-2.84444e-05,5.6912e-09,-4.00586e-13,-89616,10.9764], Tmin=(817.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-747.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCCFF) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = 'F[C]C(=C(F)F)C(F)(F)C=CF(10794)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04029,0.123833,-0.000215509,1.94944e-07,-6.87281e-11,-92420.2,34.7242], Tmin=(100,'K'), Tmax=(810.879,'K')), NASAPolynomial(coeffs=[11.8412,0.0396936,-2.17651e-05,4.33387e-09,-3.0429e-13,-93832.2,-20.5459], Tmin=(810.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-769.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=[C]C(F)(F)C(=CF)C(F)F(10795)',
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
        HinderedRotor(inertia=(1.00469,'amu*angstrom^2'), symmetry=1, barrier=(23.0997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00266,'amu*angstrom^2'), symmetry=1, barrier=(23.0531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00083,'amu*angstrom^2'), symmetry=1, barrier=(23.0111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20766,0.130168,-0.000236988,2.18676e-07,-7.72468e-11,-80395.3,36.7804], Tmin=(100,'K'), Tmax=(836.061,'K')), NASAPolynomial(coeffs=[11.2472,0.0405114,-2.21855e-05,4.37717e-09,-3.04173e-13,-81427,-14.7923], Tmin=(836.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=C(C(F)(F)F)C(F)(F)[C]=CF(10796)',
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
        HinderedRotor(inertia=(1.02321,'amu*angstrom^2'), symmetry=1, barrier=(23.5257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02391,'amu*angstrom^2'), symmetry=1, barrier=(23.5416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02315,'amu*angstrom^2'), symmetry=1, barrier=(23.5242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19971,0.126829,-0.000220286,1.95644e-07,-6.74682e-11,-87598.9,35.4166], Tmin=(100,'K'), Tmax=(820.876,'K')), NASAPolynomial(coeffs=[13.3261,0.0365742,-1.97784e-05,3.90469e-09,-2.72372e-13,-89327.6,-27.7984], Tmin=(820.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C(F)F)C(F)(F)C(F)=CF(10797)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26505,0.126974,-0.000205685,1.78191e-07,-6.09335e-11,-91367.7,35.4238], Tmin=(100,'K'), Tmax=(813.226,'K')), NASAPolynomial(coeffs=[12.7105,0.0428063,-2.19842e-05,4.27055e-09,-2.96432e-13,-93130.7,-25.9744], Tmin=(813.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-761.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(=CF)C(F)(F)F(10798)',
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
        HinderedRotor(inertia=(1.04818,'amu*angstrom^2'), symmetry=1, barrier=(24.0998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04961,'amu*angstrom^2'), symmetry=1, barrier=(24.1326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05149,'amu*angstrom^2'), symmetry=1, barrier=(24.1758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1999,0.126832,-0.000220296,1.9566e-07,-6.74763e-11,-87598.9,35.4173], Tmin=(100,'K'), Tmax=(820.83,'K')), NASAPolynomial(coeffs=[13.3264,0.0365736,-1.9778e-05,3.9046e-09,-2.72363e-13,-89327.7,-27.8004], Tmin=(820.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['FC=C=C(F)F(5101)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(27.3892,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 27.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH]C(=C(F)F)C([CH]F)=C(F)F(9827)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C(F)[C]=C(F)F(9831)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH][C]=C(F)C(=CF)C(F)(F)F(10778)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(282.564,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CF-2(3543)', 'F[CH]C([C](F)F)=C(F)F(9960)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(137)', 'FC=[C]C(F)(F)[C]=C(F)F(10779)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['FC=C1C(F)C(=C(F)F)C1(F)F(10770)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_1H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['FC=[C]C(F)(F)[C]1C(F)C1(F)F(10780)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH][C]1C(F)(F)C(=CF)C1(F)F(10781)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH]C1([C](F)F)C(=CF)C1(F)F(10782)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(98.4045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 98.0 to 98.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'F[CH]C(C(F)=C=CF)=C(F)F(10783)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.5389,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC=C=C(F)F(5101)', 'F[CH][C]=C(F)F(5078)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CC(F)(F)C([CH]F)=C(F)F(10784)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(60.8867,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'FC#CC(F)(F)C([CH]F)=C(F)F(10785)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C]=C(F)F(5078)', 'F[CH][C]=C(F)F(5078)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -3.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHF(40)', 'FC=[C]C(F)(F)[C]=C(F)F(10779)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C]=CC(F)(F)C([CH]F)=C(F)F(10786)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C]=[C]C(F)(F)C(CF)=C(F)F(10787)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH][C]=C(F)C(=C(F)F)C(F)F(10545)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(221.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH]C(F)=C(F)C([CH]F)=C(F)F(10788)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(158.255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]=C(C(F)F)C(F)(F)[C]=CF(10789)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]C(=CF)C(F)(F)C(F)=CF(10790)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(231.15,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]C(F)(F)C(=C(F)F)C(F)F(10791)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(162.648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(F)C(F)(F)C([CH]F)=C(F)F(10792)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(172.145,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]C(F)(F)C(F)(F)[C]=CF(9832)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.14642e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]F(156)', 'FC=[C]C(F)(F)[C]=CF(10063)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['FC=C1C(F)(F)C(=CF)C1(F)F(10793)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_noH]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[CH][C]1C(F)C(=C(F)F)C1(F)F(10758)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CF2(43)', 'FC=[C]C(F)(F)[C]=CF(10063)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0195237,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['F[C]C(=C(F)F)C(F)(F)C=CF(10794)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.07654e+10,'s^-1'), n=1.20849, Ea=(247.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_Cd;XH_out] + [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]=[C]C(F)(F)C(=CF)C(F)F(10795)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C(F)(F)F)C(F)(F)[C]=CF(10796)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    products = ['[CH]C(=C(F)F)C(F)(F)C(F)=CF(10797)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(240.488,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]C(F)(F)C(=CF)C(F)(F)F(10798)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(184.741,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #2746',
    isomers = [
        'F[CH]C(=C(F)F)C(F)(F)[C]=CF(9829)',
    ],
    reactants = [
        ('FC=C=C(F)F(5101)', 'FC=C=C(F)F(5101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2746',
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

