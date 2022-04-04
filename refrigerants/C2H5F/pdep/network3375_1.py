species(
    label = '[CH2]C(=CF)C[C]=CF(8639)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (55.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,615,860,1140,1343,3152,1685,370,180,1053.43,3718.88],'cm^-1')),
        HinderedRotor(inertia=(0.692699,'amu*angstrom^2'), symmetry=1, barrier=(15.9265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.71331,'amu*angstrom^2'), symmetry=1, barrier=(85.3764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22007,'amu*angstrom^2'), symmetry=1, barrier=(28.0518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79985,0.0695946,-6.58924e-05,3.30466e-08,-6.75433e-12,6785.45,27.5323], Tmin=(100,'K'), Tmax=(1166.72,'K')), NASAPolynomial(coeffs=[13.0241,0.0276853,-1.20124e-05,2.25992e-09,-1.57586e-13,3932.94,-33.3281], Tmin=(1166.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CF(5941)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(=CF)C([CH2])=CF(8637)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u0 p0 c0 {3,S} {6,S} {8,D}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {3,D} {9,S}
8  C u0 p0 c0 {2,S} {4,D} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-61.8434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150568,0.074754,-6.16217e-05,1.75068e-08,1.10278e-12,-7290.62,24.5898], Tmin=(100,'K'), Tmax=(985.1,'K')), NASAPolynomial(coeffs=[19.1504,0.0180221,-6.32504e-06,1.12386e-09,-7.88878e-14,-12024.6,-71.8163], Tmin=(985.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.8434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCFH) + group(CdCFH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'FC=[C]CC[C]=CF(8642)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,D} {13,S}
6  C u0 p0 c0 {2,S} {8,D} {14,S}
7  C u1 p0 c0 {3,S} {5,D}
8  C u1 p0 c0 {4,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (173.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,535,695,819,901,1081,1199,1292,1394,3088,3216,1670,1700,300,440,319.376,319.384,319.384,319.389,1799.05],'cm^-1')),
        HinderedRotor(inertia=(0.121065,'amu*angstrom^2'), symmetry=1, barrier=(8.76333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354067,'amu*angstrom^2'), symmetry=1, barrier=(25.6312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00225745,'amu*angstrom^2'), symmetry=1, barrier=(25.6312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3359.23,'J/mol'), sigma=(5.56579,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=524.70 K, Pc=44.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18614,0.0606823,-4.80832e-05,1.97439e-08,-3.35085e-12,20979.4,30.0822], Tmin=(100,'K'), Tmax=(1363.08,'K')), NASAPolynomial(coeffs=[12.1841,0.0284086,-1.25679e-05,2.37386e-09,-1.65061e-13,17981.2,-26.3833], Tmin=(1363.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cdj(Cs-CsHH)(Cd-F1sH)) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[C]=CF(1334)',
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
    label = '[CH2]C([CH2])=CF(6739)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,D}
3  C u1 p0 c0 {2,S} {6,S} {7,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (78.7516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.00398153,'amu*angstrom^2'), symmetry=1, barrier=(1.08133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20218,'amu*angstrom^2'), symmetry=1, barrier=(54.9214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13506,0.0345568,-1.31986e-05,-7.34369e-09,5.24512e-12,9544.42,16.2307], Tmin=(100,'K'), Tmax=(1016.69,'K')), NASAPolynomial(coeffs=[10.3402,0.0159379,-5.88635e-06,1.07937e-09,-7.62842e-14,7169.86,-26.9629], Tmin=(1016.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.7516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC=[C]C[C]=CF(9020)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {7,D} {11,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u1 p0 c0 {3,S} {5,D}
7  C u1 p0 c0 {3,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (180.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,535,695,819,901,1081,1199,1292,1394,3088,3216,1670,1700,300,440,263.628,264.273,265.358,2068.93],'cm^-1')),
        HinderedRotor(inertia=(0.203584,'amu*angstrom^2'), symmetry=1, barrier=(10.1673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214999,'amu*angstrom^2'), symmetry=1, barrier=(10.1565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59242,0.056792,-6.86453e-05,4.8239e-08,-1.41385e-11,21834.7,23.4394], Tmin=(100,'K'), Tmax=(821.247,'K')), NASAPolynomial(coeffs=[7.95199,0.02582,-1.20809e-05,2.32618e-09,-1.63319e-13,20790,-5.9901], Tmin=(821.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1CC(=CF)C1(9509)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u0 p0 c0 {1,S} {5,D} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-187.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48992,0.0372005,3.54172e-05,-7.49959e-08,3.17341e-11,-22437.9,23.0486], Tmin=(100,'K'), Tmax=(969.476,'K')), NASAPolynomial(coeffs=[16.2833,0.0204332,-7.13486e-06,1.36629e-09,-1.03518e-13,-27386.7,-58.592], Tmin=(969.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = 'CC(C=C=CF)=CF(9510)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {8,D} {12,S}
6  C u0 p0 c0 {1,S} {4,D} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u0 p0 c0 {5,D} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-168.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635389,0.0710497,-6.81927e-05,3.41445e-08,-6.86818e-12,-20158.9,23.8967], Tmin=(100,'K'), Tmax=(1196.21,'K')), NASAPolynomial(coeffs=[14.4269,0.0249323,-1.03631e-05,1.91505e-09,-1.32404e-13,-23458.4,-45.1101], Tmin=(1196.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'FC=[C]C[C]1CC1F(9511)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (157.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98749,0.0493131,-2.75475e-05,7.03303e-09,-7.16184e-13,19043.3,26.6407], Tmin=(100,'K'), Tmax=(2097.31,'K')), NASAPolynomial(coeffs=[12.8749,0.0285486,-1.26967e-05,2.31245e-09,-1.5349e-13,14476.4,-33.9485], Tmin=(2097.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[CH2][C]1CC(=CF)C1F(9512)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {3,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (81.7758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42147,0.0618579,-6.96875e-05,5.78458e-08,-2.16104e-11,9923.31,24.9421], Tmin=(100,'K'), Tmax=(741.691,'K')), NASAPolynomial(coeffs=[4.0764,0.0424124,-1.99913e-05,3.85598e-09,-2.70523e-13,9670.51,13.8777], Tmin=(741.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.7758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([CH]F)CC1=CF(9513)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {3,S} {12,S} {13,S}
7  C u1 p0 c0 {1,S} {3,S} {11,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (153.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539514,0.0679891,-5.72917e-05,2.40323e-08,-4.00653e-12,18581.7,28.8609], Tmin=(100,'K'), Tmax=(1435.08,'K')), NASAPolynomial(coeffs=[17.0339,0.0220141,-9.23685e-06,1.7084e-09,-1.17548e-13,13847.5,-56.6731], Tmin=(1435.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cd(Cd-F)-Cs-Cs) + radical(Neopentyl) + radical(CsCsF1sH)"""),
)

species(
    label = 'C=[C][CH]F(5942)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (196.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.900551,'amu*angstrom^2'), symmetry=1, barrier=(20.7054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79302e-05,1.96926e-08,-6.19852e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.897,'K')), NASAPolynomial(coeffs=[4.81186,0.0161179,-7.47821e-06,1.46101e-09,-1.03887e-13,23340.4,3.53844], Tmin=(747.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = '[CH2]C(C=C=CF)=CF(9514)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {3,S} {8,D} {9,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {3,D} {10,S}
7  C u0 p0 c0 {2,S} {8,D} {13,S}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-17.1396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,113,247,382,1207,3490,540,610,2055,229.109,697.766],'cm^-1')),
        HinderedRotor(inertia=(0.061566,'amu*angstrom^2'), symmetry=1, barrier=(21.2671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101202,'amu*angstrom^2'), symmetry=1, barrier=(34.966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355736,0.0718025,-7.07592e-05,3.48534e-08,-6.70725e-12,-1923.17,25.1833], Tmin=(100,'K'), Tmax=(1271.22,'K')), NASAPolynomial(coeffs=[17.6696,0.0173225,-6.47412e-06,1.14009e-09,-7.71011e-14,-6325.08,-62.5009], Tmin=(1271.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
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
    label = 'C#CCC([CH2])=CF(8632)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (178.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0179465,'amu*angstrom^2'), symmetry=1, barrier=(19.9507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56741,'amu*angstrom^2'), symmetry=1, barrier=(27.129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.55537,'amu*angstrom^2'), symmetry=1, barrier=(81.7449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901234,0.0610971,-5.33578e-05,2.41141e-08,-4.34938e-12,21620.1,23.0087], Tmin=(100,'K'), Tmax=(1335.59,'K')), NASAPolynomial(coeffs=[14.4496,0.0205211,-7.78728e-06,1.36749e-09,-9.16354e-14,18001.1,-46.2747], Tmin=(1335.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CF)CC#CF(9515)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {2,S} {7,T}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (79.0176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,2175,525,239,401,1367,1457.3,2698.59],'cm^-1')),
        HinderedRotor(inertia=(0.852017,'amu*angstrom^2'), symmetry=1, barrier=(19.5896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852073,'amu*angstrom^2'), symmetry=1, barrier=(19.5908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209422,'amu*angstrom^2'), symmetry=1, barrier=(51.9742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17709,0.0642037,-6.20149e-05,3.24194e-08,-7.03653e-12,9603.52,24.3433], Tmin=(100,'K'), Tmax=(1089.19,'K')), NASAPolynomial(coeffs=[10.8271,0.0287646,-1.3209e-05,2.54639e-09,-1.79812e-13,7501.4,-23.0366], Tmin=(1089.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + group(Ct-CtCs) + group(CtCF) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]F)=CC=CF(9516)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,D} {6,S} {7,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {8,D} {9,S}
6  C u1 p0 c0 {1,S} {3,S} {11,S}
7  C u1 p0 c0 {3,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-84.6254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711105,0.062441,-3.77049e-05,1.3961e-09,4.67873e-12,-10051,27.5759], Tmin=(100,'K'), Tmax=(1013.82,'K')), NASAPolynomial(coeffs=[15.6846,0.0231253,-8.77354e-06,1.59778e-09,-1.12022e-13,-14102.7,-49.8766], Tmin=(1013.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.6254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[CH2]C(=CF)CC=[C]F(9517)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {7,D}
5  C u0 p0 c0 {3,S} {8,D} {11,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {4,D} {12,S}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (67.7301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,167,640,1190,390.278,390.442],'cm^-1')),
        HinderedRotor(inertia=(0.00110572,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718866,'amu*angstrom^2'), symmetry=1, barrier=(7.77434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310004,'amu*angstrom^2'), symmetry=1, barrier=(33.537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675547,0.0694927,-6.38763e-05,3.04521e-08,-5.845e-12,8268.73,27.9945], Tmin=(100,'K'), Tmax=(1247.88,'K')), NASAPolynomial(coeffs=[14.5216,0.0251098,-1.05259e-05,1.95001e-09,-1.34848e-13,4813.12,-41.8706], Tmin=(1247.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Allyl_P) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'CC([CH][C]=CF)=CF(9518)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,D} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (5.67975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703014,0.064547,-5.35312e-05,2.31025e-08,-4.00434e-12,808.431,27.2731], Tmin=(100,'K'), Tmax=(1380.99,'K')), NASAPolynomial(coeffs=[14.8903,0.023454,-8.89689e-06,1.5555e-09,-1.03689e-13,-3110.06,-45.7518], Tmin=(1380.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.67975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'CC(=[C]F)C[C]=CF(9519)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u1 p0 c0 {3,S} {6,D}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (154.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,615,860,1140,1343,3152,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.187499,'amu*angstrom^2'), symmetry=1, barrier=(4.31097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191005,'amu*angstrom^2'), symmetry=1, barrier=(4.39159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189453,'amu*angstrom^2'), symmetry=1, barrier=(4.3559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796711,0.0753098,-9.11945e-05,6.38231e-08,-1.85546e-11,18681.3,27.7639], Tmin=(100,'K'), Tmax=(829.142,'K')), NASAPolynomial(coeffs=[9.46896,0.0334741,-1.55125e-05,2.9738e-09,-2.0826e-13,17243.2,-12.4502], Tmin=(829.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'C=C([C]F)CC=CF(9520)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {5,D} {14,S}
8  C u2 p0 c0 {2,S} {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (52.9774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14476,0.0660357,-6.08327e-05,3.04224e-08,-6.40533e-12,6471.81,25.9559], Tmin=(100,'K'), Tmax=(1108.35,'K')), NASAPolynomial(coeffs=[10.4937,0.0322953,-1.51691e-05,2.95566e-09,-2.09845e-13,4399.45,-20.109], Tmin=(1108.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.9774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'CC(=CF)C[C]=[C]F(9521)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {1,S} {5,D} {14,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (154.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,194,682,905,1196,1383,3221,1685,370,167,640,1190,236.77,236.775],'cm^-1')),
        HinderedRotor(inertia=(0.1306,'amu*angstrom^2'), symmetry=1, barrier=(5.19586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130593,'amu*angstrom^2'), symmetry=1, barrier=(5.19589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130618,'amu*angstrom^2'), symmetry=1, barrier=(5.1958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850044,0.0739806,-8.74966e-05,5.99373e-08,-1.71083e-11,18639.9,27.5744], Tmin=(100,'K'), Tmax=(842.691,'K')), NASAPolynomial(coeffs=[9.37497,0.0335181,-1.54778e-05,2.96597e-09,-2.07819e-13,17203.1,-12.095], Tmin=(842.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=C(CF)C[C]=CF(9522)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {7,D} {13,S}
7  C u1 p0 c0 {3,S} {6,D}
8  C u1 p0 c0 {5,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (167.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,278.89,278.923,278.971],'cm^-1')),
        HinderedRotor(inertia=(0.00216757,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197759,'amu*angstrom^2'), symmetry=1, barrier=(10.9186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19771,'amu*angstrom^2'), symmetry=1, barrier=(10.9183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.90568,0.0731176,-8.27049e-05,5.34284e-08,-1.44644e-11,20233.6,28.0846], Tmin=(100,'K'), Tmax=(882.77,'K')), NASAPolynomial(coeffs=[9.55371,0.0339308,-1.61172e-05,3.14034e-09,-2.22539e-13,18706.8,-12.5586], Tmin=(882.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)CC(F)=CF(9523)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {5,D} {13,S}
8  C u2 p0 c0 {4,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (45.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776884,0.0715499,-6.10516e-05,2.83004e-08,-5.51592e-12,5528.86,27.1284], Tmin=(100,'K'), Tmax=(1193.4,'K')), NASAPolynomial(coeffs=[11.5495,0.0354425,-1.56677e-05,2.94762e-09,-2.04866e-13,2957.65,-26.7481], Tmin=(1193.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC(=CF)CF(9524)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {13,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u1 p0 c0 {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (167.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,194,682,905,1196,1383,3221,1685,370,3120,650,792.5,1650,267.65,267.683],'cm^-1')),
        HinderedRotor(inertia=(0.00235121,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232951,'amu*angstrom^2'), symmetry=1, barrier=(11.8526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233216,'amu*angstrom^2'), symmetry=1, barrier=(11.8539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905683,0.0731176,-8.27048e-05,5.34283e-08,-1.44644e-11,20233.6,28.0846], Tmin=(100,'K'), Tmax=(882.777,'K')), NASAPolynomial(coeffs=[9.55371,0.0339308,-1.61172e-05,3.14033e-09,-2.22539e-13,18706.8,-12.5586], Tmin=(882.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)CC([CH2])=CF(9525)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {7,D}
5  C u0 p0 c0 {1,S} {3,S} {8,D}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {4,D} {11,S}
8  C u1 p0 c0 {5,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (59.1652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,246,474,533,1155,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,3120,650,792.5,1650,180,1019.86],'cm^-1')),
        HinderedRotor(inertia=(0.0168087,'amu*angstrom^2'), symmetry=1, barrier=(12.1796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143508,'amu*angstrom^2'), symmetry=1, barrier=(3.29954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51868,'amu*angstrom^2'), symmetry=1, barrier=(34.9174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583689,0.071897,-6.83314e-05,3.3548e-08,-6.60964e-12,7241.63,27.7009], Tmin=(100,'K'), Tmax=(1219.74,'K')), NASAPolynomial(coeffs=[14.9111,0.0249113,-1.05495e-05,1.96621e-09,-1.36539e-13,3746.5,-44.2667], Tmin=(1219.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.1652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCFH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'C=[C]C(F)C[C]=CF(8641)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u0 p0 c0 {8,D} {13,S} {14,S}
7  C u1 p0 c0 {3,S} {5,D}
8  C u1 p0 c0 {4,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (154.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,164,312,561,654,898,1207,1299,3167,615,860,1140,1343,3152,2950,3100,1380,975,1025,1650,1670,1700,300,440,225.689,231.963,4000],'cm^-1')),
        HinderedRotor(inertia=(0.408779,'amu*angstrom^2'), symmetry=1, barrier=(14.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88772,'amu*angstrom^2'), symmetry=1, barrier=(33.2153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216039,'amu*angstrom^2'), symmetry=1, barrier=(8.24934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3299.44,'J/mol'), sigma=(5.52766,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=515.36 K, Pc=44.33 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03432,0.0710804,-7.94642e-05,5.22984e-08,-1.47164e-11,18686.5,29.1103], Tmin=(100,'K'), Tmax=(843.387,'K')), NASAPolynomial(coeffs=[8.35449,0.0363629,-1.77183e-05,3.49113e-09,-2.48898e-13,17451.7,-4.95847], Tmin=(843.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-CsHH)(Cd-F1sH)) + radical(Cdj(Cs-F1sCsH)(Cd-HH))"""),
)

species(
    label = '[CH]F(837)',
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
    label = 'C=[C]C[C]=CF(6834)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {6,D} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,D} {9,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {2,S} {3,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (372.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1670,1700,300,440,310.066,312.453,317.937],'cm^-1')),
        HinderedRotor(inertia=(0.0928139,'amu*angstrom^2'), symmetry=1, barrier=(6.32657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895553,'amu*angstrom^2'), symmetry=1, barrier=(6.33696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0622,0.0448898,-3.9659e-05,2.01535e-08,-4.42268e-12,44886.4,21.5872], Tmin=(100,'K'), Tmax=(1054.8,'K')), NASAPolynomial(coeffs=[7.3524,0.0248279,-1.11292e-05,2.12134e-09,-1.48779e-13,43770.4,-4.21715], Tmin=(1054.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C1CC(=CF)C1F(9526)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-198.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30833,0.0461722,1.73737e-06,-3.32978e-08,1.50108e-11,-23728.5,21.1446], Tmin=(100,'K'), Tmax=(1034.48,'K')), NASAPolynomial(coeffs=[14.2288,0.0254883,-1.0721e-05,2.08798e-09,-1.52122e-13,-27968.2,-49.1987], Tmin=(1034.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCCFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C(C=C=CF)CF(9527)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {8,D} {11,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u0 p0 c0 {5,D} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-152.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725757,0.0698044,-6.49514e-05,3.1288e-08,-6.08751e-12,-18203.6,24.3668], Tmin=(100,'K'), Tmax=(1227.63,'K')), NASAPolynomial(coeffs=[14.2147,0.0258532,-1.12487e-05,2.12463e-09,-1.48527e-13,-21515.5,-43.4759], Tmin=(1227.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'F[CH][C]1CC(=CF)C1(9528)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {1,S} {5,S} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (80.8886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54694,0.0591319,-6.0318e-05,4.7845e-08,-1.83743e-11,9812.18,24.6925], Tmin=(100,'K'), Tmax=(677.59,'K')), NASAPolynomial(coeffs=[3.6899,0.0440212,-2.14207e-05,4.21618e-09,-3.00218e-13,9578.24,15.6048], Tmin=(677.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.8886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH2]C(=C[C]=CF)CF(9529)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {8,S} {11,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-6.6607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865051,0.068331,-6.38767e-05,3.25289e-08,-6.81275e-12,-687.687,27.0995], Tmin=(100,'K'), Tmax=(1136.51,'K')), NASAPolynomial(coeffs=[11.9698,0.0292471,-1.2292e-05,2.26943e-09,-1.56454e-13,-3211.78,-27.8952], Tmin=(1136.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.6607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
)

species(
    label = '[CH]C(=CF)CC=CF(9530)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {4,D} {12,S}
7  C u0 p0 c0 {2,S} {5,D} {13,S}
8  C u2 p0 c0 {4,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (36.7974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722893,0.0687214,-5.37919e-05,2.21914e-08,-3.79177e-12,4546.26,28.1609], Tmin=(100,'K'), Tmax=(1359.76,'K')), NASAPolynomial(coeffs=[13.1034,0.0323016,-1.36158e-05,2.49375e-09,-1.70225e-13,1179.36,-35.3726], Tmin=(1359.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.7974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(CF)C[C]=[C]F(9531)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (170.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,167,640,1190,310.67,310.745,2122.63],'cm^-1')),
        HinderedRotor(inertia=(0.198105,'amu*angstrom^2'), symmetry=1, barrier=(13.4787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195446,'amu*angstrom^2'), symmetry=1, barrier=(13.4427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197183,'amu*angstrom^2'), symmetry=1, barrier=(13.4769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935619,0.0727485,-8.41309e-05,5.6745e-08,-1.61348e-11,20595.6,28.0649], Tmin=(100,'K'), Tmax=(840.974,'K')), NASAPolynomial(coeffs=[8.89108,0.0349091,-1.66384e-05,3.24123e-09,-2.29438e-13,19257.6,-8.9377], Tmin=(840.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=[C]CC(=C)C(F)F(9532)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u1 p0 c0 {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (140.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,421.884,421.884],'cm^-1')),
        HinderedRotor(inertia=(0.0877297,'amu*angstrom^2'), symmetry=1, barrier=(11.0805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0877297,'amu*angstrom^2'), symmetry=1, barrier=(11.0805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0877297,'amu*angstrom^2'), symmetry=1, barrier=(11.0805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995739,0.0708302,-7.52891e-05,4.50977e-08,-1.13621e-11,17022.5,28.0203], Tmin=(100,'K'), Tmax=(941.898,'K')), NASAPolynomial(coeffs=[9.76357,0.0335955,-1.59919e-05,3.12762e-09,-2.22402e-13,15370.8,-13.7545], Tmin=(941.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (-61.1483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (108.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (227.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (367.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (445.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-52.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (17.0988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (101.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (65.0872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (36.8032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (92.2687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (89.7641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (194.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (174.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (275.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (95.7517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (120.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (141.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (193.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (186.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (70.5106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (226.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (176.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (197.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (122.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (132.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (470.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-52.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (17.0988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (65.0872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (394.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (148.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (333.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (86.8011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (182.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['C=C=CF(5941)', 'C=C=CF(5941)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH2]C(=CF)C([CH2])=CF(8637)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]CC[C]=CF(8642)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.49684e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=CF(1334)', '[CH2]C([CH2])=CF(6739)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(18)', 'FC=[C]C[C]=CF(9020)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['FC=C1CC(=CF)C1(9509)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['CC(C=C=CF)=CF(9510)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['FC=[C]C[C]1CC1F(9511)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.57691e+12,'s^-1'), n=0.0709135, Ea=(162.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3_D;doublebond_intra;radadd_intra_cs2H] + [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH2][C]1CC(=CF)C1F(9512)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH2]C1([CH]F)CC1=CF(9513)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(97.9515,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 97.3 to 98.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C=CF(5941)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(3.28784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[CH2]C(C=C=CF)=CF(9514)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1306.26,'m^3/(mol*s)'), n=1.34286, Ea=(11.7012,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2935929756704599, var=0.2435564102543282, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CCC([CH2])=CF(8632)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.8932,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', '[CH2]C(=CF)CC#CF(9515)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][CH]F(5942)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH2]C([CH]F)=CC=CF(9516)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=CF)CC=[C]F(9517)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['CC([CH][C]=CF)=CF(9518)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC(=[C]F)C[C]=CF(9519)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['C=C([C]F)CC=CF(9520)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.07654e+10,'s^-1'), n=1.20849, Ea=(247.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_Cd;XH_out] + [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC(=CF)C[C]=[C]F(9521)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(CF)C[C]=CF(9522)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(175.382,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH]C(=C)CC(F)=CF(9523)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(237.318,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC(=CF)CF(9524)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(146.662,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(F)CC([CH2])=CF(9525)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(180.277,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(F)C[C]=CF(8641)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]F(837)', 'C=[C]C[C]=CF(6834)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['C=C1CC(=CF)C1F(9526)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_1H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['C=C(C=C=CF)CF(9527)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['F[CH][C]1CC(=CF)C1(9528)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CHF(40)', 'C=[C]C[C]=CF(6834)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH2]C(=C[C]=CF)CF(9529)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.23617e+10,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_1H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=CF)C[C]=CF(8639)'],
    products = ['[CH]C(=CF)CC=CF(9530)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.80239e+12,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C(CF)C[C]=[C]F(9531)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]CC(=C)C(F)F(9532)'],
    products = ['[CH2]C(=CF)C[C]=CF(8639)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(158.15,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3375',
    isomers = [
        '[CH2]C(=CF)C[C]=CF(8639)',
    ],
    reactants = [
        ('C=C=CF(5941)', 'C=C=CF(5941)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3375',
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

