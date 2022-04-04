species(
    label = '[CH]=CC[C]=CC(1974)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,D} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u1 p0 c0 {1,S} {3,D}
6  C u1 p0 c0 {4,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (537.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,192.457],'cm^-1')),
        HinderedRotor(inertia=(0.399238,'amu*angstrom^2'), symmetry=1, barrier=(10.4978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398787,'amu*angstrom^2'), symmetry=1, barrier=(10.4982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399153,'amu*angstrom^2'), symmetry=1, barrier=(10.4981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59138,0.0504568,-3.32835e-05,1.1347e-08,-1.61133e-12,64754.2,24.1518], Tmin=(100,'K'), Tmax=(1581.26,'K')), NASAPolynomial(coeffs=[10.8185,0.0271158,-1.11419e-05,2.01202e-09,-1.35452e-13,61836.1,-24.5916], Tmin=(1581.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(22)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C=CC(1692)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.759585,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74634,0.0218191,8.22302e-06,-2.14762e-08,8.55596e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[6.82084,0.0192337,-7.45616e-06,1.36535e-09,-9.53184e-14,16028,-10.4336], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(24)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=CC[C]=C(1415)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {3,D}
5  C u1 p0 c0 {2,D} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (573.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,260.969],'cm^-1')),
        HinderedRotor(inertia=(0.2644,'amu*angstrom^2'), symmetry=1, barrier=(12.7891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26239,'amu*angstrom^2'), symmetry=1, barrier=(12.842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.101,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3053.36,'J/mol'), sigma=(5.32914,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.93 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08026,0.0367088,-1.90945e-05,1.3615e-09,1.39993e-12,69072.4,20.5946], Tmin=(100,'K'), Tmax=(1131.67,'K')), NASAPolynomial(coeffs=[9.36813,0.0195593,-7.77598e-06,1.41684e-09,-9.75032e-14,66871.6,-17.9026], Tmin=(1131.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])=CC(2403)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {3,D} {4,S} {5,S}
3  C u0 p0 c0 {1,S} {2,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u1 p0 c0 {4,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (418.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09672,0.0537897,-2.27408e-05,-1.19719e-08,9.36004e-12,50495.6,21.5059], Tmin=(100,'K'), Tmax=(978.95,'K')), NASAPolynomial(coeffs=[14.6624,0.0210868,-7.45452e-06,1.33248e-09,-9.36584e-14,46750.6,-49.214], Tmin=(978.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=C(1759)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 C u2 p0 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,227.33,229.304,232.839],'cm^-1')),
        HinderedRotor(inertia=(1.31513,'amu*angstrom^2'), symmetry=1, barrier=(50.4662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817973,3.34731e-05,-4.36187e-08,1.5821e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.759,'K')), NASAPolynomial(coeffs=[5.3676,0.0170742,-6.35103e-06,1.16618e-09,-8.27612e-14,44095,-3.44632], Tmin=(983.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=CC(210)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144723,-2.91516e-06,-2.04629e-09,7.91531e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.87,'K')), NASAPolynomial(coeffs=[5.67826,0.0121697,-4.94654e-06,9.00478e-10,-6.07647e-14,66718.8,-3.96163], Tmin=(1455.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'CC=C1C=CC1(2508)',
    structure = adjacencyList("""1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  C u0 p0 c0 {1,S} {6,D} {13,S}
6  C u0 p0 c0 {3,S} {5,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (183.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9823,0.0252077,6.02928e-05,-9.87089e-08,4.02776e-11,22183.9,17.5609], Tmin=(100,'K'), Tmax=(949.585,'K')), NASAPolynomial(coeffs=[15.0549,0.0184602,-5.3754e-06,9.80312e-10,-7.56161e-14,17522.7,-56.3014], Tmin=(949.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene)"""),
)

species(
    label = 'C=CC=C=CC(2509)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,D} {10,S}
3  C u0 p0 c0 {4,S} {5,D} {11,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  C u0 p0 c0 {2,D} {4,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (198.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37018,0.0484067,-1.51573e-05,-1.43602e-08,8.83648e-12,23973.9,20.1002], Tmin=(100,'K'), Tmax=(1007.76,'K')), NASAPolynomial(coeffs=[12.6957,0.0238531,-8.97447e-06,1.63645e-09,-1.14923e-13,20655.4,-39.7658], Tmin=(1007.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=[CH](58)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1227.31,2941.71],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843434,2.04094e-05,-2.51944e-08,9.47575e-12,71059.3,4.72115], Tmin=(100,'K'), Tmax=(935.367,'K')), NASAPolynomial(coeffs=[4.92141,0.00358272,-9.24433e-07,1.57271e-10,-1.19539e-14,70476.2,-2.30656], Tmin=(935.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'H(6)',
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
    label = '[CH]=CC=C=CC(2510)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,D} {10,S}
3  C u0 p0 c0 {4,S} {5,D} {12,S}
4  C u0 p0 c0 {3,S} {6,D} {11,S}
5  C u0 p0 c0 {2,D} {3,D}
6  C u1 p0 c0 {4,D} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (445.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.01797,'amu*angstrom^2'), symmetry=1, barrier=(23.4051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00853,'amu*angstrom^2'), symmetry=1, barrier=(23.188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28469,0.0520656,-3.22782e-05,3.59747e-09,2.68914e-12,53694.2,20.8103], Tmin=(100,'K'), Tmax=(1028.48,'K')), NASAPolynomial(coeffs=[13.0101,0.0208964,-7.86996e-06,1.42119e-09,-9.86978e-14,50518.9,-39.7984], Tmin=(1028.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = 'CH3(18)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CCC#C(1475)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,S} {5,T}
4  C u1 p0 c0 {2,D} {9,S}
5  C u0 p0 c0 {3,T} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (505.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.02162,'amu*angstrom^2'), symmetry=1, barrier=(23.4889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02127,'amu*angstrom^2'), symmetry=1, barrier=(23.4809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0931,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3067.07,'J/mol'), sigma=(5.28401,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=479.07 K, Pc=47.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17298,0.0348445,-1.79012e-05,-1.32328e-09,2.99839e-12,60835.1,16.9512], Tmin=(100,'K'), Tmax=(1020.6,'K')), NASAPolynomial(coeffs=[9.61241,0.0162441,-6.07927e-06,1.0894e-09,-7.51977e-14,58766.7,-21.7849], Tmin=(1020.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC#CC(2511)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,D} {12,S}
4  C u0 p0 c0 {1,S} {5,T}
5  C u0 p0 c0 {2,S} {4,T}
6  C u1 p0 c0 {3,D} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (463.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2100,2250,500,550,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0674026,'amu*angstrom^2'), symmetry=1, barrier=(9.62059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690426,'amu*angstrom^2'), symmetry=1, barrier=(15.8743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792267,'amu*angstrom^2'), symmetry=1, barrier=(73.1611,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77366,0.0428163,-1.6596e-05,-6.55821e-09,5.28504e-12,55777,21.8045], Tmin=(100,'K'), Tmax=(991.975,'K')), NASAPolynomial(coeffs=[9.51841,0.0246549,-8.89448e-06,1.54647e-09,-1.04517e-13,53597.5,-18.7379], Tmin=(991.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC(1693)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u1 p0 c0 {4,S} {9,S} {10,S}
4  C u1 p0 c0 {2,D} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12014e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.76,'K')), NASAPolynomial(coeffs=[10.7464,0.0143241,-5.20138e-06,8.69083e-10,-5.48388e-14,40045.6,-31.3798], Tmin=(2065.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C#CC[C]=CC(2512)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,D} {12,S}
4  C u1 p0 c0 {1,S} {3,D}
5  C u0 p0 c0 {1,S} {6,T}
6  C u0 p0 c0 {5,T} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (459.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,2175,525,750,770,3400,2100,210.641],'cm^-1')),
        HinderedRotor(inertia=(0.350755,'amu*angstrom^2'), symmetry=1, barrier=(11.1108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352512,'amu*angstrom^2'), symmetry=1, barrier=(11.1115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716715,'amu*angstrom^2'), symmetry=1, barrier=(22.5366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83682,0.0483386,-3.5589e-05,1.46638e-08,-2.6197e-12,55395.8,20.0329], Tmin=(100,'K'), Tmax=(1259.84,'K')), NASAPolynomial(coeffs=[8.04865,0.0286159,-1.21064e-05,2.23758e-09,-1.53839e-13,53830.6,-11.3703], Tmin=(1259.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = '[CH]C=CC=CC(1977)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u2 p0 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (392.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.11654,'amu*angstrom^2'), symmetry=1, barrier=(48.6633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11652,'amu*angstrom^2'), symmetry=1, barrier=(48.6629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11654,'amu*angstrom^2'), symmetry=1, barrier=(48.6633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32559,0.0487423,-5.87323e-06,-2.32556e-08,1.14131e-11,47318.5,22.3158], Tmin=(100,'K'), Tmax=(1013.51,'K')), NASAPolynomial(coeffs=[11.4754,0.0307363,-1.18614e-05,2.15145e-09,-1.49452e-13,44128.5,-32.3751], Tmin=(1013.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CCC=[C]C(1594)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,D} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u1 p0 c0 {2,S} {3,D}
6  C u1 p0 c0 {4,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (537.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,192.457],'cm^-1')),
        HinderedRotor(inertia=(0.399238,'amu*angstrom^2'), symmetry=1, barrier=(10.4978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398787,'amu*angstrom^2'), symmetry=1, barrier=(10.4982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399153,'amu*angstrom^2'), symmetry=1, barrier=(10.4981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59138,0.0504568,-3.32835e-05,1.1347e-08,-1.61133e-12,64754.2,24.1518], Tmin=(100,'K'), Tmax=(1581.26,'K')), NASAPolynomial(coeffs=[10.8185,0.0271158,-1.11419e-05,2.01202e-09,-1.35452e-13,61836.1,-24.5916], Tmin=(1581.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C[C]=CC(1593)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,D} {12,S}
4  C u0 p0 c0 {6,D} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {3,D}
6  C u1 p0 c0 {1,S} {4,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (528.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,273.549,275.402],'cm^-1')),
        HinderedRotor(inertia=(0.163359,'amu*angstrom^2'), symmetry=1, barrier=(8.76913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166454,'amu*angstrom^2'), symmetry=1, barrier=(8.76245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166451,'amu*angstrom^2'), symmetry=1, barrier=(8.77286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90702,0.0484235,-3.11216e-05,1.0655e-08,-1.59398e-12,63625.6,23.0809], Tmin=(100,'K'), Tmax=(1425.91,'K')), NASAPolynomial(coeffs=[7.76642,0.0319864,-1.38302e-05,2.57047e-09,-1.76534e-13,61954.7,-7.2662], Tmin=(1425.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC=CC(1596)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {2,S} {3,D} {13,S}
5  C u1 p0 c0 {1,S} {6,D}
6  C u1 p0 c0 {5,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (537.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,192.457],'cm^-1')),
        HinderedRotor(inertia=(0.399238,'amu*angstrom^2'), symmetry=1, barrier=(10.4978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398787,'amu*angstrom^2'), symmetry=1, barrier=(10.4982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399153,'amu*angstrom^2'), symmetry=1, barrier=(10.4981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59138,0.0504568,-3.32835e-05,1.1347e-08,-1.61133e-12,64754.2,24.1518], Tmin=(100,'K'), Tmax=(1581.26,'K')), NASAPolynomial(coeffs=[10.8185,0.0271158,-1.11419e-05,2.01202e-09,-1.35452e-13,61836.1,-24.5916], Tmin=(1581.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC=C[CH2](1839)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {6,D} {11,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  C u1 p0 c0 {4,D} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (451.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.781166,'amu*angstrom^2'), symmetry=1, barrier=(17.9605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0276554,'amu*angstrom^2'), symmetry=1, barrier=(17.9095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781571,'amu*angstrom^2'), symmetry=1, barrier=(17.9699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3289.69,'J/mol'), sigma=(5.72141,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.84 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4834,0.0450322,-5.9668e-06,-2.29893e-08,1.15542e-11,54380.8,24.3471], Tmin=(100,'K'), Tmax=(1010.45,'K')), NASAPolynomial(coeffs=[12.4919,0.0243955,-9.38904e-06,1.73866e-09,-1.23306e-13,50984.9,-34.6724], Tmin=(1010.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C[C]=CC(2513)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,D} {10,S}
3  C u0 p0 c0 {4,D} {5,S} {11,S}
4  C u0 p0 c0 {3,D} {6,S} {12,S}
5  C u1 p0 c0 {3,S} {13,S} {14,S}
6  C u1 p0 c0 {2,D} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (338.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53745,0.0463667,-1.1901e-05,-1.57978e-08,9.31935e-12,40858.4,21.6], Tmin=(100,'K'), Tmax=(968.726,'K')), NASAPolynomial(coeffs=[10.6019,0.0268685,-9.47277e-06,1.63757e-09,-1.11014e-13,38260.9,-26.1849], Tmin=(968.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC[C]=[C]C(1595)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {6,D}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (528.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,273.549,275.402],'cm^-1')),
        HinderedRotor(inertia=(0.163359,'amu*angstrom^2'), symmetry=1, barrier=(8.76913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166454,'amu*angstrom^2'), symmetry=1, barrier=(8.76245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166451,'amu*angstrom^2'), symmetry=1, barrier=(8.77286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90702,0.0484235,-3.11216e-05,1.0655e-08,-1.59398e-12,63625.6,23.0809], Tmin=(100,'K'), Tmax=(1425.91,'K')), NASAPolynomial(coeffs=[7.76642,0.0319864,-1.38302e-05,2.57047e-09,-1.76534e-13,61954.7,-7.2662], Tmin=(1425.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC=C(1976)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {5,S} {6,D} {10,S}
4  C u0 p0 c0 {2,D} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {3,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (442.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.1276,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50576,0.0463,-1.46267e-05,-1.06408e-08,6.41678e-12,53265.3,24.339], Tmin=(100,'K'), Tmax=(1059.35,'K')), NASAPolynomial(coeffs=[11.3117,0.0263196,-1.04712e-05,1.9333e-09,-1.3516e-13,50231.2,-28.0486], Tmin=(1059.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
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
    E0 = (127.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (582.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (297.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (530.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (135.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (205.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (325.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (262.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (267.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (264.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (209.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (275.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (541.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (244.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (361.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (232.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (302.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (288.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (319.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (232.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (259.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['C2H2(22)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH]=CC[C]=C(1415)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.707e+08,'m^3/(mol*s)'), n=-0.589251, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_Ext-7R!H-R_4CbCdCt->Cd',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_Ext-7R!H-R_4CbCdCt->Cd
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['[CH]=CC([CH2])=CC(2403)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C=C(1759)', '[C]=CC(210)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['CC=C1C=CC1(2508)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['C=CC=C=CC(2509)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[CH](58)', 'C=C=CC(1692)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.043549,'m^3/(mol*s)'), n=2.3919, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5999870202516084, var=1.7833843885607656, Tref=1000.0, N=11, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_Ext-5R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_Ext-5R!H-R_N-5R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', '[CH]=CC=C=CC(2510)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1306.26,'m^3/(mol*s)'), n=1.34286, Ea=(15.5865,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2935929756704599, var=0.2435564102543282, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(18)', '[CH]=CCC#C(1475)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(906.033,'m^3/(mol*s)'), n=1.26419, Ea=(36.4431,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.008738971226930257, var=4.370254686200645, Tref=1000.0, N=14, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_3C-u1_N-1R!H->S_N-2R!H->S_Ext-2CNO-R_Sp-4R!H-2CNO_N-2CNO-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_3R->C_3C-u1_N-1R!H->S_N-2R!H->S_Ext-2CNO-R_Sp-4R!H-2CNO_N-2CNO-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[CH]=CCC#CC(2511)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2H2(22)', '[CH2][C]=CC(1693)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.000545692,'m^3/(mol*s)'), n=2.99187, Ea=(41.7092,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.13476088593216162, var=1.8284954697859765, Tref=1000.0, N=58, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_Sp-4R!H-3R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(6)', 'C#CC[C]=CC(2512)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1108.07,'m^3/(mol*s)'), n=1.405, Ea=(14.5573,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[CH](58)', '[CH2][C]=CC(1693)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.63958e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C=CC=CC(1977)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.38e+10,'s^-1'), n=0.71, Ea=(262.755,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 155 used for R2H_S;C_rad_out_H/Cd;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CCC=[C]C(1594)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['C=[C]C[C]=CC(1593)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]CC=CC(1596)'],
    products = ['[CH]=CC[C]=CC(1974)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['[CH]=CCC=C[CH2](1839)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['[CH2]C=C[C]=CC(2513)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['C=CC[C]=[C]C(1595)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC[C]=CC(1974)'],
    products = ['[CH2]C=[C]CC=C(1976)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1193',
    isomers = [
        '[CH]=CC[C]=CC(1974)',
    ],
    reactants = [
        ('C2H2(22)', 'C=C=CC(1692)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1193',
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

