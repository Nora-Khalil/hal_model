species(
    label = 'C[CH]CCC(392)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
"""),
    E0 = (26.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2618.62,2618.74],'cm^-1')),
        HinderedRotor(inertia=(0.0965287,'amu*angstrom^2'), symmetry=1, barrier=(8.14906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963935,'amu*angstrom^2'), symmetry=1, barrier=(8.14965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964049,'amu*angstrom^2'), symmetry=1, barrier=(8.14783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441321,'amu*angstrom^2'), symmetry=1, barrier=(37.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1409,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77405,0.0441717,-1.45099e-05,-1.38035e-09,1.20136e-12,3254.59,22.595], Tmin=(100,'K'), Tmax=(1461.3,'K')), NASAPolynomial(coeffs=[8.55197,0.0345018,-1.37025e-05,2.41137e-09,-1.59035e-13,325.215,-15.9208], Tmin=(1461.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = 'C2H5(31)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.89,1642.9,3622.3,3622.32],'cm^-1')),
        HinderedRotor(inertia=(0.866818,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=CC(60)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817973,3.34731e-05,-4.36187e-08,1.5821e-11,749.326,9.54029], Tmin=(100,'K'), Tmax=(983.759,'K')), NASAPolynomial(coeffs=[5.3676,0.0170742,-6.35103e-06,1.16618e-09,-8.27612e-14,-487.157,-4.54493], Tmin=(983.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[CH]CC(285)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (52.4618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0995105,'amu*angstrom^2'), symmetry=1, barrier=(2.28794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989853,'amu*angstrom^2'), symmetry=1, barrier=(2.27587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792956,'amu*angstrom^2'), symmetry=1, barrier=(18.2316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48182,0.0270986,1.03887e-05,-2.3665e-08,8.65276e-12,6369.62,17.2943], Tmin=(100,'K'), Tmax=(1083.85,'K')), NASAPolynomial(coeffs=[6.12698,0.0282086,-1.13017e-05,2.07331e-09,-1.43451e-13,4724.1,-4.5309], Tmin=(1083.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""butyl_2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C)CC(391)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
"""),
    E0 = (31.3037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,3369.01],'cm^-1')),
        HinderedRotor(inertia=(0.599971,'amu*angstrom^2'), symmetry=1, barrier=(13.7945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00710836,'amu*angstrom^2'), symmetry=1, barrier=(4.85209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.5828,'amu*angstrom^2'), symmetry=1, barrier=(82.3757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.599974,'amu*angstrom^2'), symmetry=1, barrier=(13.7946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1409,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65932,0.0428141,2.16854e-06,-2.73902e-08,1.23517e-11,3856.87,21.7514], Tmin=(100,'K'), Tmax=(991.73,'K')), NASAPolynomial(coeffs=[9.17442,0.0322671,-1.17724e-05,2.07625e-09,-1.4196e-13,1394.35,-19.3424], Tmin=(991.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.3037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(38)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u2 p0 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82364,-0.000909745,3.21388e-05,-3.73491e-08,1.33095e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.803,'K')), NASAPolynomial(coeffs=[4.3048,0.0094308,-3.27565e-06,5.95137e-10,-4.2732e-14,40709.2,1.84239], Tmin=(960.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CC(48)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928811,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0929139,'amu*angstrom^2'), symmetry=1, barrier=(2.13627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02813,0.0147026,2.405e-05,-3.66726e-08,1.38606e-11,10512.1,12.47], Tmin=(100,'K'), Tmax=(984.472,'K')), NASAPolynomial(coeffs=[6.16552,0.0184493,-6.79021e-06,1.23047e-09,-8.6385e-14,9095.02,-6.67654], Tmin=(984.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]CCC(398)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u2 p0 c0 {2,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (304.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,307.937,307.952,308.139,2853.43],'cm^-1')),
        HinderedRotor(inertia=(0.00177492,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00178101,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0887618,'amu*angstrom^2'), symmetry=1, barrier=(5.97765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09982,0.0375761,-2.10973e-05,5.90411e-09,-6.6831e-13,36657.6,16.6014], Tmin=(100,'K'), Tmax=(1987.27,'K')), NASAPolynomial(coeffs=[11.3193,0.0190193,-7.09065e-06,1.20535e-09,-7.72077e-14,32993.3,-34.2088], Tmin=(1987.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
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
    label = 'C[C]CCC(399)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u2 p0 c0 {2,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
"""),
    E0 = (280.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,274.501,274.659,3330.42],'cm^-1')),
        HinderedRotor(inertia=(0.117778,'amu*angstrom^2'), symmetry=1, barrier=(6.2872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290435,'amu*angstrom^2'), symmetry=1, barrier=(15.5403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.455107,'amu*angstrom^2'), symmetry=1, barrier=(24.2852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454112,'amu*angstrom^2'), symmetry=1, barrier=(24.2843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49126,0.0476291,-1.77615e-05,-4.7547e-09,3.68072e-12,33788.6,20.7635], Tmin=(100,'K'), Tmax=(1127.23,'K')), NASAPolynomial(coeffs=[10.1002,0.030205,-1.20403e-05,2.19072e-09,-1.50455e-13,31013.9,-25.4992], Tmin=(1127.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC=CCC(301)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-52.5474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,275.182],'cm^-1')),
        HinderedRotor(inertia=(0.171045,'amu*angstrom^2'), symmetry=1, barrier=(9.19243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170769,'amu*angstrom^2'), symmetry=1, barrier=(9.1925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170887,'amu*angstrom^2'), symmetry=1, barrier=(9.19314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94606,0.0364835,8.3545e-06,-2.83269e-08,1.1181e-11,-6238.66,19.2593], Tmin=(100,'K'), Tmax=(1059.79,'K')), NASAPolynomial(coeffs=[8.36592,0.0318147,-1.27245e-05,2.34964e-09,-1.63921e-13,-8697.95,-17.2685], Tmin=(1059.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.5474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCCC(288)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-40.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,434.525,434.657],'cm^-1')),
        HinderedRotor(inertia=(0.0737384,'amu*angstrom^2'), symmetry=1, barrier=(9.87175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737485,'amu*angstrom^2'), symmetry=1, barrier=(9.87065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0736821,'amu*angstrom^2'), symmetry=1, barrier=(9.87085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90148,0.0362538,1.29127e-05,-3.55761e-08,1.43065e-11,-4763.1,19.8564], Tmin=(100,'K'), Tmax=(1027.62,'K')), NASAPolynomial(coeffs=[9.28072,0.0304041,-1.19375e-05,2.20663e-09,-1.55066e-13,-7487.44,-21.8216], Tmin=(1027.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]C(195)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.2686e-06,-1.44876e-08,5.7495e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93907,0.0171892,-6.69154e-06,1.21547e-09,-8.39799e-14,33151.2,-4.14876], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[CH]C(55)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3121.1],'cm^-1')),
        HinderedRotor(inertia=(0.207862,'amu*angstrom^2'), symmetry=1, barrier=(4.77915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207796,'amu*angstrom^2'), symmetry=1, barrier=(4.77763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104276,'amu*angstrom^2'), symmetry=1, barrier=(72.0127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.5147e-06,4.19238e-10,1.90524e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.57,'K')), NASAPolynomial(coeffs=[12.4232,0.018224,-7.06311e-06,1.16768e-09,-7.11811e-14,25091.4,-39.6217], Tmin=(2154.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C[CH]C(305)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  C u1 p0 c0 {1,S} {2,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (220.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,1975.35,1978.55],'cm^-1')),
        HinderedRotor(inertia=(0.119612,'amu*angstrom^2'), symmetry=1, barrier=(2.75011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116659,'amu*angstrom^2'), symmetry=1, barrier=(2.68221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119045,'amu*angstrom^2'), symmetry=1, barrier=(2.73707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119229,'amu*angstrom^2'), symmetry=1, barrier=(2.74131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5835,0.0375164,-1.14014e-05,-2.10757e-10,4.03354e-13,26602.2,21.7084], Tmin=(100,'K'), Tmax=(2050.39,'K')), NASAPolynomial(coeffs=[12.6275,0.02724,-1.07003e-05,1.7777e-09,-1.09341e-13,20524.7,-38.7366], Tmin=(2050.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CC(308)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {5,S} {14,S}
5  C u1 p0 c0 {3,S} {4,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (220.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2185.76,2186.06],'cm^-1')),
        HinderedRotor(inertia=(0.0968109,'amu*angstrom^2'), symmetry=1, barrier=(2.22587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968729,'amu*angstrom^2'), symmetry=1, barrier=(2.2273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0967461,'amu*angstrom^2'), symmetry=1, barrier=(2.22438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0967567,'amu*angstrom^2'), symmetry=1, barrier=(2.22463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82915,0.0356449,-8.78269e-06,-1.57849e-09,6.29729e-13,26591.6,21.3393], Tmin=(100,'K'), Tmax=(2058.94,'K')), NASAPolynomial(coeffs=[12.9604,0.0273789,-1.10779e-05,1.8577e-09,-1.1449e-13,19999.8,-40.731], Tmin=(2058.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC[CH]C(286)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (231.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2520.05,2520.11],'cm^-1')),
        HinderedRotor(inertia=(0.00289438,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115564,'amu*angstrom^2'), symmetry=1, barrier=(4.795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32508,'amu*angstrom^2'), symmetry=1, barrier=(54.7434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11615,'amu*angstrom^2'), symmetry=1, barrier=(4.79424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3051.59,'J/mol'), sigma=(5.73385,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.65 K, Pc=36.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95927,0.0431639,-1.8953e-05,3.28952e-09,-1.25682e-13,27931,23.7177], Tmin=(100,'K'), Tmax=(1936.26,'K')), NASAPolynomial(coeffs=[12.6501,0.0264646,-1.01888e-05,1.70861e-09,-1.07058e-13,22781.3,-37.5306], Tmin=(1936.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC(307)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {5,S} {13,S}
5  C u1 p0 c0 {4,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (231.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2520.05,2520.11],'cm^-1')),
        HinderedRotor(inertia=(0.00289438,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115564,'amu*angstrom^2'), symmetry=1, barrier=(4.795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32508,'amu*angstrom^2'), symmetry=1, barrier=(54.7434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11615,'amu*angstrom^2'), symmetry=1, barrier=(4.79424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95927,0.0431639,-1.8953e-05,3.28952e-09,-1.25682e-13,27931,23.7177], Tmin=(100,'K'), Tmax=(1936.26,'K')), NASAPolynomial(coeffs=[12.6501,0.0264646,-1.01888e-05,1.70861e-09,-1.07058e-13,22781.3,-37.5306], Tmin=(1936.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]CC(400)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
"""),
    E0 = (26.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,180,2316.45],'cm^-1')),
        HinderedRotor(inertia=(0.136457,'amu*angstrom^2'), symmetry=1, barrier=(3.13742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137521,'amu*angstrom^2'), symmetry=1, barrier=(3.16187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135752,'amu*angstrom^2'), symmetry=1, barrier=(3.1212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135992,'amu*angstrom^2'), symmetry=1, barrier=(3.12672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1409,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9332,0.043137,-1.41038e-05,-6.73484e-10,8.04661e-13,3248.32,21.1643], Tmin=(100,'K'), Tmax=(1648.18,'K')), NASAPolynomial(coeffs=[9.63731,0.033395,-1.33879e-05,2.33363e-09,-1.51515e-13,-507.59,-23.5431], Tmin=(1648.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCCC(401)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {15,S} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
"""),
    E0 = (37.1621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,183.764,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00457302,'amu*angstrom^2'), symmetry=1, barrier=(5.34156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222916,'amu*angstrom^2'), symmetry=1, barrier=(5.34155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922882,'amu*angstrom^2'), symmetry=1, barrier=(22.1141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922884,'amu*angstrom^2'), symmetry=1, barrier=(22.1141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1409,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58666,0.0452266,-8.34973e-06,-1.27366e-08,5.91002e-12,4562.96,22.3029], Tmin=(100,'K'), Tmax=(1123.45,'K')), NASAPolynomial(coeffs=[9.06514,0.0338791,-1.35994e-05,2.48445e-09,-1.70983e-13,1918.38,-18.9383], Tmin=(1123.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.1621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
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
    E0 = (17.2746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (457.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (22.9503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (262.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (272.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (323.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-0.40783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (22.3111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (224.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (223.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (264.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (264.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (275.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (275.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (17.0775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-44.5177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]CCC(392)'],
    products = ['C2H5(31)', 'C=CC(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.24e+30,'s^-1'), n=-4.86, Ea=(38050,'cal/mol'), T0=(1,'K'), comment="""Matched reaction 2923 C5H11-2 <=> propene_1 + C2H5 in R_Addition_MultipleBond/training
This reaction matched rate rule [Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R_4R!H->C_Sp-2CNO=1CCNNO_5R!H->C_Sp-5C-2CNO]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', 'C[CH]CC(285)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(154.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)CC(391)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(38)', '[CH2]CC(48)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(18)', '[CH]CCC(398)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(6)', 'C[C]CCC(399)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(6)', 'CC=CCC(301)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6531.76,'m^3/(mol*s)'), n=1.11533, Ea=(8.6234,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05436079868881756, var=0.05901623039690775, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS_5R!H->C_Ext-4CCl-R_N-Sp-6R!H=4CCClCl',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS_5R!H->C_Ext-4CCl-R_N-Sp-6R!H=4CCClCl"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', 'C=CCCC(288)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2768.44,'m^3/(mol*s)'), n=1.26631, Ea=(19.097,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3314163891060069, var=0.0061045007265244575, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing_Ext-5C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing_Ext-5C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H5(31)', '[CH2][CH]C(195)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.41902e+08,'m^3/(mol*s)'), n=-0.65, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.007291370132393764, var=0.8826412212865885, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_5R!H->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(18)', '[CH2]C[CH]C(55)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(6)', 'C[CH]C[CH]C(305)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.77266e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -0.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(6)', 'C[CH][CH]CC(308)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R
Ea raised from -0.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(6)', '[CH2]CC[CH]C(286)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(6)', '[CH2][CH]CCC(307)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC[CH]CC(400)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.352e+10,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CCCC(401)'],
    products = ['C[CH]CCC(392)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.22e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 131 C5H11 <=> C5H11-2 in intra_H_migration/training
This reaction matched rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #249',
    isomers = [
        'C[CH]CCC(392)',
    ],
    reactants = [
        ('C2H5(31)', 'C=CC(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #249',
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

