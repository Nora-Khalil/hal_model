species(
    label = '[CH2]C(=CC)C[O](2406)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (155.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,1015.66],'cm^-1')),
        HinderedRotor(inertia=(0.161482,'amu*angstrom^2'), symmetry=1, barrier=(3.7128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160572,'amu*angstrom^2'), symmetry=1, barrier=(3.69186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15986,'amu*angstrom^2'), symmetry=1, barrier=(3.67551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52271,0.0395561,1.71444e-05,-1.00883e-07,8.86473e-11,18701.4,20.0938], Tmin=(100,'K'), Tmax=(431.382,'K')), NASAPolynomial(coeffs=[3.36713,0.0408207,-1.88763e-05,3.6557e-09,-2.59099e-13,18544,15.7491], Tmin=(431.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'CH2O(19)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[CH2]C(=C)C[O](1210)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {8,S} {9,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (191.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,649.084,649.211],'cm^-1')),
        HinderedRotor(inertia=(0.0123354,'amu*angstrom^2'), symmetry=1, barrier=(3.68568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160288,'amu*angstrom^2'), symmetry=1, barrier=(3.68533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3511,'J/mol'), sigma=(5.95716,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.41 K, Pc=37.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.65316,0.0333255,-1.65856e-05,3.65186e-09,-3.087e-13,23032.8,16.9441], Tmin=(100,'K'), Tmax=(2703.16,'K')), NASAPolynomial(coeffs=[16.1268,0.0133879,-5.52211e-06,9.23337e-10,-5.63548e-14,15748.5,-61.4571], Tmin=(2703.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC[O](1988)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (251.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,180,180,2817.6],'cm^-1')),
        HinderedRotor(inertia=(0.0317877,'amu*angstrom^2'), symmetry=1, barrier=(0.730863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246844,'amu*angstrom^2'), symmetry=1, barrier=(5.67543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583169,'amu*angstrom^2'), symmetry=1, barrier=(13.4082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3713.27,'J/mol'), sigma=(6.31777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.00 K, Pc=33.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19859,0.0404045,-1.85267e-05,3.3882e-09,-2.16201e-13,30213.5,18.6557], Tmin=(100,'K'), Tmax=(2925.15,'K')), NASAPolynomial(coeffs=[38.9674,-0.00172085,-4.05301e-07,5.13529e-11,1.19541e-15,6384.08,-197.264], Tmin=(2925.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'O(7)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C([CH2])=CC(1459)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,D} {4,S} {5,S}
3  C u0 p0 c0 {1,S} {2,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (234.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.0177687,'amu*angstrom^2'), symmetry=1, barrier=(20.2237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878272,'amu*angstrom^2'), symmetry=1, barrier=(20.1932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.60896,'amu*angstrom^2'), symmetry=1, barrier=(105.969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80243,0.0390484,-7.98166e-07,-2.45991e-08,1.16042e-11,28236,18.1779], Tmin=(100,'K'), Tmax=(1004.18,'K')), NASAPolynomial(coeffs=[10.9852,0.023482,-8.93209e-06,1.63809e-09,-1.15441e-13,25332.3,-31.4369], Tmin=(1004.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CH2(T)(17)',
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
    label = 'CC=[C]C[O](1224)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {11,S}
5  C u1 p0 c0 {3,S} {4,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (280.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,180,2871.6],'cm^-1')),
        HinderedRotor(inertia=(0.255967,'amu*angstrom^2'), symmetry=1, barrier=(5.88519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255144,'amu*angstrom^2'), symmetry=1, barrier=(5.86626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48949,0.0397787,-5.33843e-05,5.63856e-08,-2.3594e-11,33785.5,19.4801], Tmin=(100,'K'), Tmax=(821.28,'K')), NASAPolynomial(coeffs=[0.0206279,0.0358523,-1.70804e-05,3.26799e-09,-2.26303e-13,34728.9,34.1799], Tmin=(821.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC=C1COC1(2410)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {4,S} {5,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-26.4653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34067,0.0242171,4.16058e-05,-6.23416e-08,2.2932e-11,-3112.61,18.7392], Tmin=(100,'K'), Tmax=(1018.08,'K')), NASAPolynomial(coeffs=[8.54096,0.0303215,-1.22747e-05,2.33386e-09,-1.67579e-13,-5953.92,-19.0386], Tmin=(1018.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.4653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
)

species(
    label = 'CC=C(C)C=O(2530)',
    structure = adjacencyList("""1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {13,S}
6  C u0 p0 c0 {1,D} {4,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-160.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13648,0.0457261,-2.41403e-05,5.35936e-09,-4.39313e-13,-19266.3,18.2889], Tmin=(100,'K'), Tmax=(2413.03,'K')), NASAPolynomial(coeffs=[23.9351,0.0134535,-6.47978e-06,1.14344e-09,-7.12488e-14,-30910.9,-108.409], Tmin=(2413.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC1C[C]1C[O](2531)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {3,S} {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (194.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2181,0.0273726,3.43695e-05,-5.60965e-08,2.10662e-11,23513.2,20.9994], Tmin=(100,'K'), Tmax=(1016.51,'K')), NASAPolynomial(coeffs=[8.91378,0.0299085,-1.19947e-05,2.26426e-09,-1.61858e-13,20659.7,-18.7531], Tmin=(1016.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1COC1C(2532)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {3,S} {6,S}
6  C u1 p0 c0 {5,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (192.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2388,0.02392,6.00388e-05,-1.00943e-07,4.41979e-11,23180.9,19.641], Tmin=(100,'K'), Tmax=(872.561,'K')), NASAPolynomial(coeffs=[11.9858,0.0197684,-2.49989e-06,7.30424e-11,3.06416e-15,19937,-34.8954], Tmin=(872.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([CH]C)CO1(2494)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {4,S} {12,S}
6  C u1 p0 c0 {2,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (241.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261731,0.0670971,-6.51388e-05,3.39042e-08,-6.6421e-12,29185.8,24.672], Tmin=(100,'K'), Tmax=(1491.99,'K')), NASAPolynomial(coeffs=[13.9118,0.0175466,-2.29794e-06,5.28439e-12,1.3211e-14,26554.5,-41.8113], Tmin=(1491.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CJC(C)OC)"""),
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
    label = '[CH2]C(C=O)=CC(1791)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,D} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u0 p0 c0 {1,D} {3,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-4.35622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.765669,'amu*angstrom^2'), symmetry=1, barrier=(17.6042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128394,'amu*angstrom^2'), symmetry=1, barrier=(6.23172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767279,'amu*angstrom^2'), symmetry=1, barrier=(17.6412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69953,0.0484286,-3.10728e-05,9.3697e-09,-1.12942e-12,-439.904,19.3608], Tmin=(100,'K'), Tmax=(1872.11,'K')), NASAPolynomial(coeffs=[14.043,0.0220551,-9.94147e-06,1.84472e-09,-1.24545e-13,-5061.59,-47.9297], Tmin=(1872.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.35622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2][O](117)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88411,-0.00363912,3.28554e-05,-4.13625e-08,1.59638e-11,23210.8,7.47975], Tmin=(100,'K'), Tmax=(933.052,'K')), NASAPolynomial(coeffs=[6.69328,0.000290109,8.61345e-07,-1.56334e-10,7.33637e-15,21991.3,-9.60391], Tmin=(933.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2]C([CH]C)=CO(2533)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {14,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u1 p0 c0 {2,S} {3,S} {10,S}
5  C u0 p0 c0 {1,S} {3,D} {11,S}
6  C u1 p0 c0 {3,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (28.2424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980699,0.0481608,1.39047e-05,-6.56182e-08,3.25637e-11,3522.37,20.7054], Tmin=(100,'K'), Tmax=(925.67,'K')), NASAPolynomial(coeffs=[20.6066,0.0103927,-1.11853e-06,9.83959e-11,-1.08883e-14,-2126.35,-83.3478], Tmin=(925.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.2424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[CH]C(C)=C[O](2534)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u1 p0 c0 {3,S} {4,S} {13,S}
6  C u0 p0 c0 {1,S} {4,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (18.2058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38406,0.0437979,8.26643e-06,-4.52147e-08,2.13463e-11,2296.24,20.3866], Tmin=(100,'K'), Tmax=(962.096,'K')), NASAPolynomial(coeffs=[15.1357,0.0200059,-6.68448e-06,1.20877e-09,-8.78593e-14,-1894.78,-53.4546], Tmin=(962.096,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.2058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = 'C[C]=C(C)C[O](2535)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
"""),
    E0 = (241.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,1685,370,275.463,3122.02],'cm^-1')),
        HinderedRotor(inertia=(0.129592,'amu*angstrom^2'), symmetry=1, barrier=(6.96935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130226,'amu*angstrom^2'), symmetry=1, barrier=(6.9629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12973,'amu*angstrom^2'), symmetry=1, barrier=(6.96186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67405,0.0381042,-1.58769e-05,2.2425e-09,-5.97156e-14,29019.4,15.7991], Tmin=(100,'K'), Tmax=(2797.11,'K')), NASAPolynomial(coeffs=[41.0501,-0.00345969,3.85301e-08,-3.16487e-11,7.76294e-15,3460.94,-211.274], Tmin=(2797.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=[C]C)CO(2536)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {14,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (167.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.169112,'amu*angstrom^2'), symmetry=1, barrier=(3.88821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0671759,'amu*angstrom^2'), symmetry=1, barrier=(1.54451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169117,'amu*angstrom^2'), symmetry=1, barrier=(3.88834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.705405,'amu*angstrom^2'), symmetry=1, barrier=(16.2186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54318,0.0567869,-5.07587e-05,2.73209e-08,-6.42604e-12,20201.5,23.0649], Tmin=(100,'K'), Tmax=(987.374,'K')), NASAPolynomial(coeffs=[7.47361,0.0327614,-1.42592e-05,2.6763e-09,-1.86003e-13,19030.4,-5.47052], Tmin=(987.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(C)C[O](1964)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u1 p0 c0 {5,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (155.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,1015.66],'cm^-1')),
        HinderedRotor(inertia=(0.161482,'amu*angstrom^2'), symmetry=1, barrier=(3.7128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160572,'amu*angstrom^2'), symmetry=1, barrier=(3.69186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15986,'amu*angstrom^2'), symmetry=1, barrier=(3.67551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52271,0.0395561,1.71444e-05,-1.00883e-07,8.86473e-11,18701.4,20.0938], Tmin=(100,'K'), Tmax=(431.382,'K')), NASAPolynomial(coeffs=[3.36713,0.0408207,-1.88763e-05,3.6557e-09,-2.59099e-13,18544,15.7491], Tmin=(431.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C([CH2])CO(1968)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {6,S} {9,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (80.9057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18601,0.0543738,-3.40959e-05,6.60754e-09,1.12999e-12,9838.4,24.1464], Tmin=(100,'K'), Tmax=(1079.44,'K')), NASAPolynomial(coeffs=[12.4208,0.0247191,-9.53141e-06,1.71575e-09,-1.17793e-13,6715.16,-34.146], Tmin=(1079.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(C)C[O](1967)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {6,D} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (255.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,257.063,257.077,1691.27],'cm^-1')),
        HinderedRotor(inertia=(0.00255168,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167596,'amu*angstrom^2'), symmetry=1, barrier=(7.85519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167559,'amu*angstrom^2'), symmetry=1, barrier=(7.85507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35794,0.042816,2.08161e-07,-6.41224e-08,5.82513e-11,30775.3,22.1369], Tmin=(100,'K'), Tmax=(465.458,'K')), NASAPolynomial(coeffs=[3.77534,0.0395797,-1.81867e-05,3.50883e-09,-2.48012e-13,30546.5,15.3419], Tmin=(465.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S)"""),
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
    label = 'C=[C]C[O](147)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {4,D} {7,S} {8,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (319.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.110137,'amu*angstrom^2'), symmetry=1, barrier=(2.53227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0632,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63077,0.0408503,-8.2829e-05,9.16991e-08,-3.63161e-11,38513.3,13.7006], Tmin=(100,'K'), Tmax=(865.187,'K')), NASAPolynomial(coeffs=[-0.468585,0.0281972,-1.41119e-05,2.70319e-09,-1.84571e-13,40059.5,34.0405], Tmin=(865.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C1COC1C(2409)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-32.7561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22975,0.0200924,6.89951e-05,-1.01222e-07,3.90834e-11,-3859.01,19.0612], Tmin=(100,'K'), Tmax=(973.819,'K')), NASAPolynomial(coeffs=[13.022,0.0230528,-8.40674e-06,1.63287e-09,-1.23592e-13,-8203.24,-44.2315], Tmin=(973.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.7561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C(C=O)CC(2537)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {1,D} {4,S} {12,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-147.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56344,0.0483204,-2.62962e-05,5.63633e-09,-2.89268e-13,-17627.3,21.484], Tmin=(100,'K'), Tmax=(1654.89,'K')), NASAPolynomial(coeffs=[14.6772,0.0238325,-1.06344e-05,1.95931e-09,-1.31443e-13,-22954.8,-51.3704], Tmin=(1654.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=C)CO(1959)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-96.8052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25696,0.0518781,-2.5319e-05,-5.59105e-09,6.40816e-12,-11536.7,22.1956], Tmin=(100,'K'), Tmax=(988.696,'K')), NASAPolynomial(coeffs=[13.1466,0.0220885,-7.90686e-06,1.40207e-09,-9.7124e-14,-14782.8,-39.5565], Tmin=(988.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.8052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH][C]1COC1(2538)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {3,S} {6,S}
6  C u1 p0 c0 {4,S} {5,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (194.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4045,0.0233936,4.90049e-05,-7.80431e-08,3.22485e-11,23424.6,20.1611], Tmin=(100,'K'), Tmax=(902.362,'K')), NASAPolynomial(coeffs=[8.52744,0.0267855,-7.38969e-06,1.12023e-09,-7.28848e-14,21076.5,-15.6376], Tmin=(902.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJ(C)CO) + radical(Cs_S)"""),
)

species(
    label = 'C=CC(=C)C[O](1963)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (128.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,262.961,268.179,1226.85],'cm^-1')),
        HinderedRotor(inertia=(0.198828,'amu*angstrom^2'), symmetry=1, barrier=(10.032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300767,'amu*angstrom^2'), symmetry=1, barrier=(14.973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75049,0.0475522,-3.17511e-05,1.10669e-08,-1.61624e-12,15585.3,21.2845], Tmin=(100,'K'), Tmax=(1535.91,'K')), NASAPolynomial(coeffs=[10.0149,0.0260291,-1.07314e-05,1.94327e-09,-1.31205e-13,13046.6,-22.1332], Tmin=(1535.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC(=C)C[O](2539)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (222.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1434.8,1434.81],'cm^-1')),
        HinderedRotor(inertia=(0.145086,'amu*angstrom^2'), symmetry=1, barrier=(3.33582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145082,'amu*angstrom^2'), symmetry=1, barrier=(3.33572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145054,'amu*angstrom^2'), symmetry=1, barrier=(3.33509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57614,0.0428805,-2.13004e-05,4.65011e-09,-3.92401e-13,26767.4,22.2694], Tmin=(100,'K'), Tmax=(2532.56,'K')), NASAPolynomial(coeffs=[13.6141,0.0254465,-1.09742e-05,1.93181e-09,-1.24061e-13,21176.6,-41.239], Tmin=(2532.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(=C[O])CC(2540)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {4,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (28.5924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19253,0.0464393,8.23164e-06,-5.06781e-08,2.47007e-11,3553.9,22.4613], Tmin=(100,'K'), Tmax=(947.574,'K')), NASAPolynomial(coeffs=[17.2526,0.0166628,-4.81409e-06,8.41188e-10,-6.26361e-14,-1196.51,-63.1604], Tmin=(947.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.5924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(CC)C[O](2541)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u1 p0 c0 {5,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (264.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,245.268,1969.32],'cm^-1')),
        HinderedRotor(inertia=(0.136421,'amu*angstrom^2'), symmetry=1, barrier=(5.86776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136001,'amu*angstrom^2'), symmetry=1, barrier=(5.85772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136809,'amu*angstrom^2'), symmetry=1, barrier=(5.87852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47681,0.0411637,1.53558e-05,-1.09166e-07,1.00486e-10,31809.3,21.505], Tmin=(100,'K'), Tmax=(430.035,'K')), NASAPolynomial(coeffs=[3.85089,0.0395888,-1.82392e-05,3.51164e-09,-2.47391e-13,31587.5,14.8307], Tmin=(430.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=CC)CO(2542)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {12,S}
6  C u2 p0 c0 {4,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (148.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60769,0.054445,-3.43396e-05,1.1813e-08,-1.80699e-12,17955.6,23.1713], Tmin=(100,'K'), Tmax=(1385.57,'K')), NASAPolynomial(coeffs=[7.54431,0.0373065,-1.57856e-05,2.88571e-09,-1.96227e-13,16310.5,-7.40542], Tmin=(1385.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-21.4268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (433.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (244.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (300.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (485.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-13.1425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (41.9733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (209.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (104.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (66.1129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (79.3886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (48.8626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (167.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (377.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (90.4988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (181.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (226.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (35.0184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (128.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (67.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (173.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (487.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-13.1425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (41.9733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (3.54644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (104.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (175.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (159.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (187.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (232.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (397.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['CH2O(19)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH2]C(=C)C[O](1210)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.07568e+09,'m^3/(mol*s)'), n=-0.779019, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_4CbCdCt->Cd',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_4CbCdCt->Cd
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC=[C]CC[O](1988)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(7)', '[CH2]C([CH2])=CC(1459)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cd;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(17)', 'CC=[C]C[O](1224)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['CC=C1COC1(2410)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['CC=C(C)C=O(2530)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['CC1C[C]1C[O](2531)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH2][C]1COC1C(2532)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH2]C1([CH]C)CO1(2494)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(87.5397,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2O(19)', '[CH2][C]=CC(1693)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(13.9251,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(6)', '[CH2]C(C=O)=CC(1791)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.5,'m^3/(mol*s)'), n=2.16, Ea=(17.9525,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_Sp-5R!H=4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_Sp-5R!H=4C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][O](117)', 'C=C=CC(1692)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(5.68222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][O](117)', '[CH2][C]=CC(1693)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH2]C([CH]C)=CO(2533)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['C[CH]C(C)=C[O](2534)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[C]=C(C)C[O](2535)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=[C]C)CO(2536)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=C(C)C[O](1964)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH2]C=C([CH2])CO(1968)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(C)C[O](1967)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(38)', 'C=[C]C[O](147)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['C=C1COC1C(2409)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['C=C(C=O)CC(2537)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['C=CC(=C)CO(1959)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['C[CH][C]1COC1(2538)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(6)', 'C=CC(=C)C[O](1963)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2768.44,'m^3/(mol*s)'), n=1.26631, Ea=(11.4193,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3314163891060069, var=0.0061045007265244575, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing_Ext-5C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing_Ext-5C-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CC(=C)C[O](2539)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH2]C(=C[O])CC(2540)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.23617e+10,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(CC)C[O](2541)'],
    products = ['[CH2]C(=CC)C[O](2406)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=CC)C[O](2406)'],
    products = ['[CH]C(=CC)CO(2542)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.4207e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1197',
    isomers = [
        '[CH2]C(=CC)C[O](2406)',
    ],
    reactants = [
        ('CH2O(19)', 'C=C=CC(1692)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1197',
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

