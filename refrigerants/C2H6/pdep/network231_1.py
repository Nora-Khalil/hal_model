species(
    label = '[CH2]OC[CH]C(317)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {3,S} {11,S}
5  C u1 p0 c0 {1,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (128.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,990.853,992.354,2395.41],'cm^-1')),
        HinderedRotor(inertia=(0.126316,'amu*angstrom^2'), symmetry=1, barrier=(2.90425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125495,'amu*angstrom^2'), symmetry=1, barrier=(2.88538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123534,'amu*angstrom^2'), symmetry=1, barrier=(2.8403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865235,'amu*angstrom^2'), symmetry=1, barrier=(19.8935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93357,0.0394074,-1.46731e-05,-4.27958e-09,3.33597e-12,15482.5,23.0053], Tmin=(100,'K'), Tmax=(1091.7,'K')), NASAPolynomial(coeffs=[8.82256,0.0249718,-9.6856e-06,1.74124e-09,-1.19024e-13,13334.5,-13.7837], Tmin=(1091.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CCJCO) + radical(CsJOCC)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[CH]C[O](125)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {3,S} {10,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (152.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1354.33,2181.55],'cm^-1')),
        HinderedRotor(inertia=(0.00394475,'amu*angstrom^2'), symmetry=1, barrier=(5.13456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222159,'amu*angstrom^2'), symmetry=1, barrier=(5.10788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10998,0.0164269,1.17436e-05,-1.70306e-08,5.13362e-12,18403.8,17.5215], Tmin=(100,'K'), Tmax=(1242.44,'K')), NASAPolynomial(coeffs=[4.47918,0.0232159,-9.97099e-06,1.87457e-09,-1.29965e-13,17199.3,7.14077], Tmin=(1242.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJCO)"""),
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
    label = '[CH2]O[CH2](115)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (177.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,393.523],'cm^-1')),
        HinderedRotor(inertia=(0.00108835,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00108835,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05159,0.0165783,-1.0043e-07,-1.04838e-08,4.7883e-12,21436.5,11.5123], Tmin=(100,'K'), Tmax=(1037.89,'K')), NASAPolynomial(coeffs=[7.39314,0.0094396,-3.64813e-06,7.0084e-10,-5.09573e-14,20018.6,-12.0839], Tmin=(1037.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOCH3) + radical(CsJOCH3)"""),
)

species(
    label = 'CC1COC1(319)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-129.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66098,0.0155577,6.65032e-05,-9.99729e-08,4.18803e-11,-15470.7,13.5541], Tmin=(100,'K'), Tmax=(882.039,'K')), NASAPolynomial(coeffs=[9.94109,0.0195445,-3.20213e-06,2.72548e-10,-1.29891e-14,-18194.3,-28.8136], Tmin=(882.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = 'CC=COC(335)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {1,S} {4,D} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-183.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72233,0.0361834,1.80278e-05,-5.43587e-08,2.50042e-11,-22033.5,16.8344], Tmin=(100,'K'), Tmax=(942.222,'K')), NASAPolynomial(coeffs=[14.904,0.0151733,-4.16405e-06,7.10665e-10,-5.28078e-14,-26068.9,-54.2076], Tmin=(942.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'C=CCOC(321)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-126.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18124,0.0373033,-1.57427e-05,1.80328e-09,2.07599e-13,-15118.6,18.8097], Tmin=(100,'K'), Tmax=(1689.22,'K')), NASAPolynomial(coeffs=[10.2866,0.0244594,-9.97542e-06,1.75219e-09,-1.14141e-13,-18762.8,-27.2244], Tmin=(1689.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]OC=CC(336)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {10,S}
5  C u1 p0 c0 {1,S} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (8.02273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,365.767,365.818],'cm^-1')),
        HinderedRotor(inertia=(0.171584,'amu*angstrom^2'), symmetry=1, barrier=(16.2943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171593,'amu*angstrom^2'), symmetry=1, barrier=(16.2939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171571,'amu*angstrom^2'), symmetry=1, barrier=(16.2941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83964,0.0340734,1.64014e-05,-5.17247e-08,2.40829e-11,1055.01,18.8989], Tmin=(100,'K'), Tmax=(939.523,'K')), NASAPolynomial(coeffs=[14.9004,0.0123741,-3.08846e-06,5.17351e-10,-3.9555e-14,-2895.65,-51.2611], Tmin=(939.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.02273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OCC=C(337)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (61.5053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,883.881,884.847],'cm^-1')),
        HinderedRotor(inertia=(0.15537,'amu*angstrom^2'), symmetry=1, barrier=(3.57227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154913,'amu*angstrom^2'), symmetry=1, barrier=(3.56176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155279,'amu*angstrom^2'), symmetry=1, barrier=(3.57018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14796,0.0354428,-1.34464e-05,-2.34155e-09,2.12167e-12,7468.41,20.5996], Tmin=(100,'K'), Tmax=(1162.85,'K')), NASAPolynomial(coeffs=[8.08549,0.0236411,-9.34537e-06,1.68379e-09,-1.14597e-13,5504.54,-11.448], Tmin=(1162.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.5053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]O[CH]CC(329)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {11,S}
5  C u1 p0 c0 {1,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (108.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4376,0.0516578,-4.26023e-05,1.91668e-08,-3.52039e-12,13160.4,21.7333], Tmin=(100,'K'), Tmax=(1297,'K')), NASAPolynomial(coeffs=[11.1425,0.0217273,-7.98713e-06,1.37431e-09,-9.08314e-14,10642.9,-27.6111], Tmin=(1297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CCsJOCs) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]CCO[CH2](338)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (133.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,983.459,983.475],'cm^-1')),
        HinderedRotor(inertia=(0.180146,'amu*angstrom^2'), symmetry=1, barrier=(4.14191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180134,'amu*angstrom^2'), symmetry=1, barrier=(4.14163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18016,'amu*angstrom^2'), symmetry=1, barrier=(4.14224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00603572,'amu*angstrom^2'), symmetry=1, barrier=(4.14233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93388,0.0435308,-2.74169e-05,9.13206e-09,-1.28728e-12,16121.4,22.2905], Tmin=(100,'K'), Tmax=(1569.37,'K')), NASAPolynomial(coeffs=[9.01099,0.0254928,-1.01763e-05,1.80831e-09,-1.20613e-13,13900.1,-15.042], Tmin=(1569.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(RCCJ) + radical(CsJOCC)"""),
)

species(
    label = 'C[CH][CH]OC(339)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {5,S} {12,S}
5  C u1 p0 c0 {1,S} {4,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (120.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73838,0.0477373,-3.55488e-05,1.44895e-08,-2.48309e-12,14606.1,21.4589], Tmin=(100,'K'), Tmax=(1344.29,'K')), NASAPolynomial(coeffs=[9.45208,0.0247847,-9.93759e-06,1.78822e-09,-1.21018e-13,12532.2,-18.0375], Tmin=(1344.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2][CH]COC(340)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {5,S} {11,S}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (145.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,343.268,1245.28,1245.28],'cm^-1')),
        HinderedRotor(inertia=(0.0054883,'amu*angstrom^2'), symmetry=1, barrier=(6.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722286,'amu*angstrom^2'), symmetry=1, barrier=(6.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722287,'amu*angstrom^2'), symmetry=1, barrier=(6.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722284,'amu*angstrom^2'), symmetry=1, barrier=(6.03949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22025,0.0398109,-2.10807e-05,5.27398e-09,-5.29631e-13,17567.4,22.0635], Tmin=(100,'K'), Tmax=(2151.17,'K')), NASAPolynomial(coeffs=[11.3373,0.0228585,-9.26003e-06,1.61073e-09,-1.0391e-13,13644.9,-28.9049], Tmin=(2151.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CCJCO) + radical(RCCJ)"""),
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
    E0 = (-0.503764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (405.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (393.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (7.78056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (62.8964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (24.4695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (76.7492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (130.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (103.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (147.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (349.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (108.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (156.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (101.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (71.1539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['CH2O(19)', 'C=CC(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(17)', 'C[CH]C[O](125)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(38)', '[CH2]O[CH2](115)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['CC1COC1(319)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['CC=COC(335)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['C=CCOC(321)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2O(19)', '[CH2][CH]C(195)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(39.5106,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][O](117)', 'C=CC(60)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.01583e+10,'m^3/(mol*s)'), n=-1.23787, Ea=(59.7063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08456738906112354, var=37.926887421795335, Tref=1000.0, N=24, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(6)', '[CH2]OC=CC(336)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(12.2345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[CH2]OCC=C(337)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(501.4,'m^3/(mol*s)'), n=1.733, Ea=(3.10828,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_N-5R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_N-5R!H->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][O](117)', '[CH2][CH]C(195)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['[CH2]O[CH]CC(329)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CCO[CH2](338)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OC[CH]C(317)'],
    products = ['C[CH][CH]OC(339)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(101.408,'s^-1'), n=3.11103, Ea=(102.21,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;XH_out] for rate rule [R3H_SS_O;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]COC(340)'],
    products = ['[CH2]OC[CH]C(317)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #231',
    isomers = [
        '[CH2]OC[CH]C(317)',
    ],
    reactants = [
        ('CH2O(19)', 'C=CC(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #231',
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

