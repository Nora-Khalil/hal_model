species(
    label = '[CH2]OC1[CH]CC1(1932)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {4,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (245.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82231,0.0317926,3.78656e-05,-7.18416e-08,2.95738e-11,29636,22.5972], Tmin=(100,'K'), Tmax=(971.461,'K')), NASAPolynomial(coeffs=[13.9366,0.021594,-7.65917e-06,1.44792e-09,-1.07716e-13,25409.8,-45.134], Tmin=(971.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CsJOCC2)"""),
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
    label = 'C1=CCC1(69)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (144.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,204.182,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,3484.69],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3170.39,'J/mol'), sigma=(5.50295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.21 K, Pc=43.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39615,-0.00343641,9.09434e-05,-1.1327e-07,4.24092e-11,17412.4,10.7746], Tmin=(100,'K'), Tmax=(950.248,'K')), NASAPolynomial(coeffs=[9.42764,0.013964,-4.06869e-06,7.74979e-10,-6.20376e-14,14334.2,-28.1812], Tmin=(950.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
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
    label = '[CH2]OC1[CH]C1(2223)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u1 p0 c0 {1,S} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (273.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25766,0.0251167,2.86897e-05,-5.84452e-08,2.50371e-11,33009.5,18.7758], Tmin=(100,'K'), Tmax=(952.077,'K')), NASAPolynomial(coeffs=[13.0725,0.0126919,-3.74556e-06,6.85772e-10,-5.28182e-14,29454,-40.7262], Tmin=(952.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CCJCO) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]O[CH]C1CC1(2224)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {2,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (238.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47017,0.0422965,9.23195e-06,-4.64038e-08,2.21801e-11,28775.7,23.1333], Tmin=(100,'K'), Tmax=(947.229,'K')), NASAPolynomial(coeffs=[14.945,0.01853,-5.60462e-06,9.68769e-10,-6.98373e-14,24736.5,-48.9914], Tmin=(947.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CCsJOCs) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]OC1CC1[CH2](1934)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (246.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41554,0.0446547,2.47085e-06,-3.90686e-08,1.95986e-11,29714.4,22.3552], Tmin=(100,'K'), Tmax=(943.033,'K')), NASAPolynomial(coeffs=[14.4732,0.0196113,-5.95765e-06,1.00874e-09,-7.1042e-14,25902.5,-47.0278], Tmin=(943.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CsJOCC2)"""),
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
    label = '[O]C1[CH]CC1(2225)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (271.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88602,0.00621909,8.18893e-05,-1.10736e-07,4.26786e-11,32730.5,17.8325], Tmin=(100,'K'), Tmax=(952.307,'K')), NASAPolynomial(coeffs=[12.0159,0.014373,-4.20099e-06,8.09161e-10,-6.53679e-14,28883,-36.8388], Tmin=(952.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C1CC2OCC12(2199)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (18.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48178,0.0172553,6.72741e-05,-9.06983e-08,3.29994e-11,2288.83,14.4361], Tmin=(100,'K'), Tmax=(1006.22,'K')), NASAPolynomial(coeffs=[9.67852,0.0295962,-1.21681e-05,2.38088e-09,-1.75333e-13,-1232.52,-30.6299], Tmin=(1006.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'COC1=CCC1(2226)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-51.4069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84886,0.0273792,5.74305e-05,-9.82653e-08,4.06326e-11,-6087.22,18.7134], Tmin=(100,'K'), Tmax=(950.246,'K')), NASAPolynomial(coeffs=[16.3504,0.0165199,-4.64517e-06,8.6219e-10,-6.84922e-14,-11109,-62.4301], Tmin=(950.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.4069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'COC1C=CC1(2200)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,D} {14,S}
6  C u0 p0 c0 {3,S} {5,D} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-7.59328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11365,0.026006,4.64079e-05,-7.50177e-08,2.92044e-11,-831.588,19.9833], Tmin=(100,'K'), Tmax=(988.502,'K')), NASAPolynomial(coeffs=[11.9338,0.0245415,-9.44679e-06,1.82026e-09,-1.34487e-13,-4642.93,-36.7378], Tmin=(988.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.59328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]OC([CH2])C=C(1847)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (225.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,212.491,212.491,212.493],'cm^-1')),
        HinderedRotor(inertia=(0.00373347,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00373361,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518955,'amu*angstrom^2'), symmetry=1, barrier=(16.629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519001,'amu*angstrom^2'), symmetry=1, barrier=(16.6291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3230.22,'J/mol'), sigma=(5.77019,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.55 K, Pc=38.15 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885469,0.06699,-6.57956e-05,3.41639e-08,-7.15716e-12,27236.3,25.0328], Tmin=(100,'K'), Tmax=(1148.23,'K')), NASAPolynomial(coeffs=[13.0988,0.0244437,-1.02153e-05,1.89403e-09,-1.31207e-13,24431.5,-35.5776], Tmin=(1148.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]CC=CO[CH2](2227)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {10,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (190.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,464.038,464.211,464.246],'cm^-1')),
        HinderedRotor(inertia=(0.000782475,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108521,'amu*angstrom^2'), symmetry=1, barrier=(16.5871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108591,'amu*angstrom^2'), symmetry=1, barrier=(16.5873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108568,'amu*angstrom^2'), symmetry=1, barrier=(16.5867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21576,0.0468261,2.99309e-06,-4.43401e-08,2.2363e-11,23040.1,25.8369], Tmin=(100,'K'), Tmax=(947.301,'K')), NASAPolynomial(coeffs=[17.0208,0.0158335,-4.53135e-06,7.87444e-10,-5.84993e-14,18441.9,-58.023], Tmin=(947.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ) + radical(C=COCJ)"""),
)

species(
    label = '[CH]1[CH]CC1(2183)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {4,S} {9,S}
4  C u1 p0 c0 {2,S} {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (389.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,180,388.393,1395.57,1395.58,1395.59,1395.6,1395.6,1395.6,1395.6,1395.61,1395.61,1802.73],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61382,-0.00133354,6.67738e-05,-7.44503e-08,2.51141e-11,46874.5,15.2545], Tmin=(100,'K'), Tmax=(1012.31,'K')), NASAPolynomial(coeffs=[3.96904,0.0238095,-9.81749e-06,1.89418e-09,-1.37336e-13,45442.4,6.81812], Tmin=(1012.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
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
    label = '[CH2]OC1=CCC1(2228)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  C u1 p0 c0 {1,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (140.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96597,0.0252722,5.579e-05,-9.56066e-08,3.96976e-11,17001.3,20.7785], Tmin=(100,'K'), Tmax=(948.557,'K')), NASAPolynomial(coeffs=[16.345,0.0137237,-3.57126e-06,6.69261e-10,-5.52709e-14,12065.1,-59.4732], Tmin=(948.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC1C=CC1(2229)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,D} {10,S}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  C u1 p0 c0 {1,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (175.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80039,0.0326925,3.00053e-05,-6.46142e-08,2.74305e-11,21260.2,20.524], Tmin=(100,'K'), Tmax=(969.866,'K')), NASAPolynomial(coeffs=[14.8351,0.0174791,-6.08008e-06,1.16791e-09,-8.8784e-14,16918.9,-51.3081], Tmin=(969.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]O[C]1CCC1(2230)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (226.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45802,0.0484298,-2.33275e-05,5.39077e-10,1.9928e-12,27331.8,21.3617], Tmin=(100,'K'), Tmax=(1168.36,'K')), NASAPolynomial(coeffs=[10.9933,0.0273482,-1.11073e-05,2.03707e-09,-1.40273e-13,24314.4,-29.5019], Tmin=(1168.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(C2CsJOCs) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2]OC1C[CH]C1(2231)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {4,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (233.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78358,0.03766,1.09462e-05,-3.62768e-08,1.50388e-11,28181.6,22.5894], Tmin=(100,'K'), Tmax=(1028.83,'K')), NASAPolynomial(coeffs=[11.0592,0.0271521,-1.09924e-05,2.0824e-09,-1.4895e-13,24920.5,-28.997], Tmin=(1028.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(CsJOCC2)"""),
)

species(
    label = 'CO[C]1[CH]CC1(2232)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {2,S} {6,S}
6  C u1 p0 c0 {3,S} {5,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (242.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00557,0.0360925,3.0697e-06,-2.08788e-08,8.14742e-12,29274.9,23.1679], Tmin=(100,'K'), Tmax=(1113.58,'K')), NASAPolynomial(coeffs=[8.22902,0.0308073,-1.28042e-05,2.38979e-09,-1.66736e-13,26830.4,-12.2783], Tmin=(1113.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(C2CsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'COC1[CH]C[CH]1(2233)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {3,S} {13,S}
6  C u1 p0 c0 {2,S} {3,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (261.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33002,0.019904,6.2795e-05,-9.15063e-08,3.50451e-11,31580.8,23.8545], Tmin=(100,'K'), Tmax=(976.868,'K')), NASAPolynomial(coeffs=[11.5038,0.0245196,-9.06048e-06,1.73284e-09,-1.2868e-13,27775.9,-30.4903], Tmin=(976.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = 'COC1[CH][CH]C1(2234)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {6,S} {14,S}
6  C u1 p0 c0 {3,S} {5,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (249.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29884,0.0256766,3.62359e-05,-5.6452e-08,2.0748e-11,30126.1,23.8199], Tmin=(100,'K'), Tmax=(1029.61,'K')), NASAPolynomial(coeffs=[8.61962,0.0300913,-1.24022e-05,2.36945e-09,-1.70099e-13,27288.9,-14.3158], Tmin=(1029.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(cyclobutane)"""),
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
    E0 = (30.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (631.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (190.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (191.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (438.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (39.0304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (94.1462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (55.7193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (134.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (98.7439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (133.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (174.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (151.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (178.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (367.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (185.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (189.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (132.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (130.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (63.6575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['CH2O(19)', 'C1=CCC1(69)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH2]OC1[CH]C1(2223)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(152.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['[CH2]O[CH]C1CC1(2224)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]OC1CC1[CH2](1934)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(17)', '[O]C1[CH]CC1(2225)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['C1CC2OCC12(2199)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['COC1=CCC1(2226)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['COC1C=CC1(2200)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC([CH2])C=C(1847)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC=CO[CH2](2227)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2O(19)', '[CH]1[CH]CC1(2183)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(77.6593,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][O](117)', 'C1=CCC1(69)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.60317e+11,'m^3/(mol*s)'), n=-1.23787, Ea=(52.3208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08456738906112354, var=37.926887421795335, Tref=1000.0, N=24, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(6)', '[CH2]OC1=CCC1(2228)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(14.4581,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(6)', '[CH2]OC1C=CC1(2229)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1020.43,'m^3/(mol*s)'), n=1.51234, Ea=(5.13068,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08797594653728344, var=6.974501284715509, Tref=1000.0, N=21, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][O](117)', '[CH]1[CH]CC1(2183)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['[CH2]O[C]1CCC1(2230)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['[CH2]OC1C[CH]C1(2231)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['CO[C]1[CH]CC1(2232)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(50.7042,'s^-1'), n=3.11103, Ea=(102.21,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;XH_out] for rate rule [R3H_SS_O;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['COC1[CH]C[CH]1(2233)'],
    products = ['[CH2]OC1[CH]CC1(1932)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.43333e+06,'s^-1'), n=1.5, Ea=(83.8474,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 436 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC1[CH]CC1(1932)'],
    products = ['COC1[CH][CH]C1(2234)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.13447e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS_OCC;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1130',
    isomers = [
        '[CH2]OC1[CH]CC1(1932)',
    ],
    reactants = [
        ('CH2O(19)', 'C1=CCC1(69)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1130',
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

