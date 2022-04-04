species(
    label = '[CH2]CC=[C]C(1118)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {3,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (390.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,216.692],'cm^-1')),
        HinderedRotor(inertia=(0.14184,'amu*angstrom^2'), symmetry=1, barrier=(4.67359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13653,'amu*angstrom^2'), symmetry=1, barrier=(4.6726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132688,'amu*angstrom^2'), symmetry=1, barrier=(4.68272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04972,0.0405227,-2.17348e-05,5.57612e-09,-5.73512e-13,47043.1,21.1549], Tmin=(100,'K'), Tmax=(2150.77,'K')), NASAPolynomial(coeffs=[12.3921,0.021288,-8.32004e-06,1.41801e-09,-9.01843e-14,42594.3,-36.6614], Tmin=(2150.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C2H4(29)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C#CC(211)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (172.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.639],'cm^-1')),
        HinderedRotor(inertia=(0.630823,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30524,0.0109264,1.31992e-05,-2.25168e-08,8.87439e-12,20738.2,7.1916], Tmin=(100,'K'), Tmax=(969.981,'K')), NASAPolynomial(coeffs=[5.80052,0.0113536,-4.03494e-06,7.1907e-10,-5.02191e-14,19750,-7.36949], Tmin=(969.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C=[C]C(1189)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {4,D} {8,S}
3  C u1 p0 c0 {2,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,D}
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
    label = 'CC1=CCC1(1185)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (105.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55862,0.0164924,6.06482e-05,-8.6523e-08,3.31332e-11,12743.8,15.2307], Tmin=(100,'K'), Tmax=(968.638,'K')), NASAPolynomial(coeffs=[10.359,0.0225466,-7.98447e-06,1.49742e-09,-1.10505e-13,9437.45,-31.4195], Tmin=(968.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=CC=CC(382)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (57.8956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00727,0.0328459,1.55855e-05,-4.25745e-08,1.84259e-11,7044.82,16.9534], Tmin=(100,'K'), Tmax=(972.319,'K')), NASAPolynomial(coeffs=[11.2869,0.0212416,-7.50361e-06,1.3618e-09,-9.72233e-14,3984.25,-34.0139], Tmin=(972.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.8956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=CCC(1120)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {5,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (124.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00142,0.0359534,-1.58974e-07,-2.04778e-08,9.05633e-12,14999.4,18.1131], Tmin=(100,'K'), Tmax=(1043.78,'K')), NASAPolynomial(coeffs=[9.31671,0.0253743,-1.004e-05,1.85451e-09,-1.299e-13,12521.5,-22.0471], Tmin=(1043.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=[C]C(212)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (490.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859242,'amu*angstrom^2'), symmetry=1, barrier=(1.97557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04072e-06,1.20991e-09,-3.01638e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.72,'K')), NASAPolynomial(coeffs=[6.96329,0.0100049,-3.70699e-06,6.32786e-10,-4.05637e-14,57370,-10.4656], Tmin=(1807.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
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
    label = 'C=CC=[C]C(1122)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {4,D} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {3,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (295.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.779761,'amu*angstrom^2'), symmetry=1, barrier=(17.9282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779364,'amu*angstrom^2'), symmetry=1, barrier=(17.9191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96302,0.0375743,-9.62994e-06,-1.27883e-08,7.2586e-12,35648.7,17.5862], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[10.1968,0.0205755,-7.68659e-06,1.38864e-09,-9.67247e-14,33193.3,-26.1507], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC=C=C(1190)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (329.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,597.041],'cm^-1')),
        HinderedRotor(inertia=(0.411052,'amu*angstrom^2'), symmetry=1, barrier=(9.45089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818847,'amu*angstrom^2'), symmetry=1, barrier=(18.8269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02941,0.0367789,-1.08449e-05,-8.05942e-09,4.61387e-12,39682.5,19.8002], Tmin=(100,'K'), Tmax=(1084.8,'K')), NASAPolynomial(coeffs=[9.24557,0.0228621,-9.15063e-06,1.68413e-09,-1.17037e-13,37370.1,-19.0429], Tmin=(1084.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH2](39)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,971.411,1436.01,2688.07,2689.99],'cm^-1')),
        HinderedRotor(inertia=(0.0117579,'amu*angstrom^2'), symmetry=1, barrier=(17.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100313,6.40915e-06,-1.41289e-08,5.92664e-12,38288.2,6.11704], Tmin=(100,'K'), Tmax=(954.263,'K')), NASAPolynomial(coeffs=[5.5225,0.00856171,-2.90742e-06,5.02351e-10,-3.4457e-14,37547.7,-5.75282], Tmin=(954.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CC#CC(1191)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {5,T}
5  C u0 p0 c0 {2,S} {4,T}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (312.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2100,2250,500,550,2382.99],'cm^-1')),
        HinderedRotor(inertia=(0.131705,'amu*angstrom^2'), symmetry=1, barrier=(8.89945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08036,'amu*angstrom^2'), symmetry=1, barrier=(72.8278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07449,'amu*angstrom^2'), symmetry=1, barrier=(72.8341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33337,0.031624,1.00005e-06,-2.19535e-08,1.10847e-11,37667.4,19.1115], Tmin=(100,'K'), Tmax=(898.439,'K')), NASAPolynomial(coeffs=[7.3815,0.0228612,-7.26335e-06,1.16577e-09,-7.5492e-14,36206.9,-7.78214], Tmin=(898.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ)"""),
)

species(
    label = 'C[C]=C[CH]C(1192)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {5,D} {13,S}
5  C u1 p0 c0 {2,S} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (326.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,310.455],'cm^-1')),
        HinderedRotor(inertia=(0.126891,'amu*angstrom^2'), symmetry=1, barrier=(8.57706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00479519,'amu*angstrom^2'), symmetry=1, barrier=(8.59683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347395,'amu*angstrom^2'), symmetry=1, barrier=(23.5091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0954,0.0359647,-7.4598e-06,-8.15017e-09,3.70225e-12,39331.1,18.0956], Tmin=(100,'K'), Tmax=(1187.73,'K')), NASAPolynomial(coeffs=[7.89777,0.0277523,-1.13954e-05,2.08941e-09,-1.43358e-13,37153.7,-14.2596], Tmin=(1187.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[CH2](1193)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {5,S} {9,S}
4  C u1 p0 c0 {1,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (304.198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00401063,'amu*angstrom^2'), symmetry=1, barrier=(7.11751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.986217,'amu*angstrom^2'), symmetry=1, barrier=(22.6751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395285,'amu*angstrom^2'), symmetry=1, barrier=(22.6808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01079,0.0344277,7.32934e-06,-3.03227e-08,1.30186e-11,36666.4,21.0919], Tmin=(100,'K'), Tmax=(1010.22,'K')), NASAPolynomial(coeffs=[10.0278,0.0241677,-9.33579e-06,1.72607e-09,-1.22041e-13,33950.3,-23.0929], Tmin=(1010.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C[C]=CC(1194)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,D} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {3,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (390.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,216.692],'cm^-1')),
        HinderedRotor(inertia=(0.14184,'amu*angstrom^2'), symmetry=1, barrier=(4.67359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13653,'amu*angstrom^2'), symmetry=1, barrier=(4.6726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132688,'amu*angstrom^2'), symmetry=1, barrier=(4.68272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04972,0.0405227,-2.17348e-05,5.57612e-09,-5.73512e-13,47043.1,21.1549], Tmin=(100,'K'), Tmax=(2150.77,'K')), NASAPolynomial(coeffs=[12.3921,0.021288,-8.32004e-06,1.41801e-09,-9.01843e-14,42594.3,-36.6614], Tmin=(2150.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=[C]CC(1195)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {5,D}
5  C u1 p0 c0 {3,S} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (423.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.303071,'amu*angstrom^2'), symmetry=1, barrier=(6.9682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302867,'amu*angstrom^2'), symmetry=1, barrier=(6.96351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30297,'amu*angstrom^2'), symmetry=1, barrier=(6.96587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42429,0.0396229,-2.13175e-05,5.66535e-09,-6.45071e-13,50943,18.4619], Tmin=(100,'K'), Tmax=(1705.61,'K')), NASAPolynomial(coeffs=[6.54364,0.0299622,-1.28214e-05,2.3445e-09,-1.58316e-13,49537.8,-3.611], Tmin=(1705.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[CH]C(389)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u1 p0 c0 {1,S} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u1 p0 c0 {4,S} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (240.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1818,0.0283565,2.70959e-05,-5.14697e-08,2.05698e-11,28948.9,17.5847], Tmin=(100,'K'), Tmax=(990.207,'K')), NASAPolynomial(coeffs=[10.2368,0.0240426,-9.12523e-06,1.70245e-09,-1.22295e-13,25969.9,-28.1839], Tmin=(990.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCC(1196)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u1 p0 c0 {5,S} {12,S} {13,S}
5  C u1 p0 c0 {3,D} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (336.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91756,0.0385394,-7.73085e-06,-1.25885e-08,6.29119e-12,40588.1,20.1152], Tmin=(100,'K'), Tmax=(1066.88,'K')), NASAPolynomial(coeffs=[9.32924,0.0254946,-1.0119e-05,1.85674e-09,-1.28994e-13,38167.5,-20.054], Tmin=(1066.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
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
    E0 = (90.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (442.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (98.8933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (154.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (115.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (279.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (207.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (253.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (300.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (241.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (509.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (226.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (211.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (324.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (275.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (240.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (144.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['C2H4(29)', 'C#CC(211)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(17)', '[CH2]C=[C]C(1189)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['CC1=CCC1(1185)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['C=CC=CC(382)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['C=C=CCC(1120)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H4(29)', '[CH]=[C]C(212)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.80516e+17,'m^3/(mol*s)'), n=-3.08137, Ea=(46.7466,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.2182269765482571, var=280.02273246586066, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_N-Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_N-Sp-4R!H-3R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(6)', 'C=CC=[C]C(1122)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000216213,'m^3/(mol*s)'), n=3.29826, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.02158166639280058, var=0.8466995369935014, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', '[CH2]CC=C=C(1190)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.63842,'m^3/(mol*s)'), n=2.47011, Ea=(12.3908,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.16475677084270915, var=0.2152759280376351, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_N-Sp-4CCl-2C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_N-Sp-4CCl-2C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH2](39)', 'C#CC(211)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.75673e+24,'m^3/(mol*s)'), n=-5.41269, Ea=(109.874,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=2.33529332423437, var=201.40057549922756, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R_4R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[CH2]CC#CC(1191)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2440,'m^3/(mol*s)'), n=1.64, Ea=(17.131,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH2](39)', '[CH]=[C]C(212)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.83804e+08,'m^3/(mol*s)'), n=-0.65, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.007291370132393764, var=0.8826412212865885, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[C]=C[CH]C(1192)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.18e+09,'s^-1'), n=1.22, Ea=(199.995,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 157 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C=CC[CH2](1193)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C[C]=CC(1194)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[C]=[C]CC(1195)'],
    products = ['[CH2]CC=[C]C(1118)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['[CH2]C=C[CH]C(389)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC=[C]C(1118)'],
    products = ['[CH2][C]=CCC(1196)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #647',
    isomers = [
        '[CH2]CC=[C]C(1118)',
    ],
    reactants = [
        ('C2H4(29)', 'C#CC(211)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #647',
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

