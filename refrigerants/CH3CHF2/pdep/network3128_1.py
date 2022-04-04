species(
    label = '[CH]=CC1[CH]C1F(9070)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u1 p0 c0 {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (418.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,1133.74,1133.78,1133.79,1133.82,1133.83],'cm^-1')),
        HinderedRotor(inertia=(0.245038,'amu*angstrom^2'), symmetry=1, barrier=(5.63391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30918,0.029499,4.58749e-06,-2.17922e-08,8.80503e-12,50422.6,21.7322], Tmin=(100,'K'), Tmax=(1076.9,'K')), NASAPolynomial(coeffs=[8.86719,0.0223887,-9.53359e-06,1.8224e-09,-1.29713e-13,48010,-15.036], Tmin=(1076.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(23)',
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
    label = 'FC1C2C=CC12(9095)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (105.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49368,0.0171045,5.40038e-05,-8.18845e-08,3.19833e-11,12780.6,13.466], Tmin=(100,'K'), Tmax=(975.899,'K')), NASAPolynomial(coeffs=[12.5614,0.0159676,-5.92823e-06,1.19217e-09,-9.27284e-14,8904.74,-44.6499], Tmin=(975.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_4_ene_1)"""),
)

species(
    label = 'C=CC1=CC1F(9096)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {8,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (117.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83812,0.0357189,5.3241e-06,-3.64356e-08,1.74322e-11,14185,16.7211], Tmin=(100,'K'), Tmax=(968.239,'K')), NASAPolynomial(coeffs=[14.1119,0.0135724,-4.60998e-06,8.67712e-10,-6.52526e-14,10469.5,-49.0095], Tmin=(968.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cd-Cs(F)-Cd)"""),
)

species(
    label = 'C=CC1C=C1F(9071)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {2,S} {3,D} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {8,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (182.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01185,0.0380547,-2.06897e-05,3.4162e-09,4.54876e-13,21986.3,17.5782], Tmin=(100,'K'), Tmax=(1248.32,'K')), NASAPolynomial(coeffs=[9.7414,0.0205686,-8.42803e-06,1.54075e-09,-1.05367e-13,19489.2,-23.6993], Tmin=(1248.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = '[CH]C=CC=CF(7661)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,D} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {6,S} {9,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u2 p0 c0 {4,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (236.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,194,682,905,1196,1383,3221,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13607,'amu*angstrom^2'), symmetry=1, barrier=(49.1124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14321,'amu*angstrom^2'), symmetry=1, barrier=(49.2767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58409,0.0441825,-1.15193e-05,-1.68531e-08,9.73732e-12,28574.1,20.4358], Tmin=(100,'K'), Tmax=(988.152,'K')), NASAPolynomial(coeffs=[12.0433,0.0219236,-8.21074e-06,1.47854e-09,-1.0316e-13,25526.7,-34.8595], Tmin=(988.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[CH](264)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843366,2.04091e-05,-2.5194e-08,9.47559e-12,71059.3,4.72117], Tmin=(100,'K'), Tmax=(935.371,'K')), NASAPolynomial(coeffs=[4.92143,0.00358269,-9.24415e-07,1.57267e-10,-1.19536e-14,70476.2,-2.30666], Tmin=(935.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CC1=CC1F(9097)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {8,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (364.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,282.763,283.18,283.245,283.599,283.737],'cm^-1')),
        HinderedRotor(inertia=(0.513742,'amu*angstrom^2'), symmetry=1, barrier=(29.4191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7644,0.0392423,-1.13421e-05,-1.90308e-08,1.1501e-11,43904.7,17.3888], Tmin=(100,'K'), Tmax=(968.403,'K')), NASAPolynomial(coeffs=[14.3511,0.0107409,-3.57661e-06,6.69071e-10,-5.03938e-14,40365.5,-48.6169], Tmin=(968.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cd-Cs(F)-Cd) + radical(Cds_P)"""),
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
    label = '[CH]=CC1C=C1(9098)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {1,S} {2,D} {8,S}
4  C u0 p0 c0 {1,S} {5,D} {9,S}
5  C u1 p0 c0 {4,D} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (597.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180,229.809,900.899,901.105,904.396,2545.3],'cm^-1')),
        HinderedRotor(inertia=(0.0034164,'amu*angstrom^2'), symmetry=1, barrier=(1.97481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0931,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28089,0.0319613,-1.04946e-05,-8.14908e-09,5.17931e-12,71984.2,15.4561], Tmin=(100,'K'), Tmax=(1014.21,'K')), NASAPolynomial(coeffs=[9.28063,0.0167425,-6.30769e-06,1.14197e-09,-7.95259e-14,69927.2,-21.5533], Tmin=(1014.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1C=C1F(9099)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {2,S} {3,D} {8,S}
5  C u0 p0 c0 {2,S} {6,D} {9,S}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (429.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,1154.35,1154.68,2195.41],'cm^-1')),
        HinderedRotor(inertia=(0.219039,'amu*angstrom^2'), symmetry=1, barrier=(5.03614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14227,0.0393027,-2.99644e-05,1.19566e-08,-1.97555e-12,51696.9,17.504], Tmin=(100,'K'), Tmax=(1397.65,'K')), NASAPolynomial(coeffs=[9.2927,0.0188386,-8.00168e-06,1.48064e-09,-1.01696e-13,49698.1,-19.3866], Tmin=(1397.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cd(F)-Cd) + radical(Cds_P)"""),
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
    label = 'C#CC1[CH]C1F(9100)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,T}
6  C u0 p0 c0 {5,T} {10,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (336.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2175,525,750,770,3400,2100,292.685,819.854,823.65,827.565,829.123,830.371],'cm^-1')),
        HinderedRotor(inertia=(0.00523741,'amu*angstrom^2'), symmetry=1, barrier=(2.52365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40202,0.0293025,-2.88239e-06,-1.38051e-08,6.56965e-12,40538,18.8001], Tmin=(100,'K'), Tmax=(1022.95,'K')), NASAPolynomial(coeffs=[8.04207,0.0198965,-7.63641e-06,1.38016e-09,-9.5464e-14,38722.3,-11.7725], Tmin=(1022.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane)"""),
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
    label = '[CH]=C[C]1C=C1(9101)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {1,S} {2,D} {7,S}
4 C u0 p0 c0 {1,S} {5,D} {8,S}
5 C u1 p0 c0 {4,D} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (729.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180,1053.24,1053.3,1053.3,1053.3],'cm^-1')),
        HinderedRotor(inertia=(0.00668586,'amu*angstrom^2'), symmetry=1, barrier=(5.26321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0851,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56856,0.0264513,-4.65221e-06,-1.13808e-08,6.02513e-12,87846.8,12.6766], Tmin=(100,'K'), Tmax=(986.555,'K')), NASAPolynomial(coeffs=[7.98618,0.0158294,-5.74985e-06,1.01616e-09,-6.98423e-14,86225.8,-16.1848], Tmin=(986.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Allyl_T) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]1CC1F(7765)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u1 p0 c0 {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (324.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43436,0.0234297,2.83216e-05,-5.04765e-08,2.00321e-11,39124,17.8784], Tmin=(100,'K'), Tmax=(987.319,'K')), NASAPolynomial(coeffs=[9.86162,0.0202428,-7.71052e-06,1.45274e-09,-1.05393e-13,36346.1,-24.4997], Tmin=(987.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(Allyl_T) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1C[C]1F(9102)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u1 p0 c0 {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (383.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19604,0.0328574,-5.00768e-06,-1.22085e-08,5.58319e-12,46163.6,22.4955], Tmin=(100,'K'), Tmax=(1110.04,'K')), NASAPolynomial(coeffs=[9.0732,0.0220785,-9.36416e-06,1.77203e-09,-1.24856e-13,43774.1,-15.2867], Tmin=(1110.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C1[CH]C1F(8780)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (409.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2950,3100,1380,975,1025,1650,1685,370,180,180,180,1233.24,1233.28,1233.59,1233.59],'cm^-1')),
        HinderedRotor(inertia=(0.00285399,'amu*angstrom^2'), symmetry=1, barrier=(3.07976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34261,0.0305819,-3.16973e-06,-1.09787e-08,4.47581e-12,49306.7,21.688], Tmin=(100,'K'), Tmax=(1183.47,'K')), NASAPolynomial(coeffs=[8.05562,0.0237444,-1.03112e-05,1.94899e-09,-1.3616e-13,47081,-10.5263], Tmin=(1183.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C1CC1F(9103)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {6,D}
6  C u1 p0 c0 {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (430.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,250,417,511,1155,1315,1456,3119,1685,370,3120,650,792.5,1650,180,180,1038.58,1038.83,1039.05,1041.14,1041.75],'cm^-1')),
        HinderedRotor(inertia=(0.00676002,'amu*angstrom^2'), symmetry=1, barrier=(5.17813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10845,0.0335788,-2.3313e-06,-1.81299e-08,8.37476e-12,51865,20.5771], Tmin=(100,'K'), Tmax=(1048.18,'K')), NASAPolynomial(coeffs=[10.1351,0.0203809,-8.39208e-06,1.592e-09,-1.13503e-13,49224.6,-23.093], Tmin=(1048.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C1[CH]C1F(8684)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {2,S} {3,S} {8,S}
5  C u0 p0 c0 {3,D} {6,S} {9,S}
6  C u1 p0 c0 {5,S} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (257.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2950,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,672.272,672.275,672.276,672.278,672.283],'cm^-1')),
        HinderedRotor(inertia=(0.000372998,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16392,0.0293292,1.67852e-05,-4.08921e-08,1.69291e-11,31056.8,18.0024], Tmin=(100,'K'), Tmax=(1004.49,'K')), NASAPolynomial(coeffs=[11.2291,0.0195252,-7.84059e-06,1.51223e-09,-1.1063e-13,27909.1,-32.3758], Tmin=(1004.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd(Cd)-Cs-Cs(F)) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC1[CH][C]1F(9104)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u1 p0 c0 {1,S} {2,S} {4,S}
4  C u1 p0 c0 {2,S} {3,S} {8,S}
5  C u0 p0 c0 {2,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (362.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46403,0.029406,-4.01973e-06,-7.71511e-09,2.91769e-12,43603.9,23.489], Tmin=(100,'K'), Tmax=(1310.14,'K')), NASAPolynomial(coeffs=[7.56739,0.0245966,-1.08462e-05,2.03411e-09,-1.40142e-13,41342.2,-6.03843], Tmin=(1310.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(cyclopropane)"""),
)

species(
    label = 'FC=CC1[CH][CH]1(9105)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u1 p0 c0 {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u0 p0 c0 {1,S} {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (374.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51809,0.024007,1.77141e-05,-3.48522e-08,1.34364e-11,45085.6,22.501], Tmin=(100,'K'), Tmax=(1030.09,'K')), NASAPolynomial(coeffs=[8.36678,0.0217924,-8.90781e-06,1.69382e-09,-1.21291e-13,42793.2,-11.1673], Tmin=(1030.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cyclopropane) + radical(cyclopropane) + radical(cyclopropane)"""),
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
    E0 = (97.3219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (105.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (160.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (122.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (146.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (322.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (260.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (393.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (319.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (220.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (226.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (578.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (374.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (258.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (258.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (202.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (250.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (289.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (141.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (303.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['C2H2(23)', 'FC1C=C1(6880)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['FC1C2C=CC12(9095)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['C=CC1=CC1F(9096)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['C=CC1C=C1F(9071)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C=CC=CF(7661)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[CH](264)', 'FC1C=C1(6880)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0205056,'m^3/(mol*s)'), n=2.42326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7837337452902932, var=8.455757054931007, Tref=1000.0, N=46, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', '[CH]=CC1=CC1F(9097)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(5.26325,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', '[CH]=CC1C=C1(9098)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(44.3293,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH]=CC1C=C1F(9099)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C2H2(23)', 'FC1[CH][CH]1(8666)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.1632,'m^3/(mol*s)'), n=2.24, Ea=(15.6832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C#CC1[CH]C1F(9100)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[CH](264)', 'FC1[CH][CH]1(8666)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', '[CH]=C[C]1C=C1(9101)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(247.032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['[CH]=C[C]1CC1F(7765)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['[CH]=CC1C[C]1F(9102)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['C=[C]C1[CH]C1F(8780)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C1CC1F(9103)'],
    products = ['[CH]=CC1[CH]C1F(9070)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_23cy3;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['[CH2]C=C1[CH]C1F(8684)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['C=CC1[CH][C]1F(9104)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC1[CH]C1F(9070)'],
    products = ['FC=CC1[CH][CH]1(9105)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(206.198,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

network(
    label = 'PDepNetwork #3128',
    isomers = [
        '[CH]=CC1[CH]C1F(9070)',
    ],
    reactants = [
        ('C2H2(23)', 'FC1C=C1(6880)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3128',
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

