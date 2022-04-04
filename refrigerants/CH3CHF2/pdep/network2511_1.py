species(
    label = 'C=C([O])C=[C]C(6342)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (219.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.62738,'amu*angstrom^2'), symmetry=1, barrier=(14.4247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626539,'amu*angstrom^2'), symmetry=1, barrier=(14.4054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26113,0.0566476,-5.64874e-05,3.00478e-08,-6.3653e-12,26519.4,20.8486], Tmin=(100,'K'), Tmax=(1149.08,'K')), NASAPolynomial(coeffs=[12.2279,0.018472,-6.65322e-06,1.1353e-09,-7.49594e-14,23999,-33.5835], Tmin=(1149.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C#CC(5603)',
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.64],'cm^-1')),
        HinderedRotor(inertia=(0.630825,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30522,0.0109266,1.31984e-05,-2.25157e-08,8.87389e-12,20738.2,7.19167], Tmin=(100,'K'), Tmax=(969.992,'K')), NASAPolynomial(coeffs=[5.80059,0.0113535,-4.03487e-06,7.19054e-10,-5.02177e-14,19749.9,-7.36988], Tmin=(969.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C=[C]C(6461)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,D} {5,S} {9,S}
3  C u0 p0 c0 {5,D} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {2,D}
5  C u1 p0 c0 {2,S} {3,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (494.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(1.02157,'amu*angstrom^2'), symmetry=1, barrier=(23.488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02245,'amu*angstrom^2'), symmetry=1, barrier=(23.508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92364,0.0436709,-3.71914e-05,1.81261e-08,-3.67553e-12,59579,17.3236], Tmin=(100,'K'), Tmax=(1170.44,'K')), NASAPolynomial(coeffs=[8.65616,0.0206623,-7.7042e-06,1.33056e-09,-8.80665e-14,58003,-16.2164], Tmin=(1170.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1C=C(C)O1(6489)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {4,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {3,D} {4,S} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (66.9629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09119,0.0321592,8.56317e-06,-3.26107e-08,1.42801e-11,8131.19,15.3081], Tmin=(100,'K'), Tmax=(995.857,'K')), NASAPolynomial(coeffs=[10.9371,0.0197137,-7.4627e-06,1.39529e-09,-1.00399e-13,5224.63,-33.0788], Tmin=(995.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.9629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Oxetene)"""),
)

species(
    label = 'C=C=CC(=C)O(6352)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (20.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00849,0.0551578,-3.13435e-05,-1.05157e-08,1.15277e-11,2595.99,18.8355], Tmin=(100,'K'), Tmax=(926.311,'K')), NASAPolynomial(coeffs=[17.9641,0.00942042,-1.7798e-06,2.34015e-10,-1.70725e-14,-1724.24,-68.0317], Tmin=(926.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C[C]=C[C]1CO1(6517)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {2,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (382.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57696,0.0483526,-4.31194e-05,2.32378e-08,-5.02899e-12,46054.2,22.4697], Tmin=(100,'K'), Tmax=(1271.13,'K')), NASAPolynomial(coeffs=[8.67984,0.0213428,-5.74932e-06,7.5542e-10,-4.02291e-14,44624.8,-12.0213], Tmin=(1271.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Cds_S)"""),
)

species(
    label = 'C[C]1C=C([O])C1(6518)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {3,S} {6,S}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {4,S} {5,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (167.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35678,0.0184828,6.19909e-05,-9.79295e-08,4.03204e-11,20274.2,17.0094], Tmin=(100,'K'), Tmax=(931.536,'K')), NASAPolynomial(coeffs=[14.0046,0.0135928,-2.79796e-06,4.39834e-10,-3.54717e-14,16146.3,-48.8675], Tmin=(931.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]C1([O])C=C1C(6508)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {10,S}
6  C u1 p0 c0 {2,S} {11,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (447.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24582,0.05429,-4.48456e-05,1.85237e-08,-3.05153e-12,53917.2,19.6139], Tmin=(100,'K'), Tmax=(1447.12,'K')), NASAPolynomial(coeffs=[14.1979,0.0184887,-7.73567e-06,1.4276e-09,-9.80291e-14,50168.6,-47.6592], Tmin=(1447.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = 'C=C=CC(=C)[O](6519)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (158.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13088,'amu*angstrom^2'), symmetry=1, barrier=(26.0011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976072,0.0558197,-5.47946e-05,2.69001e-08,-5.05273e-12,19170.5,20.4568], Tmin=(100,'K'), Tmax=(1432.27,'K')), NASAPolynomial(coeffs=[15.4409,0.0107348,-2.66824e-06,3.52093e-10,-1.99491e-14,15507.8,-52.8458], Tmin=(1432.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=[C]C(5604)',
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
        HinderedRotor(inertia=(0.0859241,'amu*angstrom^2'), symmetry=1, barrier=(1.97556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04067e-06,1.20986e-09,-3.01512e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.69,'K')), NASAPolynomial(coeffs=[6.96312,0.0100051,-3.7071e-06,6.32809e-10,-4.05654e-14,57370.1,-10.4646], Tmin=(1807.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=O(136)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C([O])C#CC(6520)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {6,T}
6  C u0 p0 c0 {3,S} {5,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (152.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08153,'amu*angstrom^2'), symmetry=1, barrier=(24.8666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57099,'amu*angstrom^2'), symmetry=1, barrier=(36.1201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73389,0.0434681,-2.31597e-05,-5.30301e-09,6.84063e-12,18421.5,19.1879], Tmin=(100,'K'), Tmax=(923.178,'K')), NASAPolynomial(coeffs=[11.9877,0.0155676,-4.68152e-06,7.46383e-10,-4.93388e-14,15824,-33.276], Tmin=(923.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C=CC(=C)[O](6521)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,D} {4,S} {8,S}
3  C u0 p0 c0 {2,D} {5,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (99.8633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,301.17,301.577],'cm^-1')),
        HinderedRotor(inertia=(0.348285,'amu*angstrom^2'), symmetry=1, barrier=(22.5317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15482,'amu*angstrom^2'), symmetry=1, barrier=(74.4997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60393,0.0427774,-7.06754e-06,-2.71099e-08,1.55237e-11,12106.2,20.7035], Tmin=(100,'K'), Tmax=(925.169,'K')), NASAPolynomial(coeffs=[13.8391,0.015431,-4.15964e-06,6.48301e-10,-4.42371e-14,8748.7,-43.2826], Tmin=(925.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C([O])[C]=CC(6522)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {6,D} {10,S}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  C u1 p0 c0 {3,D} {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (180.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25,0.0581805,-5.93187e-05,3.17209e-08,-6.34328e-12,21846.7,20.0126], Tmin=(100,'K'), Tmax=(918.938,'K')), NASAPolynomial(coeffs=[11.598,0.0195302,-6.66454e-06,1.09228e-09,-7.02862e-14,19674.9,-30.5045], Tmin=(918.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C(O)[C]=[C]C(6523)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {12,S}
2  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {6,D}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (280.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.970892,'amu*angstrom^2'), symmetry=1, barrier=(22.3227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.971298,'amu*angstrom^2'), symmetry=1, barrier=(22.332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970786,'amu*angstrom^2'), symmetry=1, barrier=(22.3203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89414,0.0662015,-7.44397e-05,4.07451e-08,-7.7449e-12,33891.1,20.2651], Tmin=(100,'K'), Tmax=(872.049,'K')), NASAPolynomial(coeffs=[14.0485,0.0158441,-4.98753e-06,7.73839e-10,-4.82197e-14,31217.4,-43.5718], Tmin=(872.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)C=[C]C(6524)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {11,S}
2  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {6,D}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {3,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (328.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.852926,'amu*angstrom^2'), symmetry=1, barrier=(19.6104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.850499,'amu*angstrom^2'), symmetry=1, barrier=(19.5547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.849567,'amu*angstrom^2'), symmetry=1, barrier=(19.5332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635451,0.0663373,-7.33077e-05,4.03043e-08,-8.50385e-12,39690.1,21.9957], Tmin=(100,'K'), Tmax=(1218.92,'K')), NASAPolynomial(coeffs=[16.6925,0.0114635,-3.09605e-06,4.35212e-10,-2.56207e-14,35937.7,-57.9842], Tmin=(1218.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=CC(6525)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u1 p0 c0 {5,D} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (228.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.730971,'amu*angstrom^2'), symmetry=1, barrier=(16.8065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.731578,'amu*angstrom^2'), symmetry=1, barrier=(16.8204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01625,0.057943,-5.64881e-05,2.8466e-08,-5.60934e-12,27644.7,21.6588], Tmin=(100,'K'), Tmax=(1270.71,'K')), NASAPolynomial(coeffs=[14.4178,0.0148873,-4.63657e-06,7.23922e-10,-4.53766e-14,24309,-45.9302], Tmin=(1270.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC(=C)O(6526)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u1 p0 c0 {6,S} {10,S} {11,S}
6  C u1 p0 c0 {3,D} {5,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (199.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21845,0.0511783,-2.36727e-05,-1.59099e-08,1.30676e-11,24151.9,21.0603], Tmin=(100,'K'), Tmax=(915.418,'K')), NASAPolynomial(coeffs=[16.3454,0.0116446,-2.42252e-06,3.15274e-10,-2.0936e-14,20269.3,-56.6613], Tmin=(915.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=CC[C]=O(6346)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u1 p0 c0 {3,S} {4,D}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (268.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,375.019],'cm^-1')),
        HinderedRotor(inertia=(0.107017,'amu*angstrom^2'), symmetry=1, barrier=(10.6546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109231,'amu*angstrom^2'), symmetry=1, barrier=(10.6431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107272,'amu*angstrom^2'), symmetry=1, barrier=(10.6515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.61,'J/mol'), sigma=(5.8495,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.13 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96819,0.0421755,-2.52846e-05,6.64221e-09,-6.69285e-13,32385,22.0889], Tmin=(100,'K'), Tmax=(2314.41,'K')), NASAPolynomial(coeffs=[19.081,0.0125998,-6.11645e-06,1.12089e-09,-7.28852e-14,24463.7,-74.8311], Tmin=(2314.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[C]=C[C]=O(6527)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3 C u0 p0 c0 {4,D} {5,S} {9,S}
4 C u1 p0 c0 {2,S} {3,D}
5 C u1 p0 c0 {1,D} {3,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (281.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.732973,'amu*angstrom^2'), symmetry=1, barrier=(16.8525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.733633,'amu*angstrom^2'), symmetry=1, barrier=(16.8677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98528,0.0509116,-8.52265e-05,8.20597e-08,-3.10125e-11,33874.5,17.4691], Tmin=(100,'K'), Tmax=(804.997,'K')), NASAPolynomial(coeffs=[4.58402,0.0251455,-1.32652e-05,2.62986e-09,-1.85053e-13,33872.6,8.08204], Tmin=(804.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC1=CC(=O)C1(6471)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {1,D} {2,S} {6,S}
6  C u0 p0 c0 {4,D} {5,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-0.0136073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45851,0.0226891,3.2671e-05,-4.90191e-08,1.68379e-11,63.1711,15.5769], Tmin=(100,'K'), Tmax=(1114.14,'K')), NASAPolynomial(coeffs=[9.60024,0.027029,-1.35353e-05,2.78137e-09,-2.04963e-13,-3388.95,-28.0003], Tmin=(1114.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.0136073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C=CC(C)=O(6353)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,D}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,D} {2,S} {4,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (14.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39174,0.0393619,-2.13869e-05,5.05012e-09,-4.56278e-13,1747.43,18.031], Tmin=(100,'K'), Tmax=(2530.23,'K')), NASAPolynomial(coeffs=[17.9004,0.0148439,-6.8515e-06,1.22024e-09,-7.78565e-14,-6100.46,-71.1857], Tmin=(2530.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C1[CH][C](C)O1(6528)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (225.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43605,0.0393706,2.16356e-05,-7.02133e-08,3.43057e-11,27226.2,15.4093], Tmin=(100,'K'), Tmax=(908.38,'K')), NASAPolynomial(coeffs=[19.1002,0.00673011,9.91525e-07,-3.54737e-10,2.31069e-14,22154.6,-78.364], Tmin=(908.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C2CsJOC(O)) + radical(C=CCJCO)"""),
)

species(
    label = 'C[C]=C=C(C)[O](6529)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {3,S} {6,D}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (271.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,1685,370,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(0.488802,'amu*angstrom^2'), symmetry=1, barrier=(11.2385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489848,'amu*angstrom^2'), symmetry=1, barrier=(11.2626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55337,0.0563875,-6.37837e-05,4.24019e-08,-1.17223e-11,32771.6,20.6527], Tmin=(100,'K'), Tmax=(871.158,'K')), NASAPolynomial(coeffs=[8.22218,0.0257652,-1.10534e-05,2.04674e-09,-1.40651e-13,31609.7,-10.6002], Tmin=(871.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C(C)[O](6530)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {6,S} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (180.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987207,0.057773,-5.63128e-05,2.87614e-08,-5.71049e-12,21781.2,21.7692], Tmin=(100,'K'), Tmax=(1331.23,'K')), NASAPolynomial(coeffs=[14.0151,0.0150176,-4.06926e-06,5.61249e-10,-3.20618e-14,18632.5,-43.6088], Tmin=(1331.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
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
    E0 = (-6.20955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (511.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (2.07477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (18.7637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (180.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (131.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (221.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (156.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (204.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (123.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (157.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (425.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (81.7358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (190.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (185.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (295.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (406.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (53.6844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (288.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (436.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (2.07477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (18.7637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (131.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (197.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (47.9732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['CH2CO(28)', 'C#CC(5603)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', 'C=[C]C=[C]C(6461)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=C1C=C(C)O1(6489)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=C=CC(=C)O(6352)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C[C]=C[C]1CO1(6517)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C[C]1C=C([O])C1(6518)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['[CH2]C1([O])C=C1C(6508)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.69117e+12,'s^-1'), n=0.1675, Ea=(227.777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_2H;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra_2H;radadd_intra_cd_Cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 227.6 to 227.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', 'C=C=CC(=C)[O](6519)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.63842,'m^3/(mol*s)'), n=2.47011, Ea=(12.3908,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.16475677084270915, var=0.2152759280376351, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_N-Sp-4CCl-2C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_N-Sp-4CCl-2C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2CO(28)', '[CH]=[C]C(5604)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.25713e-13,'m^3/(mol*s)'), n=5.31456, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.372338763425909, var=4.180287361148109, Tref=1000.0, N=384, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=O(136)', 'C#CC(5603)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0645,'m^3/(mol*s)'), n=2.41, Ea=(17.0232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-2R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-2R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C=C([O])C#CC(6520)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(763,'m^3/(mol*s)'), n=1.64, Ea=(19.5327,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(136)', '[CH]=[C]C(5604)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C=CC(=C)[O](6521)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=C([O])[C]=CC(6522)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)[C]=[C]C(6523)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.85843e+07,'s^-1'), n=1.27338, Ea=(130.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(O)C=[C]C(6524)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C([O])C=CC(6525)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['[CH2][C]=CC(=C)O(6526)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[C]=CC[C]=O(6346)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(T)(18)', 'C[C]=C[C]=O(6527)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['CC1=CC(=O)C1(6471)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=C=CC(C)=O(6353)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=C1[CH][C](C)O1(6528)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[C]=C=C(C)[O](6529)'],
    products = ['C=C([O])C=[C]C(6342)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])C=[C]C(6342)'],
    products = ['C=[C]C=C(C)[O](6530)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2511',
    isomers = [
        'C=C([O])C=[C]C(6342)',
    ],
    reactants = [
        ('CH2CO(28)', 'C#CC(5603)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2511',
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

