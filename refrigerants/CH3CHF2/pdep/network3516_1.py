species(
    label = 'C=C(CC[CH]F)OF(10087)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {1,S} {5,S} {13,S}
8  C u0 p0 c0 {6,D} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-93.7838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,347.424,348.39,348.649,2097.81],'cm^-1')),
        HinderedRotor(inertia=(0.00139167,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130129,'amu*angstrom^2'), symmetry=1, barrier=(11.1662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130541,'amu*angstrom^2'), symmetry=1, barrier=(11.1623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249822,'amu*angstrom^2'), symmetry=1, barrier=(21.5376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634788,0.0794802,-9.12998e-05,5.98886e-08,-1.64237e-11,-11163.1,29.0231], Tmin=(100,'K'), Tmax=(873.136,'K')), NASAPolynomial(coeffs=[10.0416,0.0363853,-1.72641e-05,3.35933e-09,-2.37787e-13,-12805.8,-15.083], Tmin=(873.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.7838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    label = 'C=C[CH]F(1341)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {8,S}
4 C u0 p0 c0 {2,D} {6,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-39.4934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.618875,'amu*angstrom^2'), symmetry=1, barrier=(28.5117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2734.35,'J/mol'), sigma=(4.78057,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=427.10 K, Pc=56.79 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95088,0.00313686,8.55963e-05,-1.7635e-07,1.14205e-10,-4747,9.3909], Tmin=(10,'K'), Tmax=(483.604,'K')), NASAPolynomial(coeffs=[2.28776,0.0264346,-1.62621e-05,4.86442e-09,-5.65127e-13,-4697.72,15.0527], Tmin=(483.604,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-39.4934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CDC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45064e-06,-3.57997e-09,7.56176e-13,50400.6,-0.411759], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62651,0.00394757,-1.49921e-06,2.54533e-10,-1.62951e-14,50691.7,6.78356], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(C[CH]F)OF(607)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {1,S} {4,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-68.3779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,335.458,337.38,1313.75],'cm^-1')),
        HinderedRotor(inertia=(0.00149164,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139785,'amu*angstrom^2'), symmetry=1, barrier=(11.2594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278579,'amu*angstrom^2'), symmetry=1, barrier=(22.3753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41057,0.0599284,-6.14946e-05,3.42528e-08,-7.91388e-12,-8133.17,24.1181], Tmin=(100,'K'), Tmax=(1027.8,'K')), NASAPolynomial(coeffs=[10.0491,0.0263089,-1.2429e-05,2.427e-09,-1.72592e-13,-9908.89,-17.7945], Tmin=(1027.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.3779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-CdHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)CC(=C)OF(10095)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {6,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-102.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,386.417,386.417,1645.32],'cm^-1')),
        HinderedRotor(inertia=(0.087959,'amu*angstrom^2'), symmetry=1, barrier=(9.32007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0879589,'amu*angstrom^2'), symmetry=1, barrier=(9.32007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228869,'amu*angstrom^2'), symmetry=1, barrier=(24.2509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0879589,'amu*angstrom^2'), symmetry=1, barrier=(9.32007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3240.62,'J/mol'), sigma=(5.57238,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.18 K, Pc=42.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584694,0.0804666,-9.17011e-05,5.87902e-08,-1.56879e-11,-12222.4,28.4271], Tmin=(100,'K'), Tmax=(897.405,'K')), NASAPolynomial(coeffs=[10.6062,0.0357976,-1.70369e-05,3.32297e-09,-2.35657e-13,-14021.1,-18.8357], Tmin=(897.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'O=C(CF)CC[CH]F(10165)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  C u1 p0 c0 {2,S} {5,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-442.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420448,0.0871508,-0.000129138,1.14737e-07,-4.12864e-11,-53108.2,29.0179], Tmin=(100,'K'), Tmax=(798.78,'K')), NASAPolynomial(coeffs=[6.76218,0.0424341,-2.083e-05,4.0313e-09,-2.80954e-13,-53707.8,2.43545], Tmin=(798.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    label = '[CH]CCC(=C)OF(10166)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  C u2 p0 c0 {4,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (351.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02703,0.0661629,-6.18326e-05,3.11017e-08,-6.46029e-12,42341.3,26.5254], Tmin=(100,'K'), Tmax=(1139.15,'K')), NASAPolynomial(coeffs=[11.6172,0.0289768,-1.28671e-05,2.44567e-09,-1.71406e-13,39928.5,-25.9459], Tmin=(1139.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(230)',
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
    label = '[CH2]CC(=C)OF(602)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (131.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,956.622],'cm^-1')),
        HinderedRotor(inertia=(0.00520283,'amu*angstrom^2'), symmetry=1, barrier=(3.34741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146565,'amu*angstrom^2'), symmetry=1, barrier=(3.36983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671117,'amu*angstrom^2'), symmetry=1, barrier=(15.4303,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0882,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72886,0.049493,-4.05168e-05,1.76615e-08,-3.20419e-12,15954.6,22.5402], Tmin=(100,'K'), Tmax=(1280.65,'K')), NASAPolynomial(coeffs=[9.98563,0.0237037,-1.03103e-05,1.93696e-09,-1.34554e-13,13839.8,-19.3364], Tmin=(1280.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ)"""),
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
    label = 'C=C(CC[C]F)OF-2(10167)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  C u2 p0 c0 {1,S} {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (84.6839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768087,0.0765107,-8.90898e-05,5.80688e-08,-1.57314e-11,10296.8,27.0207], Tmin=(100,'K'), Tmax=(884.871,'K')), NASAPolynomial(coeffs=[10.2731,0.0335435,-1.62527e-05,3.19247e-09,-2.27206e-13,8614.71,-17.6729], Tmin=(884.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'FO[C]1CCC(F)C1(10117)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-161.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05974,0.0541148,-1.78033e-05,-1.24803e-08,7.60416e-12,-19347.1,25.0139], Tmin=(100,'K'), Tmax=(1067.51,'K')), NASAPolynomial(coeffs=[13.4112,0.0288456,-1.18213e-05,2.22217e-09,-1.57058e-13,-23181.4,-40.9892], Tmin=(1067.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + ring(Cyclopentane) + radical(C2CsJO)"""),
)

species(
    label = '[CH2]C1(OF)CCC1F(10136)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
8  C u1 p0 c0 {4,S} {14,S} {15,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-73.2487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261211,0.0726981,-6.08592e-05,2.55042e-08,-4.24615e-12,-8667.3,25.5508], Tmin=(100,'K'), Tmax=(1439.45,'K')), NASAPolynomial(coeffs=[17.9813,0.0234567,-9.54625e-06,1.73902e-09,-1.18652e-13,-13768.7,-66.3928], Tmin=(1439.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.2487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsOs) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs-Cs) + radical(CJC(C)2O)"""),
)

species(
    label = 'CH2CHF(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(=C)OF(605)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 C u0 p0 c0 {3,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (109.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,343.543],'cm^-1')),
        HinderedRotor(inertia=(0.153653,'amu*angstrom^2'), symmetry=1, barrier=(12.8689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0184304,'amu*angstrom^2'), symmetry=1, barrier=(33.8385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0616,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87336,0.00851378,0.000122333,-3.04869e-07,2.24407e-10,13187.1,10.7552], Tmin=(10,'K'), Tmax=(466.39,'K')), NASAPolynomial(coeffs=[4.30812,0.0285537,-1.8564e-05,5.80221e-09,-6.95957e-13,12888,6.21793], Tmin=(466.39,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(109.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[CH2]C(DC)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(CC=CF)OF(9275)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  C u0 p0 c0 {1,S} {6,D} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-159.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,194,682,905,1196,1383,3221,286.911,665.835,4000],'cm^-1')),
        HinderedRotor(inertia=(0.417869,'amu*angstrom^2'), symmetry=1, barrier=(24.4093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.417849,'amu*angstrom^2'), symmetry=1, barrier=(24.4093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170684,'amu*angstrom^2'), symmetry=1, barrier=(9.97009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3301.05,'J/mol'), sigma=(5.40021,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=515.62 K, Pc=47.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935366,0.0669879,-6.1109e-05,2.9322e-08,-5.76137e-12,-19108.1,27.204], Tmin=(100,'K'), Tmax=(1204.61,'K')), NASAPolynomial(coeffs=[12.7173,0.0278641,-1.23904e-05,2.35909e-09,-1.65477e-13,-21946.6,-31.8303], Tmin=(1204.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH)"""),
)

species(
    label = 'C=C([O])CC[CH]F(10168)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {1,S} {4,S} {12,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-113.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,247.124,247.157,247.304,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0027595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454221,'amu*angstrom^2'), symmetry=1, barrier=(19.6676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453656,'amu*angstrom^2'), symmetry=1, barrier=(19.667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919726,0.0664497,-6.38355e-05,3.28989e-08,-6.88965e-12,-13573.2,27.0175], Tmin=(100,'K'), Tmax=(1144.6,'K')), NASAPolynomial(coeffs=[12.51,0.0259451,-1.07535e-05,1.98124e-09,-1.3664e-13,-16226.4,-30.4642], Tmin=(1144.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O]F(126)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=[C]CC[CH]F(10169)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {11,S}
5  C u0 p0 c0 {6,D} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (200.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,1685,370,347.976,348.176,348.433],'cm^-1')),
        HinderedRotor(inertia=(0.00138854,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11335,'amu*angstrom^2'), symmetry=1, barrier=(9.74566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113378,'amu*angstrom^2'), symmetry=1, barrier=(9.74375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73507,0.0536195,-4.8517e-05,2.64549e-08,-6.39094e-12,24230.3,23.0176], Tmin=(100,'K'), Tmax=(953.953,'K')), NASAPolynomial(coeffs=[6.86107,0.0321244,-1.47157e-05,2.83133e-09,-1.99571e-13,23252.3,-1.47034], Tmin=(953.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-CsHH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]F(233)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (116.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0021411,'amu*angstrom^2'), symmetry=1, barrier=(6.24533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46308e-06,3.53003e-10,6.36796e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95462,0.00674548,-2.74609e-06,4.76054e-10,-2.99987e-14,11484.3,-14.7487], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(Cs_P)"""),
)

species(
    label = 'C=[C]OF(202)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (235.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2950,3100,1380,975,1025,1650,1685,370,1858.56],'cm^-1')),
        HinderedRotor(inertia=(0.190472,'amu*angstrom^2'), symmetry=1, barrier=(8.45471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27376,0.0120477,-4.21846e-06,6.03973e-11,9.31488e-14,28267.1,10.4556], Tmin=(100,'K'), Tmax=(2667.03,'K')), NASAPolynomial(coeffs=[19.9694,-0.00396407,5.52785e-07,-7.38522e-11,6.52159e-15,17217.5,-85.6832], Tmin=(2667.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=CJO)"""),
)

species(
    label = '[CH2]C[CH]F(229)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 C u1 p0 c0 {2,S} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (90.3064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,1987.25],'cm^-1')),
        HinderedRotor(inertia=(0.0822317,'amu*angstrom^2'), symmetry=1, barrier=(1.89067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0825916,'amu*angstrom^2'), symmetry=1, barrier=(1.89894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11718,0.0233784,-1.02074e-05,1.57251e-09,-3.06009e-14,10889.4,15.695], Tmin=(100,'K'), Tmax=(2177.16,'K')), NASAPolynomial(coeffs=[12.3643,0.0108442,-4.64116e-06,8.07975e-10,-5.07364e-14,5807.02,-38.5364], Tmin=(2177.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.3064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(F1s)(H)) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(=CC[CH]F)OF(9345)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {3,S} {5,D} {8,S}
7  C u1 p0 c0 {1,S} {4,S} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (54.5131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,180,180,695.647],'cm^-1')),
        HinderedRotor(inertia=(0.0456428,'amu*angstrom^2'), symmetry=1, barrier=(15.6992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108889,'amu*angstrom^2'), symmetry=1, barrier=(2.50357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683313,'amu*angstrom^2'), symmetry=1, barrier=(15.7107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00843729,'amu*angstrom^2'), symmetry=1, barrier=(47.4578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6724,0.0766643,-8.727e-05,5.43069e-08,-1.38156e-11,6673.32,29.266], Tmin=(100,'K'), Tmax=(947.371,'K')), NASAPolynomial(coeffs=[11.6512,0.0303102,-1.38772e-05,2.66101e-09,-1.87071e-13,4593.09,-23.107], Tmin=(947.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.5131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(C[CH][CH]F)OF(10170)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {4,S} {8,S} {11,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  C u1 p0 c0 {1,S} {6,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (100.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,334,575,1197,1424,3202,426.13,427.361,428.301,1845.31],'cm^-1')),
        HinderedRotor(inertia=(0.0590828,'amu*angstrom^2'), symmetry=1, barrier=(7.74913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0593367,'amu*angstrom^2'), symmetry=1, barrier=(7.75759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0598897,'amu*angstrom^2'), symmetry=1, barrier=(7.74316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100516,'amu*angstrom^2'), symmetry=1, barrier=(24.2442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08684,0.0737311,-7.71992e-05,2.26773e-08,1.59304e-11,12213.9,30.0336], Tmin=(100,'K'), Tmax=(537.688,'K')), NASAPolynomial(coeffs=[8.11121,0.0370549,-1.83463e-05,3.59635e-09,-2.53944e-13,11233.3,-1.59076], Tmin=(537.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH]=C(CC[CH]F)OF(9348)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {1,S} {5,S} {13,S}
8  C u1 p0 c0 {6,D} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (153.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,334,575,1197,1424,3202,3120,650,792.5,1650,302.641,304.89,305.281],'cm^-1')),
        HinderedRotor(inertia=(0.0651098,'amu*angstrom^2'), symmetry=1, barrier=(4.31985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0649199,'amu*angstrom^2'), symmetry=1, barrier=(4.30461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181335,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264896,'amu*angstrom^2'), symmetry=1, barrier=(17.3242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.412825,0.0848544,-0.000114954,8.71896e-08,-2.70008e-11,18562.9,30.2159], Tmin=(100,'K'), Tmax=(784.994,'K')), NASAPolynomial(coeffs=[10.5838,0.0330279,-1.59233e-05,3.08736e-09,-2.16804e-13,16966.1,-16.3913], Tmin=(784.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-CsHH)(F1s)(H)) + radical(Cds_P)"""),
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
    label = 'C=C(CC[C]F)OF(9362)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  C u0 p1 c0 {1,S} {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (80.5982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,617,898,1187,264.216,264.216,264.225],'cm^-1')),
        HinderedRotor(inertia=(0.756191,'amu*angstrom^2'), symmetry=1, barrier=(37.4637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280262,'amu*angstrom^2'), symmetry=1, barrier=(13.8838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00241472,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155972,'amu*angstrom^2'), symmetry=1, barrier=(13.8838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900054,0.0720099,-7.69841e-05,4.55865e-08,-1.12098e-11,9802.16,28.2351], Tmin=(100,'K'), Tmax=(970.241,'K')), NASAPolynomial(coeffs=[10.6538,0.0317981,-1.48161e-05,2.86981e-09,-2.03033e-13,7909.48,-18.5262], Tmin=(970.241,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.5982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'C=C(C[CH]CF)OF(10130)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {5,S} {13,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-91.3478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,307.721,307.732,307.739,307.753],'cm^-1')),
        HinderedRotor(inertia=(0.00801071,'amu*angstrom^2'), symmetry=1, barrier=(9.45141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00178029,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00177998,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245868,'amu*angstrom^2'), symmetry=1, barrier=(16.5242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969389,0.0726686,-7.57786e-05,4.63659e-08,-1.22475e-11,-10882.6,29.2411], Tmin=(100,'K'), Tmax=(889.962,'K')), NASAPolynomial(coeffs=[8.48549,0.038886,-1.88375e-05,3.71043e-09,-2.64778e-13,-12220.4,-6.14335], Tmin=(889.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.3478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C(=CCCF)OF(10143)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {3,S} {6,D} {8,S}
8  C u1 p0 c0 {7,S} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-139.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639096,0.0765955,-7.9769e-05,4.55832e-08,-1.07375e-11,-16625.3,28.1097], Tmin=(100,'K'), Tmax=(1014.83,'K')), NASAPolynomial(coeffs=[11.8198,0.0325263,-1.4631e-05,2.79242e-09,-1.96089e-13,-18894.6,-25.9951], Tmin=(1014.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(CCCF)OF(10171)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-38.7938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,3120,650,792.5,1650,201.688,201.69,201.69],'cm^-1')),
        HinderedRotor(inertia=(0.00414409,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00414401,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577479,'amu*angstrom^2'), symmetry=1, barrier=(16.6694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577433,'amu*angstrom^2'), symmetry=1, barrier=(16.6694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643576,0.0782234,-8.4288e-05,5.02529e-08,-1.24416e-11,-4548.63,28.2511], Tmin=(100,'K'), Tmax=(963.8,'K')), NASAPolynomial(coeffs=[11.1874,0.0344648,-1.6186e-05,3.14733e-09,-2.23145e-13,-6581.1,-22.2284], Tmin=(963.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.7938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])CCC(F)F(10172)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-527.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467067,0.0740901,-6.96518e-05,3.41022e-08,-6.70967e-12,-63301,28.2879], Tmin=(100,'K'), Tmax=(1221.2,'K')), NASAPolynomial(coeffs=[15.0823,0.0262185,-1.08512e-05,2.00234e-09,-1.38307e-13,-66870.6,-45.1425], Tmin=(1221.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    E0 = (-116.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (443.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (94.3562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (11.4866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (366.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (289.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (239.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-94.7606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-88.4846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-61.2004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-5.44963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-63.7294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (246.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (175.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (268.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (214.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (255.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (307.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (213.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (238.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (14.1977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-23.6326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-63.2144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-103.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['HF(38)', 'CH2CO(28)', 'C=C[CH]F(1341)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(34.9309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', 'C=C(C[CH]F)OF(607)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(150.208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['[CH2]C(F)CC(=C)OF(10095)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['O=C(CF)CC[CH]F(10165)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(162.731,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', '[CH]CCC(=C)OF(10166)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(230)', '[CH2]CC(=C)OF(602)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'C=C(CC[C]F)OF-2(10167)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['FO[C]1CCC(F)C1(10117)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.23e+07,'s^-1'), n=1, Ea=(56.484,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_HH_D;doublebond_intra_secNd_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['[CH2]C1(OF)CCC1F(10136)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH2CHF(56)', '[CH2]C(=C)OF(605)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.002829,'m^3/(mol*s)'), n=2.61489, Ea=(39.7183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4622630802209839, var=1.4932365737855162, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_N-Sp-5R!H-4R!H_Ext-2R!H-R_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_N-Sp-5R!H-4R!H_Ext-2R!H-R_Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C=C(CC=CF)OF(9275)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C=C([O])CC[CH]F(10168)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(34.625,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]F(126)', 'C=[C]CC[CH]F(10169)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]F(233)', '[CH2]C(=C)OF(605)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(6.80044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OF(202)', '[CH2]C[CH]F(229)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -22.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH2]C(=CC[CH]F)OF(9345)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(5.42832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'C=C(C[CH][CH]F)OF(10170)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[CH]=C(CC[CH]F)OF(9348)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -4.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CHF(40)', '[CH2]CC(=C)OF(602)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'C=C(CC[C]F)OF(9362)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(3.74367,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['C=C(C[CH]CF)OF(10130)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(186.798,'s^-1'), n=3.1775, Ea=(165.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['[CH2]C(=CCCF)OF(10143)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(CCCF)OF(10171)'],
    products = ['C=C(CC[CH]F)OF(10087)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 3.1622776601683795
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(CC[CH]F)OF(10087)'],
    products = ['C=C([O])CCC(F)F(10172)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(47.904,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #3516',
    isomers = [
        'C=C(CC[CH]F)OF(10087)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', 'C=C[CH]F(1341)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3516',
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

