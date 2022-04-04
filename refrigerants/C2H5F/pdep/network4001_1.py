species(
    label = '[CH]=C(CF)C[CH]F(8719)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {3,S} {12,S}
7  C u1 p0 c0 {5,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (20.9153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,334,575,1197,1424,3202,3120,650,792.5,1650,308.521,315.182],'cm^-1')),
        HinderedRotor(inertia=(0.00170623,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13907,'amu*angstrom^2'), symmetry=1, barrier=(9.94522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352902,'amu*angstrom^2'), symmetry=1, barrier=(24.4294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35705,0.0622928,-6.1872e-05,3.43929e-08,-8.10546e-12,2607.24,24.5982], Tmin=(100,'K'), Tmax=(996.892,'K')), NASAPolynomial(coeffs=[9.15803,0.030991,-1.4772e-05,2.8944e-09,-2.06129e-13,1051.92,-13.0126], Tmin=(996.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(Cds_P)"""),
)

species(
    label = 'CH2CHF(55)',
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
    label = 'C#CCF(4379)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (8.12032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([319,1023,1071,1259,1317,1409,3054,3019,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.60388,'amu*angstrom^2'), symmetry=1, barrier=(36.8763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2742,'J/mol'), sigma=(4.732,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.29 K, Pc=58.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95512,0.00280826,7.03906e-05,-1.43186e-07,9.06266e-11,979.285,8.14914], Tmin=(10,'K'), Tmax=(503.334,'K')), NASAPolynomial(coeffs=[2.71827,0.0217978,-1.34993e-05,4.0832e-09,-4.79139e-13,987.76,12.1145], Tmin=(503.334,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(8.12032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""C#CCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]=CC[CH]F(928)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u1 p0 c0 {1,S} {2,S} {9,S}
5  C u1 p0 c0 {3,D} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (239.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.285427,'amu*angstrom^2'), symmetry=1, barrier=(6.56254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.017796,'amu*angstrom^2'), symmetry=1, barrier=(9.60029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2893.3,'J/mol'), sigma=(5.10367,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=451.93 K, Pc=49.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50437,0.0364087,-3.12328e-05,1.77903e-08,-4.83487e-12,28849.6,17.8965], Tmin=(100,'K'), Tmax=(825.466,'K')), NASAPolynomial(coeffs=[4.67061,0.0259124,-1.21609e-05,2.38841e-09,-1.70625e-13,28492,7.861], Tmin=(825.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(Cds_P)"""),
)

species(
    label = 'CH2(S)(25)',
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
    label = '[CH]=C(F)C[CH]F(2561)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u1 p0 c0 {2,S} {3,S} {9,S}
6  C u1 p0 c0 {4,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (43.0035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,246,474,533,1155,334,575,1197,1424,3202,3120,650,792.5,1650,180,1194.57],'cm^-1')),
        HinderedRotor(inertia=(0.0749775,'amu*angstrom^2'), symmetry=1, barrier=(1.72388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0773715,'amu*angstrom^2'), symmetry=1, barrier=(1.77892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0713,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2911.53,'J/mol'), sigma=(5.01807,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=454.77 K, Pc=52.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15142,0.0449382,-4.75799e-05,3.09966e-08,-8.88676e-12,5234.95,20.6129], Tmin=(100,'K'), Tmax=(817.822,'K')), NASAPolynomial(coeffs=[6.04574,0.0258907,-1.26439e-05,2.51751e-09,-1.80902e-13,4597.98,2.60827], Tmin=(817.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C(CF)C([CH2])F(8755)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  C u1 p0 c0 {5,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (22.7327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.289418,'amu*angstrom^2'), symmetry=1, barrier=(6.65429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289281,'amu*angstrom^2'), symmetry=1, barrier=(6.65115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992575,'amu*angstrom^2'), symmetry=1, barrier=(22.8212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09172,0.0699952,-8.84028e-05,6.65061e-08,-2.11968e-11,2833.33,24.803], Tmin=(100,'K'), Tmax=(752.926,'K')), NASAPolynomial(coeffs=[7.82671,0.0342158,-1.71241e-05,3.39525e-09,-2.42158e-13,1819.11,-5.7783], Tmin=(752.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(H)(H)) + radical(Cds_P)"""),
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
    label = '[CH]C(=C)CF(8267)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {8,S} {9,S}
5  C u2 p0 c0 {3,S} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (157.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,452.601,452.824,462.097,462.988],'cm^-1')),
        HinderedRotor(inertia=(0.351513,'amu*angstrom^2'), symmetry=1, barrier=(51.8374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344079,'amu*angstrom^2'), symmetry=1, barrier=(51.7383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24297,0.0347638,-1.39344e-05,4.52096e-10,6.29281e-13,19059.5,17.8685], Tmin=(100,'K'), Tmax=(1493.19,'K')), NASAPolynomial(coeffs=[8.99165,0.0236942,-9.85524e-06,1.77445e-09,-1.18442e-13,16262.8,-20.0121], Tmin=(1493.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'FCC1=CC(F)C1(11141)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-278.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68611,0.0418107,-7.03715e-06,-1.54311e-08,7.11556e-12,-33400.6,20.3897], Tmin=(100,'K'), Tmax=(1120.59,'K')), NASAPolynomial(coeffs=[11.2301,0.0264975,-1.16435e-05,2.24443e-09,-1.59746e-13,-36717.1,-31.9953], Tmin=(1120.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd)"""),
)

species(
    label = 'C=C(C=CF)CF(8720)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {5,D} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-310.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8937,0.0606342,-5.18555e-05,2.24692e-08,-3.86995e-12,-37281.7,22.256], Tmin=(100,'K'), Tmax=(1395.95,'K')), NASAPolynomial(coeffs=[15.2326,0.0195472,-7.70609e-06,1.38475e-09,-9.39499e-14,-41285,-51.7039], Tmin=(1395.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCFH)"""),
)

species(
    label = '[CH]=[C]CF(4380)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {4,D}
4 C u1 p0 c0 {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (311.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,377.76,379.216,382.52,1900.36,1902.09,1902.61,1902.73,3346.48],'cm^-1')),
        HinderedRotor(inertia=(0.11674,'amu*angstrom^2'), symmetry=1, barrier=(12.0012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74483,0.0317272,-5.02124e-05,4.76765e-08,-1.76459e-11,37561.4,13.6835], Tmin=(100,'K'), Tmax=(834.526,'K')), NASAPolynomial(coeffs=[4.0464,0.0168537,-7.95764e-06,1.52224e-09,-1.05087e-13,37644.8,9.44116], Tmin=(834.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=CF)CF(11142)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {2,S} {5,D} {11,S}
7  C u1 p0 c0 {4,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-63.8645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.899093,'amu*angstrom^2'), symmetry=1, barrier=(20.6719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897373,'amu*angstrom^2'), symmetry=1, barrier=(20.6324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0319,0.0617593,-6.06178e-05,3.03236e-08,-6.02638e-12,-7571.3,22.1569], Tmin=(100,'K'), Tmax=(1217.49,'K')), NASAPolynomial(coeffs=[14.0481,0.0189955,-7.93115e-06,1.47388e-09,-1.02379e-13,-10740.7,-43.2002], Tmin=(1217.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.8645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]F(839)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46308e-06,3.53001e-10,6.36801e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95459,0.00674552,-2.7461e-06,4.76057e-10,-2.99989e-14,11484.3,-14.7486], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(Cs_P)"""),
)

species(
    label = 'CH2F(46)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-42.5685,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(33.0141,'amu')),
        NonlinearRotor(inertia=([1.91548,16.2277,17.9803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([576.418,1180.5,1217.62,1485.55,3118.23,3268.88],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03338,-0.00262849,2.74227e-05,-3.89096e-08,1.85259e-11,-5119.82,5.20374], Tmin=(10,'K'), Tmax=(594.366,'K')), NASAPolynomial(coeffs=[2.59024,0.00857266,-4.60348e-06,1.22743e-09,-1.29255e-13,-4974.57,11.194], Tmin=(594.366,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-42.5685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[CH2]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C#CC[CH]F(932)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3 C u1 p0 c0 {1,S} {2,S} {8,S}
4 C u0 p0 c0 {2,S} {5,T}
5 C u0 p0 c0 {4,T} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (163.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.056649,'amu*angstrom^2'), symmetry=1, barrier=(7.75292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60745,'amu*angstrom^2'), symmetry=1, barrier=(36.9584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0728,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85366,0.0151478,0.000182468,-9.31572e-07,1.41579e-09,19717.1,9.81645], Tmin=(10,'K'), Tmax=(227.72,'K')), NASAPolynomial(coeffs=[4.21994,0.027048,-1.6686e-05,5.02114e-09,-5.85397e-13,19652.9,7.54758], Tmin=(227.72,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(163.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C#CC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(=CCF)CF(11143)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {12,S}
7  C u2 p0 c0 {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-49.2732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59182,0.0587375,-4.58124e-05,2.14283e-08,-4.73262e-12,-5844.33,22.7963], Tmin=(100,'K'), Tmax=(984.145,'K')), NASAPolynomial(coeffs=[5.77307,0.0417431,-1.99101e-05,3.8818e-09,-2.75331e-13,-6667.32,2.69098], Tmin=(984.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.2732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C[CH]F)CF(6956)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {7,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-137.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38214,0.0582969,-4.7124e-05,1.95606e-08,-3.37169e-12,-16407.9,22.171], Tmin=(100,'K'), Tmax=(1334.02,'K')), NASAPolynomial(coeffs=[11.5,0.027959,-1.30112e-05,2.51293e-09,-1.76873e-13,-19107.4,-29.5575], Tmin=(1334.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(=CF)C[CH]F(8514)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {7,D}
5  C u1 p0 c0 {1,S} {3,S} {10,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {4,D} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-90.9722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,228.357,1052.75],'cm^-1')),
        HinderedRotor(inertia=(0.570697,'amu*angstrom^2'), symmetry=1, barrier=(21.1108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135039,'amu*angstrom^2'), symmetry=1, barrier=(5.19399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555203,'amu*angstrom^2'), symmetry=1, barrier=(43.2623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3148.36,'J/mol'), sigma=(5.41053,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=491.77 K, Pc=45.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05117,0.0610607,-5.27246e-05,2.34526e-08,-4.21759e-12,-10832.2,24.7675], Tmin=(100,'K'), Tmax=(1320.93,'K')), NASAPolynomial(coeffs=[13.505,0.0233487,-9.90057e-06,1.83965e-09,-1.27138e-13,-14122.4,-38.7813], Tmin=(1320.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=CF)CCF(8722)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {12,S}
7  C u2 p0 c0 {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-65.5188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,194,682,905,1196,1383,3221,356.076,356.076,356.076,356.076,356.076],'cm^-1')),
        HinderedRotor(inertia=(0.568802,'amu*angstrom^2'), symmetry=1, barrier=(51.1767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568802,'amu*angstrom^2'), symmetry=1, barrier=(51.1767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568801,'amu*angstrom^2'), symmetry=1, barrier=(51.1767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18124,0.0615028,-4.67683e-05,1.89618e-08,-3.24801e-12,-7778.46,23.8518], Tmin=(100,'K'), Tmax=(1331.11,'K')), NASAPolynomial(coeffs=[10.8097,0.0325693,-1.41638e-05,2.63231e-09,-1.81131e-13,-10341.8,-25.3539], Tmin=(1331.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.5188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)CC(F)F(6434)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  C u2 p0 c0 {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-95.2334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,2950,3100,1380,975,1025,1650,350.685,350.689,350.689,350.693,350.713],'cm^-1')),
        HinderedRotor(inertia=(0.582534,'amu*angstrom^2'), symmetry=1, barrier=(50.8358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582522,'amu*angstrom^2'), symmetry=1, barrier=(50.8357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582497,'amu*angstrom^2'), symmetry=1, barrier=(50.8357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11656,0.0619794,-4.70209e-05,1.87348e-08,-3.12344e-12,-11349.2,24.0167], Tmin=(100,'K'), Tmax=(1373.43,'K')), NASAPolynomial(coeffs=[11.6115,0.0314137,-1.36382e-05,2.53075e-09,-1.73863e-13,-14232,-29.9455], Tmin=(1373.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.2334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-51.2667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (386.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (541.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (109.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (300.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-42.9823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (12.1335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (91.5202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (75.7583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (71.1143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (69.7999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (355.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (224.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (57.1453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (140.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (93.9183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (18.3241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (136.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['CH2CHF(55)', 'C#CCF(4379)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[CH]=CC[CH]F(928)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.44768e+28,'m^3/(mol*s)'), n=-6.4458, Ea=(80.7879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.6591608693425184, var=5.4995193120720405, Tref=1000.0, N=19, data_mean=0.0, correlation='CH',), comment="""Estimated from node CH"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', '[CH]=C(F)C[CH]F(2561)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(151.292,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(CF)C([CH2])F(8755)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(837)', '[CH]C(=C)CF(8267)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['FCC1=CC(F)C1(11141)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_1H;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['C=C(C=CF)CF(8720)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CHF(55)', '[CH]=[C]CF(4380)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(4.79179,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH]=C(C=CF)CF(11142)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]F(839)', 'C#CCF(4379)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(18.9765,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2F(46)', 'C#CC[CH]F(932)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(20.6103,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]F(839)', '[CH]=[C]CF(4380)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -19.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHF(40)', '[CH]C(=C)CF(8267)'],
    products = ['[CH]=C(CF)C[CH]F(8719)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['[CH]C(=CCF)CF(11143)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.23666e+09,'s^-1'), n=1.04335, Ea=(108.412,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['[CH2]C(=C[CH]F)CF(6956)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['[CH2]C(=CF)C[CH]F(8514)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['[CH]C(=CF)CCF(8722)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.67145,'s^-1'), n=3.19634, Ea=(69.5908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_1H;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(CF)C[CH]F(8719)'],
    products = ['[CH]C(=C)CC(F)F(6434)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(187.99,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #4001',
    isomers = [
        '[CH]=C(CF)C[CH]F(8719)',
    ],
    reactants = [
        ('CH2CHF(55)', 'C#CCF(4379)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4001',
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

