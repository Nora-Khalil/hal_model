species(
    label = 'CC([O])(F)C[CH]F(5144)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-323.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,334,437,422,648,895,1187,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.202804,'amu*angstrom^2'), symmetry=1, barrier=(4.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203138,'amu*angstrom^2'), symmetry=1, barrier=(4.67054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203239,'amu*angstrom^2'), symmetry=1, barrier=(4.67287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902227,0.0737471,-0.00010127,8.21527e-08,-2.76684e-11,-38789.9,24.7025], Tmin=(100,'K'), Tmax=(739.093,'K')), NASAPolynomial(coeffs=[8.24481,0.0328604,-1.59592e-05,3.10003e-09,-2.17621e-13,-39843.9,-8.28916], Tmin=(739.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(O2sj(Cs-F1sCsCs)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    label = 'CC(=O)F(499)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-453.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926],'cm^-1')),
        HinderedRotor(inertia=(0.163766,'amu*angstrom^2'), symmetry=1, barrier=(3.76529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3044.49,'J/mol'), sigma=(4.92747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.54 K, Pc=57.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95746,0.00861016,1.67316e-05,-2.39571e-08,8.75423e-12,-54563.9,8.39893], Tmin=(10,'K'), Tmax=(963.489,'K')), NASAPolynomial(coeffs=[3.63904,0.0176469,-9.3479e-06,2.39871e-09,-2.40789e-13,-54860.6,8.06499], Tmin=(963.489,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-453.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C[CH]F(1986)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
6  C u1 p0 c0 {2,S} {4,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-281.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.350077,'amu*angstrom^2'), symmetry=1, barrier=(8.04895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351447,'amu*angstrom^2'), symmetry=1, barrier=(8.08046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3363.16,'J/mol'), sigma=(5.63685,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.32 K, Pc=42.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33407,0.0638981,-0.000100712,8.57603e-08,-2.86915e-11,-33766.9,20.7484], Tmin=(100,'K'), Tmax=(832.244,'K')), NASAPolynomial(coeffs=[8.567,0.0213618,-1.00368e-05,1.90348e-09,-1.3057e-13,-34701.6,-11.2009], Tmin=(832.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)C(C)([O])F(5142)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-317.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,334,437,422,648,895,1187,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.190603,'amu*angstrom^2'), symmetry=1, barrier=(4.38233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189579,'amu*angstrom^2'), symmetry=1, barrier=(4.35879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.573661,'amu*angstrom^2'), symmetry=1, barrier=(13.1896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04898,0.069958,-8.33188e-05,5.66193e-08,-1.60421e-11,-38090.8,23.7819], Tmin=(100,'K'), Tmax=(847.431,'K')), NASAPolynomial(coeffs=[9.21494,0.031414,-1.50947e-05,2.94871e-09,-2.08989e-13,-39474.9,-14.2624], Tmin=(847.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(O2sj(Cs-F1sCsCs)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[C](F)C[CH]F(886)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  C u1 p0 c0 {2,S} {3,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-163.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,334,575,1197,1424,3202,180,2237.08],'cm^-1')),
        HinderedRotor(inertia=(0.32658,'amu*angstrom^2'), symmetry=1, barrier=(7.50872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.326769,'amu*angstrom^2'), symmetry=1, barrier=(7.51306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32622,'amu*angstrom^2'), symmetry=1, barrier=(7.50045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0872,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36119,0.0670235,-0.000109268,1.04372e-07,-3.87059e-11,-19607.1,21.7813], Tmin=(100,'K'), Tmax=(834.684,'K')), NASAPolynomial(coeffs=[4.12958,0.0346958,-1.69178e-05,3.2532e-09,-2.2476e-13,-19405.3,12.9031], Tmin=(834.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(Cs-HHH)(F1s)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    label = '[CH2]C(C)([O])F(1498)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-94.4025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859729,'amu*angstrom^2'), symmetry=1, barrier=(1.97669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865757,'amu*angstrom^2'), symmetry=1, barrier=(1.99055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0696,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90274,0.0518036,-8.00086e-05,7.24717e-08,-2.5814e-11,-11284,18.2012], Tmin=(100,'K'), Tmax=(838.692,'K')), NASAPolynomial(coeffs=[5.28502,0.024564,-1.14231e-05,2.16143e-09,-1.48215e-13,-11460.6,4.80775], Tmin=(838.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.4025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCsCs)) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'CC1(F)CC(F)O1(5147)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {6,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-611.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65962,0.0404686,2.00586e-06,-2.73088e-08,1.16906e-11,-73425.7,18.925], Tmin=(100,'K'), Tmax=(1071.31,'K')), NASAPolynomial(coeffs=[12.2055,0.0255071,-1.12289e-05,2.19882e-09,-1.5914e-13,-77086.3,-39.2179], Tmin=(1071.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + group(Cs-CsHHH) + ring(O2s-Cs-Cs-Cs(F))"""),
)

species(
    label = 'CC(O)(F)C=CF(5200)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {6,D} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-625.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23107,0.0635112,-6.24283e-05,3.30594e-08,-7.26182e-12,-75164.7,20.9653], Tmin=(100,'K'), Tmax=(1076.85,'K')), NASAPolynomial(coeffs=[10.7512,0.0281484,-1.31698e-05,2.56407e-09,-1.82093e-13,-77215,-25.6688], Tmin=(1076.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CC(=O)C[CH]F(1210)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u1 p0 c0 {1,S} {3,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-251.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,334,575,1197,1424,3202,282.122,287.998],'cm^-1')),
        HinderedRotor(inertia=(0.0796702,'amu*angstrom^2'), symmetry=1, barrier=(4.59594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0852569,'amu*angstrom^2'), symmetry=1, barrier=(4.76147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0177968,'amu*angstrom^2'), symmetry=1, barrier=(15.9307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0882,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01645,0.0468123,-3.59546e-05,1.5075e-08,-2.76124e-12,-30217.5,20.2906], Tmin=(100,'K'), Tmax=(1214.43,'K')), NASAPolynomial(coeffs=[7.63575,0.0283039,-1.30942e-05,2.52576e-09,-1.77892e-13,-31582.4,-7.91104], Tmin=(1214.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(Csj(Cs-COHH)(F1s)(H))"""),
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
    label = 'CH3(19)',
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
    label = 'O=C(F)C[CH]F(2002)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {4,S} {9,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-461.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,486,617,768,1157,1926,180,2379.27],'cm^-1')),
        HinderedRotor(inertia=(0.510628,'amu*angstrom^2'), symmetry=1, barrier=(11.7403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78978,'amu*angstrom^2'), symmetry=1, barrier=(41.1505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53445,0.042379,-7.36494e-05,9.74719e-08,-5.58992e-11,-55523.9,11.9088], Tmin=(10,'K'), Tmax=(484.303,'K')), NASAPolynomial(coeffs=[4.88444,0.0277354,-1.74745e-05,5.24973e-09,-6.04863e-13,-55613.7,6.79761], Tmin=(484.303,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-461.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""ODC(F)C[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C[C]([O])F(1344)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-67.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,4000],'cm^-1')),
        HinderedRotor(inertia=(0.378092,'amu*angstrom^2'), symmetry=1, barrier=(8.69308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91286,0.0277959,-4.27461e-05,4.09524e-08,-1.51725e-11,-8085.5,14.0056], Tmin=(100,'K'), Tmax=(854.905,'K')), NASAPolynomial(coeffs=[3.47688,0.0163951,-7.36921e-06,1.37672e-09,-9.36531e-14,-7861.75,13.2456], Tmin=(854.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'CC([O])(F)C=CF(5201)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {6,D} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-395.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23047,'amu*angstrom^2'), symmetry=1, barrier=(28.2909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23767,'amu*angstrom^2'), symmetry=1, barrier=(28.4566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3894,0.0596035,-5.69521e-05,2.87113e-08,-5.98535e-12,-47464.1,20.2026], Tmin=(100,'K'), Tmax=(1130.46,'K')), NASAPolynomial(coeffs=[10.8822,0.026015,-1.23845e-05,2.4289e-09,-1.73116e-13,-49610.3,-26.7588], Tmin=(1130.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-395.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(CC(C)OJ)"""),
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
    label = 'CC([O])=C[CH]F(1177)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-150.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,180],'cm^-1')),
        HinderedRotor(inertia=(0.0257446,'amu*angstrom^2'), symmetry=1, barrier=(11.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493684,'amu*angstrom^2'), symmetry=1, barrier=(11.3508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63579,0.0487383,-4.26136e-05,1.94177e-08,-3.57407e-12,-18060.9,20.1357], Tmin=(100,'K'), Tmax=(1295.87,'K')), NASAPolynomial(coeffs=[11.4466,0.0184548,-7.55952e-06,1.38389e-09,-9.49594e-14,-20603.6,-29.7386], Tmin=(1295.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = 'C=C([O])C[CH]F(1114)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {1,S} {3,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-88.3792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,397.777,397.811,397.873],'cm^-1')),
        HinderedRotor(inertia=(0.00106553,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11138,'amu*angstrom^2'), symmetry=1, barrier=(12.4991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3507.39,'J/mol'), sigma=(5.84493,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=547.85 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38616,0.0503588,-4.53015e-05,2.08297e-08,-3.7768e-12,-10529.5,23.2349], Tmin=(100,'K'), Tmax=(1340.87,'K')), NASAPolynomial(coeffs=[13.4055,0.0145039,-5.19188e-06,8.87851e-10,-5.87482e-14,-13752.8,-38.2771], Tmin=(1340.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.3792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-CdHH)(F1s)(H))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'CC([O])(F)[CH]CF(5202)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {6,S} {13,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-315.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2800,2850,1350,1500,750,1050,1375,1000,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,323.729,323.73,323.732],'cm^-1')),
        HinderedRotor(inertia=(0.00160851,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107463,'amu*angstrom^2'), symmetry=1, barrier=(7.99181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107462,'amu*angstrom^2'), symmetry=1, barrier=(7.99181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36955,0.0613619,-5.97726e-05,3.21617e-08,-7.29029e-12,-37866.7,24.3872], Tmin=(100,'K'), Tmax=(1036.82,'K')), NASAPolynomial(coeffs=[9.56981,0.0297256,-1.40033e-05,2.7323e-09,-1.94195e-13,-39567.1,-15.4708], Tmin=(1036.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(O2sj(Cs-F1sCsCs)) + radical(CCJCO)"""),
)

species(
    label = 'CC(O)(F)[CH][CH]F(5203)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {7,S} {11,S}
7  C u1 p0 c0 {2,S} {6,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-357.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04347,0.068837,-7.63417e-05,4.58004e-08,-1.12951e-11,-42888,25.5589], Tmin=(100,'K'), Tmax=(971.565,'K')), NASAPolynomial(coeffs=[10.94,0.028092,-1.34349e-05,2.6348e-09,-1.87724e-13,-44811,-21.9005], Tmin=(971.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(CCJCO) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(O)(F)C[CH]F(3319)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {5,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
"""),
    E0 = (-346.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,334,575,1197,1424,3202,211.774,213.075],'cm^-1')),
        HinderedRotor(inertia=(0.195264,'amu*angstrom^2'), symmetry=1, barrier=(6.20793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194458,'amu*angstrom^2'), symmetry=1, barrier=(6.2032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194954,'amu*angstrom^2'), symmetry=1, barrier=(6.20584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.944275,'amu*angstrom^2'), symmetry=1, barrier=(30.314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388961,0.0868976,-0.000136157,1.17372e-07,-4.0261e-11,-41551.6,25.6483], Tmin=(100,'K'), Tmax=(806.082,'K')), NASAPolynomial(coeffs=[9.5294,0.0319042,-1.58915e-05,3.07712e-09,-2.13877e-13,-42712.1,-14.5368], Tmin=(806.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(Csj(Cs-F1sO2sCs)(H)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C([O])(F)CCF(1528)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-304.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,316,385,515,654,689,1295,528,1116,1182,1331,1402,1494,3075,3110,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.199598,'amu*angstrom^2'), symmetry=1, barrier=(4.58916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199329,'amu*angstrom^2'), symmetry=1, barrier=(4.58296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.515615,'amu*angstrom^2'), symmetry=1, barrier=(11.855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703193,0.0795224,-0.000119786,1.03859e-07,-3.62949e-11,-36529.7,24.522], Tmin=(100,'K'), Tmax=(797.445,'K')), NASAPolynomial(coeffs=[7.97665,0.0338678,-1.66588e-05,3.22317e-09,-2.24474e-13,-37398.1,-7.09367], Tmin=(797.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(O2sj(Cs-F1sCsCs)) + radical(Csj(Cs-F1sO2sCs)(H)(H))"""),
)

species(
    label = 'C[C](C[CH]F)OF(1958)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {1,S} {4,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (1.10871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,334,575,1197,1424,3202,282.7,282.706,2871.21],'cm^-1')),
        HinderedRotor(inertia=(0.138317,'amu*angstrom^2'), symmetry=1, barrier=(7.8455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138325,'amu*angstrom^2'), symmetry=1, barrier=(7.84556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138322,'amu*angstrom^2'), symmetry=1, barrier=(7.84558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614051,'amu*angstrom^2'), symmetry=1, barrier=(34.8272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679153,0.0808329,-0.000126028,1.11132e-07,-3.88298e-11,245.415,26.5517], Tmin=(100,'K'), Tmax=(824.608,'K')), NASAPolynomial(coeffs=[7.67025,0.0337476,-1.64161e-05,3.14194e-09,-2.16712e-13,-459.701,-3.11246], Tmin=(824.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.10871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(C2CsJO) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'C[C]([O])CC(F)F(4705)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {4,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-320.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,274.115,275.879,275.904],'cm^-1')),
        HinderedRotor(inertia=(0.174385,'amu*angstrom^2'), symmetry=1, barrier=(9.52842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179812,'amu*angstrom^2'), symmetry=1, barrier=(9.52539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218344,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971675,0.0715165,-8.92564e-05,6.3786e-08,-1.88703e-11,-38397.7,23.8027], Tmin=(100,'K'), Tmax=(816.144,'K')), NASAPolynomial(coeffs=[9.25814,0.0309066,-1.46245e-05,2.82735e-09,-1.98832e-13,-39750.4,-14.4918], Tmin=(816.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsCsFFH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
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
    E0 = (-176.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (284.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-12.3674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (225.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (267.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-168.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-113.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-10.7707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-126.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-153.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-67.6315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-37.0426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (195.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (35.1277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (37.3762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (190.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-43.7656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-101.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-176.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-74.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (221.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (52.4594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])(F)C[CH]F(5144)'],
    products = ['CH2CHF(55)', 'CC(=O)F(499)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[O]C(F)C[CH]F(1986)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(486642,'m^3/(mol*s)'), n=0.225671, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)C(C)([O])F(5142)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'C[C](F)C[CH]F(886)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(837)', '[CH2]C(C)([O])F(1498)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC([O])(F)C[CH]F(5144)'],
    products = ['CC1(F)CC(F)O1(5147)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC([O])(F)C[CH]F(5144)'],
    products = ['CC(O)(F)C=CF(5200)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'CC(=O)C[CH]F(1210)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(21.5933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]F(839)', 'CC(=O)F(499)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(64.7682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(19)', 'O=C(F)C[CH]F(2002)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.000120113,'m^3/(mol*s)'), n=2.79634, Ea=(25.4492,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_Ext-1R!H-R_4R!H->F_3C-u1',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_Ext-1R!H-R_4R!H->F_3C-u1"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2CHF(55)', 'C[C]([O])F(1344)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00126,'m^3/(mol*s)'), n=2.41, Ea=(6.38187,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_4R!H->C_Sp-2C=1C_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_4R!H->C_Sp-2C=1C_Ext-3C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'CC([O])(F)C=CF(5201)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]F(839)', 'C[C]([O])F(1344)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'CC([O])=C[CH]F(1177)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(320.58,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'C=C([O])C[CH]F(1114)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(260.311,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHF(40)', '[CH2]C(C)([O])F(1498)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC([O])(F)[CH]CF(5202)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00689639,'s^-1'), n=4.185, Ea=(125.282,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])(F)C[CH]F(5144)'],
    products = ['CC(O)(F)[CH][CH]F(5203)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(F)C[CH]F(3319)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.249e+09,'s^-1'), n=0.385, Ea=(23.1,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 443 used for R3H_SS_Cs;C_rad_out_2H;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 22.2 to 23.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])(F)CCF(1528)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(262028,'s^-1'), n=1.73125, Ea=(84.0042,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[C](C[CH]F)OF(1958)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(74.2675,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C[C]([O])CC(F)F(4705)'],
    products = ['CC([O])(F)C[CH]F(5144)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(226.029,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1860',
    isomers = [
        'CC([O])(F)C[CH]F(5144)',
    ],
    reactants = [
        ('CH2CHF(55)', 'CC(=O)F(499)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1860',
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

