species(
    label = '[O]C(=O)C([O])(F)CF(4283)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-576.782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,365.689,365.845,365.864,365.91,365.937,366.064,366.5,366.907],'cm^-1')),
        HinderedRotor(inertia=(0.0012574,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520569,'amu*angstrom^2'), symmetry=1, barrier=(49.5176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25413,0.067488,-0.000106149,9.76911e-08,-3.64209e-11,-69278.9,23.6923], Tmin=(100,'K'), Tmax=(769.148,'K')), NASAPolynomial(coeffs=[6.17364,0.0318373,-1.69908e-05,3.39638e-09,-2.40995e-13,-69737.9,3.18542], Tmin=(769.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-576.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(CCOJ)"""),
)

species(
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C(F)CF(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(=O)C([O])F(991)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-365.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,180,180,1403.72,1403.76,1403.82,1404.04],'cm^-1')),
        HinderedRotor(inertia=(0.207779,'amu*angstrom^2'), symmetry=1, barrier=(4.77724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4306.78,'J/mol'), sigma=(6.22203,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.71 K, Pc=40.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43783,0.0452706,-9.27281e-05,1.00362e-07,-3.93626e-11,-43879.1,17.7905], Tmin=(100,'K'), Tmax=(857.506,'K')), NASAPolynomial(coeffs=[0.391352,0.0277341,-1.46777e-05,2.85094e-09,-1.95986e-13,-42532.4,33.1549], Tmin=(857.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + radical(O2sj(Cs-F1sCOH)) + radical(CCOJ)"""),
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
    label = '[O]C(=O)C([O])(F)F(4260)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
"""),
    E0 = (-548.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,180,180,180,199.969,1311.18,1311.2,1311.22,1311.22,1311.26,1311.27],'cm^-1')),
        HinderedRotor(inertia=(0.0460789,'amu*angstrom^2'), symmetry=1, barrier=(56.2195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4314.47,'J/mol'), sigma=(6.43546,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=673.91 K, Pc=36.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.71603,0.0199589,-8.12262e-06,3.65602e-10,1.33898e-13,-66021.9,11.5686], Tmin=(100,'K'), Tmax=(2713.47,'K')), NASAPolynomial(coeffs=[36.7775,-0.0132893,2.50971e-06,-3.4324e-10,2.384e-14,-88580.8,-184.621], Tmin=(2713.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsOs) + radical(O2sj(Cs-F1sF1sCO)) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CF)=C([O])OF(4331)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {3,S} {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-293.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,311.866,312.319,312.862,2694.92],'cm^-1')),
        HinderedRotor(inertia=(0.188761,'amu*angstrom^2'), symmetry=1, barrier=(13.0242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459161,'amu*angstrom^2'), symmetry=1, barrier=(31.8607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626483,0.0812754,-0.000137318,1.18422e-07,-3.96336e-11,-35131.1,25.4929], Tmin=(100,'K'), Tmax=(834.693,'K')), NASAPolynomial(coeffs=[10.7285,0.0220061,-1.12933e-05,2.18092e-09,-1.50143e-13,-36439.2,-19.1523], Tmin=(834.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C](F)C(=O)OCF(2014)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {7,D}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {8,S}
8  C u1 p0 c0 {2,S} {5,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-648.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61005,0.0631257,-0.000111134,1.10985e-07,-4.19741e-11,-77906.2,25.0542], Tmin=(100,'K'), Tmax=(848.605,'K')), NASAPolynomial(coeffs=[2.43437,0.0347023,-1.7519e-05,3.3672e-09,-2.31339e-13,-77162.6,26.4184], Tmin=(848.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-648.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHHO) + group(Cds-OdCsOs) + radical(O2sj(Cs-F1sCOH)) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C(=O)[C](F)CF(532)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6 C u1 p0 c0 {2,S} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-446.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,234,347,1316,1464,180,180,180,912.417,912.428,912.526,912.547],'cm^-1')),
        HinderedRotor(inertia=(0.909528,'amu*angstrom^2'), symmetry=1, barrier=(54.4468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0921229,'amu*angstrom^2'), symmetry=1, barrier=(54.4447,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75198,0.0346705,4.71953e-06,-9.18392e-08,9.51983e-11,-53703.2,18.4967], Tmin=(100,'K'), Tmax=(411.608,'K')), NASAPolynomial(coeffs=[4.34519,0.029713,-1.55704e-05,3.14807e-09,-2.27235e-13,-53923.5,11.1416], Tmin=(411.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O]C(F)([C]=O)CF(4332)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-376.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17749,'amu*angstrom^2'), symmetry=1, barrier=(27.0728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18103,'amu*angstrom^2'), symmetry=1, barrier=(27.1542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34264,0.0624342,-9.28453e-05,7.15729e-08,-2.18622e-11,-45231.4,22.7275], Tmin=(100,'K'), Tmax=(803.644,'K')), NASAPolynomial(coeffs=[10.4147,0.0172816,-8.57213e-06,1.66702e-09,-1.16733e-13,-46689.6,-19.0573], Tmin=(803.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1OOC1(F)CF(4287)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-709.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46877,0.0373099,2.29966e-05,-6.29313e-08,2.74922e-11,-85240.4,21.2261], Tmin=(100,'K'), Tmax=(989.843,'K')), NASAPolynomial(coeffs=[18.7872,0.0115738,-5.05692e-06,1.12421e-09,-9.21135e-14,-90836.6,-73.0983], Tmin=(989.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-709.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(CsCCFO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + longDistanceInteraction_cyclic(Cs(F)-CO) + ring(Cyclobutane)"""),
)

species(
    label = '[O]C(F)(CF)[C]1OO1(4333)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-221.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62745,0.0549613,-5.88166e-05,3.29745e-08,-7.56578e-12,-26527.7,25.4836], Tmin=(100,'K'), Tmax=(1040.38,'K')), NASAPolynomial(coeffs=[10.3129,0.0215683,-1.06719e-05,2.12407e-09,-1.52615e-13,-28335,-16.7625], Tmin=(1040.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(O2s-O2s-Cs(C-F)) + radical(O2sj(Cs-F1sCsCs)) + radical(Cs_P)"""),
)

species(
    label = '[O][C]1OOC1(F)CF(4334)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-243.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72099,0.0561322,-6.95432e-05,5.14786e-08,-1.65213e-11,-29251.3,24.2445], Tmin=(100,'K'), Tmax=(739.808,'K')), NASAPolynomial(coeffs=[6.60056,0.0297524,-1.60628e-05,3.29104e-09,-2.39319e-13,-29973.3,2.1734], Tmin=(739.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C1([O])OC1(F)CF(4324)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
8  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-436.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63466,0.0622156,-0.000107047,1.0918e-07,-4.28145e-11,-52465.1,22.6952], Tmin=(100,'K'), Tmax=(820.539,'K')), NASAPolynomial(coeffs=[1.95458,0.0375426,-1.96897e-05,3.87487e-09,-2.71015e-13,-51739.5,25.9564], Tmin=(820.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsOs) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]C(=O)C(=O)CF(4335)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {6,D}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6 C u0 p0 c0 {2,D} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-447.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,375,552.5,462.5,1710,292.391,292.396,292.398,292.401,781.915,1376.86,4000],'cm^-1')),
        HinderedRotor(inertia=(0.434039,'amu*angstrom^2'), symmetry=1, barrier=(26.3336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434078,'amu*angstrom^2'), symmetry=1, barrier=(26.3336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88321,0.0540514,-6.69751e-05,3.10706e-08,4.85411e-12,-53722.1,20.945], Tmin=(100,'K'), Tmax=(531.471,'K')), NASAPolynomial(coeffs=[7.35326,0.0239322,-1.31548e-05,2.67965e-09,-1.92882e-13,-54459.6,-3.45542], Tmin=(531.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-447.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ)"""),
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
    label = '[O]C(=O)C(=O)F(1357)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
"""),
    E0 = (-483.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,619,818,1246,1924,378.786,378.802,378.932,379.061,379.161,2122.23],'cm^-1')),
        HinderedRotor(inertia=(0.00117326,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50218,0.0367996,-5.2221e-05,3.85633e-08,-1.1598e-11,-58127.1,16.9862], Tmin=(100,'K'), Tmax=(804.936,'K')), NASAPolynomial(coeffs=[7.32144,0.0128509,-7.59241e-06,1.60056e-09,-1.17879e-13,-58902.9,-5.21817], Tmin=(804.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-483.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO) + radical(C=OC=OOJ)"""),
)

species(
    label = '[O][C]=O(517)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.00186761,1.03198e-05,-1.52369e-08,5.8052e-12,3804.5,8.40409], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.36181,0.000422681,-4.06569e-07,1.52493e-10,-1.51968e-14,2816.74,-6.4394], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[O][C](F)CF(882)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-241.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.77787,'amu*angstrom^2'), symmetry=1, barrier=(17.8848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43603,0.0332746,-3.37965e-05,1.70828e-08,-3.40461e-12,-28931.2,16.3061], Tmin=(100,'K'), Tmax=(1215.45,'K')), NASAPolynomial(coeffs=[9.73799,0.00924406,-4.14003e-06,8.16317e-10,-5.88299e-14,-30706.2,-20.3463], Tmin=(1215.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F2(78)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(=O)C([O])=O(4336)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,D} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-71.7875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,180,180,180,180,913.768,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00415723,'amu*angstrom^2'), symmetry=1, barrier=(2.46304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0492705,'amu*angstrom^2'), symmetry=1, barrier=(29.1886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01751,0.0479637,-6.98525e-05,5.60541e-08,-1.84278e-11,-8566.69,20.8696], Tmin=(100,'K'), Tmax=(738.791,'K')), NASAPolynomial(coeffs=[7.46474,0.0184721,-9.9768e-06,2.02576e-09,-1.45833e-13,-9371.59,-3.76126], Tmin=(738.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.7875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(CJCC=O)"""),
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
    label = '[O]C(=O)C(=O)[CH]F(4337)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,D}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,D} {6,S} {7,S}
6 C u1 p0 c0 {1,S} {5,S} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-288.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,235,1215,1347,1486,3221,358.102,358.351,358.362,358.366,358.463,1312.23,1312.64],'cm^-1')),
        HinderedRotor(inertia=(0.00131302,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295876,'amu*angstrom^2'), symmetry=1, barrier=(27.0003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74505,0.0547435,-8.28963e-05,6.78809e-08,-2.26517e-11,-34592.5,21.1954], Tmin=(100,'K'), Tmax=(729.212,'K')), NASAPolynomial(coeffs=[8.10483,0.0198599,-1.11449e-05,2.28778e-09,-1.65438e-13,-35520.1,-7.47894], Tmin=(729.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(CsCOF1sH)"""),
)

species(
    label = '[O]C(=O)C(O)(F)[CH]F(4338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u1 p0 c0 {2,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-625.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.877537,0.077539,-0.000132697,1.2325e-07,-4.49402e-11,-75098.7,24.9416], Tmin=(100,'K'), Tmax=(802.005,'K')), NASAPolynomial(coeffs=[7.25128,0.0306355,-1.67042e-05,3.33231e-09,-2.34819e-13,-75634.9,-1.37126], Tmin=(802.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CsCsF1sH)"""),
)

species(
    label = '[O]C(F)([CH]F)C(=O)O(4339)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {5,D} {6,S}
8  C u1 p0 c0 {2,S} {6,S} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-607.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77558,0.0770693,-0.000122845,1.03954e-07,-3.50828e-11,-72924,25.9246], Tmin=(100,'K'), Tmax=(771.314,'K')), NASAPolynomial(coeffs=[10.2599,0.0241631,-1.27199e-05,2.51555e-09,-1.77031e-13,-74276.4,-16.6515], Tmin=(771.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-607.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(CsCsF1sH)"""),
)

species(
    label = '[O]C([O])=C(CF)OF(4340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-289.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,357.98,363.56,366.836,1195.24],'cm^-1')),
        HinderedRotor(inertia=(0.109448,'amu*angstrom^2'), symmetry=1, barrier=(10.4598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343029,'amu*angstrom^2'), symmetry=1, barrier=(31.0604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9221,0.0748761,-0.000111919,8.38457e-08,-2.32018e-11,-34701.7,24.3675], Tmin=(100,'K'), Tmax=(636.94,'K')), NASAPolynomial(coeffs=[10.5716,0.0234776,-1.2542e-05,2.50887e-09,-1.78244e-13,-36117.6,-19.298], Tmin=(636.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(F)(OF)C([O])=O(4341)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-300.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54229,0.0657458,-6.29142e-05,-2.65944e-08,6.59939e-11,-36036,21.5497], Tmin=(100,'K'), Tmax=(473.405,'K')), NASAPolynomial(coeffs=[8.35264,0.0296145,-1.6276e-05,3.27004e-09,-2.319e-13,-36920.8,-8.74768], Tmin=(473.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])(F)C(=O)OF(4342)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {5,D} {6,S}
8  C u1 p0 c0 {6,S} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-282.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835202,0.0755586,-0.00011131,8.59986e-08,-2.66543e-11,-33835.6,24.5409], Tmin=(100,'K'), Tmax=(788.176,'K')), NASAPolynomial(coeffs=[11.0508,0.0237134,-1.26399e-05,2.53891e-09,-1.81428e-13,-35445.9,-22.3117], Tmin=(788.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(CJCO)"""),
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
    E0 = (-299.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (75.1541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (297.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (168.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-172.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (73.1313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (143.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-291.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (55.6659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (33.0754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-159.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-87.2107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-234.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-205.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-282.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (67.4335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (205.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-33.7679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-215.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-205.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (72.1057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (89.3163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (76.6682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['CO2(14)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[O]C(=O)C([O])F(991)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.79957e-05,'m^3/(mol*s)'), n=3.10993, Ea=(24.6887,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', '[O]C(=O)C([O])(F)F(4260)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(149.502,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(CF)=C([O])OF(4331)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(184.907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O][C](F)C(=O)OCF(2014)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(126.993,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', '[O]C(=O)[C](F)CF(532)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', '[O]C(F)([C]=O)CF(4332)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['O=C1OOC1(F)CF(4287)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O]C(F)(CF)[C]1OO1(4333)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(355.527,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 354.3 to 355.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O][C]1OOC1(F)CF(4334)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(332.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 331.8 to 332.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O]C1([O])OC1(F)CF(4324)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(139.938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 138.9 to 139.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(=O)C(=O)CF(4335)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(10.2198,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2F(46)', '[O]C(=O)C(=O)F(1357)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4291e-06,'m^3/(mol*s)'), n=3.20779, Ea=(14.7579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=O(517)', 'O=C(F)CF(879)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(98.2606,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CO2(14)', '[O][C](F)CF(882)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.14506e-05,'m^3/(mol*s)'), n=2.93163, Ea=(84.7045,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.051016206943592955, var=0.2066755265279035, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=O(517)', '[O][C](F)CF(882)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(78)', '[CH2]C(=O)C([O])=O(4336)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(8.77517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[O]C(=O)C(=O)[CH]F(4337)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(258.678,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O]C(=O)C(O)(F)[CH]F(4338)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(18000,'s^-1'), n=2.287, Ea=(84.6842,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(=O)C([O])(F)CF(4283)'],
    products = ['[O]C(F)([CH]F)C(=O)O(4339)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.31548e+07,'s^-1'), n=1.28185, Ea=(93.8674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C([O])=C(CF)OF(4340)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(84.5763,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)(OF)C([O])=O(4341)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(112.651,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])(F)C(=O)OF(4342)'],
    products = ['[O]C(=O)C([O])(F)CF(4283)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(81.9745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1053',
    isomers = [
        '[O]C(=O)C([O])(F)CF(4283)',
    ],
    reactants = [
        ('CO2(14)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1053',
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

