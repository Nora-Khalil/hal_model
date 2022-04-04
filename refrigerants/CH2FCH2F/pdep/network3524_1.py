species(
    label = 'FCC(F)OOC=C(F)F(5256)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-837.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,207.984,207.984],'cm^-1')),
        HinderedRotor(inertia=(0.560321,'amu*angstrom^2'), symmetry=1, barrier=(17.1998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21247,'amu*angstrom^2'), symmetry=1, barrier=(37.2183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21247,'amu*angstrom^2'), symmetry=1, barrier=(37.2183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21247,'amu*angstrom^2'), symmetry=1, barrier=(37.2183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374692,0.087511,-0.00010053,6.20278e-08,-1.59258e-11,-100549,30.4043], Tmin=(100,'K'), Tmax=(927.245,'K')), NASAPolynomial(coeffs=[11.9293,0.0376654,-1.98941e-05,4.05135e-09,-2.94206e-13,-102691,-24.4668], Tmin=(927.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-837.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'O=C(F)CF(476)',
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
    label = 'O=CC(F)F(990)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.595284,'amu*angstrom^2'), symmetry=1, barrier=(13.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93683,0.00748989,0.000222302,-1.39306e-06,2.75418e-09,-67109.3,9.32225], Tmin=(10,'K'), Tmax=(169.939,'K')), NASAPolynomial(coeffs=[3.97443,0.0203448,-1.24424e-05,3.60772e-09,-4.02215e-13,-67130.4,8.62375], Tmin=(169.939,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-557.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2CH(73)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-92.0165,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([43.4896,46.129,89.6186],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([420.994,524.562,552.014,640.304,746.121,978.508,1261.79,1779.04,3347.18],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95511,0.00251604,4.75146e-05,-8.97685e-08,4.98372e-11,-11065.2,8.58309], Tmin=(10,'K'), Tmax=(611.191,'K')), NASAPolynomial(coeffs=[3.6485,0.0155109,-1.13452e-05,3.84899e-09,-4.87703e-13,-11233,8.2324], Tmin=(611.191,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-92.0165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[CH]DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)CF(495)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-420.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,180],'cm^-1')),
        HinderedRotor(inertia=(0.439397,'amu*angstrom^2'), symmetry=1, barrier=(10.1026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10367,'amu*angstrom^2'), symmetry=1, barrier=(25.3756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.65,'J/mol'), sigma=(5.51737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.71 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88929,0.0316345,-2.19993e-05,7.29201e-09,-9.3545e-13,-50578.4,11.3001], Tmin=(10,'K'), Tmax=(2088.97,'K')), NASAPolynomial(coeffs=[7.11816,0.0206133,-1.06111e-05,2.54887e-09,-2.35115e-13,-50871.7,-4.12904], Tmin=(2088.97,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-420.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]OC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CF(125)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p1 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-105.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,734,1109,1255,1358,2983,3011,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.0460479,'amu*angstrom^2'), symmetry=1, barrier=(35.5232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4418.31,'J/mol'), sigma=(4.687e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8659,0.0107238,1.79401e-05,-3.81644e-08,1.95617e-11,-12729.5,8.48308], Tmin=(10,'K'), Tmax=(672.698,'K')), NASAPolynomial(coeffs=[3.4045,0.0188633,-1.22417e-05,3.67097e-09,-4.17395e-13,-12789.5,9.6187], Tmin=(672.698,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-105.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OOC=C(F)F(784)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {8,S}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-417.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.24886,'amu*angstrom^2'), symmetry=1, barrier=(5.72177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810155,'amu*angstrom^2'), symmetry=1, barrier=(18.6271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67377,0.0342444,-3.39692e-05,1.78405e-08,-3.9069e-12,-50168.9,10.7464], Tmin=(10,'K'), Tmax=(1013.01,'K')), NASAPolynomial(coeffs=[7.34757,0.0197379,-1.24891e-05,3.70435e-09,-4.1826e-13,-50913.2,-7.0251], Tmin=(1013.01,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-417.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OOCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FCOOC=C(F)F(4999)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-613.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.213475,'amu*angstrom^2'), symmetry=1, barrier=(4.90821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55396,'amu*angstrom^2'), symmetry=1, barrier=(35.7287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06234,'amu*angstrom^2'), symmetry=1, barrier=(47.4172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3151.58,'J/mol'), sigma=(5.06004,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.27 K, Pc=55.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3617,0.0642518,-8.1117e-05,5.97414e-08,-1.87046e-11,-73753.3,25.149], Tmin=(100,'K'), Tmax=(763.752,'K')), NASAPolynomial(coeffs=[7.64185,0.0313578,-1.65079e-05,3.34009e-09,-2.41023e-13,-74712.5,-3.45603], Tmin=(763.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'H2(8)',
    structure = adjacencyList("""1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-8.60349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3765.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (2.01594,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(496.376,'J/mol'), sigma=(2.8327,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212711,-2.78626e-07,3.40268e-10,-7.76035e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.07,'K')), NASAPolynomial(coeffs=[2.78817,0.000587632,1.59015e-07,-5.52748e-11,4.34317e-15,-596.149,0.11269], Tmin=(1959.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C]C(F)OOC=C(F)F-2(4623)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 C u0 p1 c0 {4,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-485.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.41161,'amu*angstrom^2'), symmetry=1, barrier=(32.4556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40911,'amu*angstrom^2'), symmetry=1, barrier=(32.3983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39839,'amu*angstrom^2'), symmetry=1, barrier=(32.1518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40578,'amu*angstrom^2'), symmetry=1, barrier=(32.3217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240842,0.0909294,-0.000132507,1.04149e-07,-3.33613e-11,-58231.3,31.8149], Tmin=(100,'K'), Tmax=(758.798,'K')), NASAPolynomial(coeffs=[11.209,0.033107,-1.81957e-05,3.70994e-09,-2.67631e-13,-59895.7,-18.0719], Tmin=(758.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CJ2_singlet-FCs)"""),
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
    label = 'FC(F)=COOC(F)F(5017)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-846.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.295909,'amu*angstrom^2'), symmetry=1, barrier=(6.80354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57203,'amu*angstrom^2'), symmetry=1, barrier=(36.144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3170.8,'J/mol'), sigma=(4.97696,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.27 K, Pc=58.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03607,0.0740732,-9.7785e-05,6.06602e-08,-8.91296e-12,-101700,27.3076], Tmin=(100,'K'), Tmax=(570.062,'K')), NASAPolynomial(coeffs=[8.90091,0.0310775,-1.67266e-05,3.37671e-09,-2.41949e-13,-102795,-7.95279], Tmin=(570.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-846.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'FC(F)=COOCC(F)F(5023)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-847.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,217.347,218.906],'cm^-1')),
        HinderedRotor(inertia=(0.138875,'amu*angstrom^2'), symmetry=1, barrier=(4.74018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23369,'amu*angstrom^2'), symmetry=1, barrier=(42.1775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2232,'amu*angstrom^2'), symmetry=1, barrier=(42.1777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292655,'amu*angstrom^2'), symmetry=1, barrier=(9.92851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.85,'J/mol'), sigma=(5.29792,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.55 K, Pc=50.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0928902,0.0984951,-0.00014399,1.18329e-07,-3.99546e-11,-101792,32.3113], Tmin=(100,'K'), Tmax=(719.868,'K')), NASAPolynomial(coeffs=[10.5512,0.0393423,-2.07151e-05,4.14939e-09,-2.96086e-13,-103325,-15.5402], Tmin=(719.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-847.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'O=CC(F)(F)OC(F)CF(10725)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
10 C u0 p0 c0 {6,D} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1191.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.278575,0.0886971,-0.000105375,6.63918e-08,-1.71618e-11,-143153,28.7868], Tmin=(100,'K'), Tmax=(928.408,'K')), NASAPolynomial(coeffs=[12.8307,0.0346142,-1.79911e-05,3.64043e-09,-2.63449e-13,-145484,-30.8371], Tmin=(928.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1191.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH)"""),
)

species(
    label = 'FC[C](F)OOC[C](F)F(10666)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
9  C u1 p0 c0 {2,S} {6,S} {8,S}
10 C u1 p0 c0 {3,S} {4,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-635.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,190,488,555,1236,1407,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.801512,0.117271,-0.000198424,1.77753e-07,-6.26021e-11,-76317.3,34.7406], Tmin=(100,'K'), Tmax=(801.56,'K')), NASAPolynomial(coeffs=[11.3952,0.0395367,-2.13854e-05,4.25026e-09,-2.98773e-13,-77730.6,-18.0229], Tmin=(801.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-635.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + radical(CsCsF1sO2s) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'FC[C](F)OO[CH]C(F)F(10662)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {13,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u1 p0 c0 {6,S} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-651.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,522,611,926,1093,1137,1374,1416,3112,395,473,707,1436,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3245.02,'J/mol'), sigma=(5.45309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.86 K, Pc=45.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.733978,0.115325,-0.000192016,1.71307e-07,-6.04344e-11,-78246.8,33.9558], Tmin=(100,'K'), Tmax=(794.031,'K')), NASAPolynomial(coeffs=[11.1598,0.0402848,-2.16886e-05,4.31216e-09,-3.03609e-13,-79658.8,-17.6793], Tmin=(794.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-651.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + radical(CsCsF1sO2s) + radical(CCsJOOC)"""),
)

species(
    label = 'F[CH]C(F)OOC[C](F)F(10717)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-631.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,487,638,688,1119,1325,1387,3149,190,488,555,1236,1407,334,575,1197,1424,3202,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.845522,0.116202,-0.000186685,1.59157e-07,-5.427e-11,-75822.1,34.4521], Tmin=(100,'K'), Tmax=(754.611,'K')), NASAPolynomial(coeffs=[13.322,0.0366787,-1.98151e-05,3.96477e-09,-2.81029e-13,-77834.3,-29.0744], Tmin=(754.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-631.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)OO[CH]C(F)F(5135)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 C u1 p0 c0 {4,S} {7,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-647.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,334,575,1197,1424,3202,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.462391,0.109642,-0.000158755,1.13763e-07,-2.82584e-11,-77765.1,32.5849], Tmin=(100,'K'), Tmax=(610.823,'K')), NASAPolynomial(coeffs=[12.8509,0.0378815,-2.04043e-05,4.09854e-09,-2.92104e-13,-79679.2,-27.4364], Tmin=(610.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-647.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCsJOOC) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'FC[CH]OOC=C(F)F(5246)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {4,S} {6,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-439.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,855.614],'cm^-1')),
        HinderedRotor(inertia=(0.310385,'amu*angstrom^2'), symmetry=1, barrier=(7.13635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77275,'amu*angstrom^2'), symmetry=1, barrier=(40.759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77348,'amu*angstrom^2'), symmetry=1, barrier=(40.7759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896879,0.0802225,-8.53507e-05,1.24645e-08,3.32882e-11,-52767.5,28.393], Tmin=(100,'K'), Tmax=(510.43,'K')), NASAPolynomial(coeffs=[8.86489,0.0378556,-1.98403e-05,3.95236e-09,-2.8059e-13,-53842.5,-7.25122], Tmin=(510.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-439.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]C(F)OOC=C(F)F(4983)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
7  C u1 p0 c0 {6,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-448.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.143754,'amu*angstrom^2'), symmetry=1, barrier=(3.30518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884035,'amu*angstrom^2'), symmetry=1, barrier=(20.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89975,'amu*angstrom^2'), symmetry=1, barrier=(43.679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90064,'amu*angstrom^2'), symmetry=1, barrier=(43.6994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3321.75,'J/mol'), sigma=(5.41148,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.85 K, Pc=47.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627947,0.0812749,-9.76684e-05,6.3474e-08,-1.70946e-11,-53876.7,30.1803], Tmin=(100,'K'), Tmax=(888.89,'K')), NASAPolynomial(coeffs=[11.1356,0.0339902,-1.78747e-05,3.62816e-09,-2.62813e-13,-55744.7,-19.2751], Tmin=(888.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-448.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CJCOOH)"""),
)

species(
    label = 'F[C]=COOC(F)CF(4677)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-385.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.884042,'amu*angstrom^2'), symmetry=1, barrier=(20.3259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72392,'amu*angstrom^2'), symmetry=1, barrier=(39.6363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.882332,'amu*angstrom^2'), symmetry=1, barrier=(20.2865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71887,'amu*angstrom^2'), symmetry=1, barrier=(39.5202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.378964,0.0873626,-0.000116929,8.58e-08,-2.59139e-11,-46263,29.6701], Tmin=(100,'K'), Tmax=(800.096,'K')), NASAPolynomial(coeffs=[10.8872,0.0348263,-1.8433e-05,3.72738e-09,-2.68549e-13,-47944.4,-18.6821], Tmin=(800.096,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'O=C[C](F)F(787)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-391.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([179,346,818,1406,1524,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.55681,'amu*angstrom^2'), symmetry=1, barrier=(31.8674,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92143,0.00544306,7.55479e-05,-1.88818e-07,1.40422e-10,-47057,10.1435], Tmin=(10,'K'), Tmax=(452.346,'K')), NASAPolynomial(coeffs=[3.74847,0.0196353,-1.35046e-05,4.31321e-09,-5.19454e-13,-47170.9,9.4087], Tmin=(452.346,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-391.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)CF(490)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-437.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([528,1116,1182,1331,1402,1494,3075,3110,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.906018,'amu*angstrom^2'), symmetry=1, barrier=(20.8311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85154,0.0137524,4.01343e-05,-8.47153e-08,4.68171e-11,-52581.1,10.8499], Tmin=(10,'K'), Tmax=(631.493,'K')), NASAPolynomial(coeffs=[3.88776,0.0255121,-1.62765e-05,4.90162e-09,-5.62996e-13,-52824.7,8.7991], Tmin=(631.493,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-437.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[O]C(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC=C(F)F(785)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-270.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.857306,'amu*angstrom^2'), symmetry=1, barrier=(19.7111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3441.78,'J/mol'), sigma=(5.35945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.60 K, Pc=50.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84,0.010514,0.000109777,-2.90322e-07,2.13215e-10,-32466.3,11.2834], Tmin=(10,'K'), Tmax=(492.518,'K')), NASAPolynomial(coeffs=[6.26831,0.0212777,-1.585e-05,5.39864e-09,-6.82851e-13,-33075.2,-2.46567], Tmin=(492.518,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-270.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2F-CHF(66)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-267.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,334,575,1197,1424,3202,700.651],'cm^-1')),
        HinderedRotor(inertia=(0.57132,'amu*angstrom^2'), symmetry=1, barrier=(13.1358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.042,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(2,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93382,0.0142501,1.77758e-06,-8.62923e-09,3.21423e-12,-32220.6,9.76326], Tmin=(10,'K'), Tmax=(1157.03,'K')), NASAPolynomial(coeffs=[6.23595,0.0134751,-6.53103e-06,1.52429e-09,-1.39119e-13,-33234.2,-3.75694], Tmin=(1157.03,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-267.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[CH]CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[CH]OOC=C(F)F(1677)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,D} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {6,D}
8  C u1 p0 c0 {3,S} {5,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-414.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.25081,'amu*angstrom^2'), symmetry=1, barrier=(5.76662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66128,'amu*angstrom^2'), symmetry=1, barrier=(38.1961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998036,0.07318,-0.000119119,1.06928e-07,-3.84225e-11,-49799.7,26.6702], Tmin=(100,'K'), Tmax=(774.496,'K')), NASAPolynomial(coeffs=[8.01183,0.0282627,-1.52892e-05,3.06084e-09,-2.16921e-13,-50625.4,-3.69169], Tmin=(774.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsF1sHO2s)"""),
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
    label = 'FC[C](F)OOC=C(F)F(4624)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {5,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-642.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.993603,'amu*angstrom^2'), symmetry=1, barrier=(22.8449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424255,'amu*angstrom^2'), symmetry=1, barrier=(9.75445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91473,'amu*angstrom^2'), symmetry=1, barrier=(44.0233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90707,'amu*angstrom^2'), symmetry=1, barrier=(43.8474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20727,0.101774,-0.000163382,1.43481e-07,-5.06709e-11,-77130.7,33.2663], Tmin=(100,'K'), Tmax=(760.997,'K')), NASAPolynomial(coeffs=[10.5581,0.0371293,-2.00767e-05,4.02472e-09,-2.8573e-13,-78535.8,-14.1971], Tmin=(760.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-642.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[CH]C(F)OOC=C(F)F(4579)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {7,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-638.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67729,'amu*angstrom^2'), symmetry=1, barrier=(38.5643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67858,'amu*angstrom^2'), symmetry=1, barrier=(38.5938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4336,'amu*angstrom^2'), symmetry=1, barrier=(9.96932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67989,'amu*angstrom^2'), symmetry=1, barrier=(38.624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.85,'J/mol'), sigma=(5.29792,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.55 K, Pc=50.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0936353,0.098574,-0.00014272,1.10721e-07,-3.48592e-11,-76642.2,32.4285], Tmin=(100,'K'), Tmax=(772.882,'K')), NASAPolynomial(coeffs=[12.244,0.0347172,-1.87792e-05,3.80642e-09,-2.73737e-13,-78549.2,-23.914], Tmin=(772.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-638.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'FCC(F)OO[C]=C(F)F(4625)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-597.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.506293,'amu*angstrom^2'), symmetry=1, barrier=(11.6407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0142669,'amu*angstrom^2'), symmetry=1, barrier=(38.3463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66733,'amu*angstrom^2'), symmetry=1, barrier=(38.3352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66797,'amu*angstrom^2'), symmetry=1, barrier=(38.35,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06976,0.0792848,-6.14431e-05,-6.57658e-08,1.07498e-10,-71748.2,31.2555], Tmin=(100,'K'), Tmax=(458.767,'K')), NASAPolynomial(coeffs=[8.62283,0.040303,-2.18533e-05,4.38823e-09,-3.1179e-13,-72724,-2.38039], Tmin=(458.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C=CJO)"""),
)

species(
    label = 'FCC(F)O[O+]=[C-]C(F)F(10608)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p1 c+1 {5,S} {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u0 p1 c-1 {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-709.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.370435,0.109104,-0.000185172,1.73442e-07,-6.37028e-11,-85142.4,43.3846], Tmin=(100,'K'), Tmax=(807.64,'K')), NASAPolynomial(coeffs=[7.53812,0.0463502,-2.48178e-05,4.92023e-09,-3.45588e-13,-85650.7,11.682], Tmin=(807.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-709.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'F[C]C(F)OOC(F)CF(10610)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {6,S} {10,S} {14,S}
10 C u0 p1 c0 {4,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-670.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,2750,2850,1437.5,1250,1305,750,350,617,898,1187,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49519,'amu*angstrom^2'), symmetry=1, barrier=(34.3773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4949,'amu*angstrom^2'), symmetry=1, barrier=(34.3706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49535,'amu*angstrom^2'), symmetry=1, barrier=(34.3811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49608,'amu*angstrom^2'), symmetry=1, barrier=(34.3979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49428,'amu*angstrom^2'), symmetry=1, barrier=(34.3563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.137276,0.0922484,-0.000106514,6.38192e-08,-1.56714e-11,-80523.7,30.4519], Tmin=(100,'K'), Tmax=(973.616,'K')), NASAPolynomial(coeffs=[13.9647,0.0354403,-1.89931e-05,3.89159e-09,-2.83622e-13,-83216.2,-35.8879], Tmin=(973.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CJ2_singlet-FCs)"""),
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
    label = 'C=COOC=C(F)F(4698)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-294.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,182,240,577,636,1210,1413,240.062],'cm^-1')),
        HinderedRotor(inertia=(0.315415,'amu*angstrom^2'), symmetry=1, barrier=(12.9001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57616,'amu*angstrom^2'), symmetry=1, barrier=(23.5684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793051,'amu*angstrom^2'), symmetry=1, barrier=(32.4436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3304.43,'J/mol'), sigma=(5.27966,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=516.14 K, Pc=50.95 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14602,0.0635932,-6.19324e-05,3.08471e-08,-6.24095e-12,-35288.2,26.8148], Tmin=(100,'K'), Tmax=(1175.03,'K')), NASAPolynomial(coeffs=[12.7329,0.0241493,-1.15798e-05,2.27892e-09,-1.62774e-13,-38011.2,-30.9544], Tmin=(1175.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    label = 'C=C(F)OOC=C(F)F(1750)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-479.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,326,540,652,719,1357,182,240,577,636,1210,1413,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.148573,'amu*angstrom^2'), symmetry=1, barrier=(3.41598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14843,'amu*angstrom^2'), symmetry=1, barrier=(3.41269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28631,'amu*angstrom^2'), symmetry=1, barrier=(52.5667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717383,0.0801711,-0.000124309,1.12571e-07,-4.15644e-11,-57554.8,29.7561], Tmin=(100,'K'), Tmax=(757.959,'K')), NASAPolynomial(coeffs=[7.15173,0.036726,-1.9553e-05,3.91505e-09,-2.78368e-13,-58257.7,2.29525], Tmin=(757.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'FC=COOC=C(F)F(2505)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-459.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,194,682,905,1196,1383,3221,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.0836088,'amu*angstrom^2'), symmetry=1, barrier=(21.7354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22505,'amu*angstrom^2'), symmetry=1, barrier=(5.17434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58333,'amu*angstrom^2'), symmetry=1, barrier=(59.3959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918625,0.0747774,-9.67782e-05,7.10186e-08,-2.18264e-11,-55130.2,29.1101], Tmin=(100,'K'), Tmax=(781.67,'K')), NASAPolynomial(coeffs=[8.97024,0.0335749,-1.77109e-05,3.5833e-09,-2.58452e-13,-56388.9,-7.75092], Tmin=(781.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'FC#COOC(F)CF(10611)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {9,T}
9  C u0 p0 c0 {3,S} {8,T}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-347.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,2175,525,239,401,1367,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47126,'amu*angstrom^2'), symmetry=1, barrier=(33.8271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47058,'amu*angstrom^2'), symmetry=1, barrier=(33.8115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47017,'amu*angstrom^2'), symmetry=1, barrier=(33.8021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4716,'amu*angstrom^2'), symmetry=1, barrier=(33.8349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632893,0.082099,-0.000106653,7.47145e-08,-2.16375e-11,-41727.6,26.87], Tmin=(100,'K'), Tmax=(830.249,'K')), NASAPolynomial(coeffs=[10.7834,0.0331955,-1.82993e-05,3.76887e-09,-2.74648e-13,-43413.1,-20.2119], Tmin=(830.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(CtCF)"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = '[C-]=[O+]OC(F)CF(10726)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p1 c+1 {3,S} {7,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p1 c-1 {4,D} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-268.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,180,180,180,385.596,782.824,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154655,'amu*angstrom^2'), symmetry=1, barrier=(3.55583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154655,'amu*angstrom^2'), symmetry=1, barrier=(3.55583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154655,'amu*angstrom^2'), symmetry=1, barrier=(3.55583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1267,0.0718778,-0.000119698,1.13279e-07,-4.23168e-11,-32250.5,23.8477], Tmin=(100,'K'), Tmax=(798.865,'K')), NASAPolynomial(coeffs=[5.67919,0.0333345,-1.77557e-05,3.53007e-09,-2.48838e-13,-32475.4,6.05217], Tmin=(798.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-461.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-124.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (0.247794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-64.5949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (82.0061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-190.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-396.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-279.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-295.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-259.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-275.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-19.0649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-28.4167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (34.8198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-445.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-190.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-164.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-109.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-83.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-78.9884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-37.8849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-302.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-178.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (50.6215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-203.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-175.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-12.2182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-125.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FCC(F)OOC=C(F)F(5256)'],
    products = ['O=C(F)CF(476)', 'O=CC(F)F(990)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(28.1932,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]CF(125)', 'OOC=C(F)F(784)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(50.6963,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'FCOOC=C(F)F(4999)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.12553e-07,'m^3/(mol*s)'), n=3.34134, Ea=(127.848,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2(8)', 'F[C]C(F)OOC=C(F)F-2(4623)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(S)(25)', 'FC(F)=COOC(F)F(5017)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(161.708,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FCC(F)OOC=C(F)F(5256)'],
    products = ['FC(F)=COOCC(F)F(5023)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FCC(F)OOC=C(F)F(5256)'],
    products = ['O=CC(F)(F)OC(F)CF(10725)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(92.6483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FC[C](F)OOC[C](F)F(10666)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC[C](F)OO[CH]C(F)F(10662)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(F)OOC[C](F)F(10717)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(F)OO[CH]C(F)F(5135)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'FC[CH]OOC=C(F)F(5246)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[CH2]C(F)OOC=C(F)F(4983)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[C]=COOC(F)CF(4677)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C[C](F)F(787)', '[O]C(F)CF(490)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(35.0126,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC=C(F)F(785)', 'CH2F-CHF(66)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CF2CH(73)', '[O]OC(F)CF(495)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2F(46)', 'F[CH]OOC=C(F)F(1677)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'FC[C](F)OOC=C(F)F(4624)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.57075e+24,'m^3/(mol*s)'), n=-6.76739, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-2C-R
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'F[CH]C(F)OOC=C(F)F(4579)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -1.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', 'FCC(F)OO[C]=C(F)F(4625)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Sp-3C=2C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Sp-3C=2C
Ea raised from -4.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['FCC(F)O[O+]=[C-]C(F)F(10608)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.18269e+17,'s^-1'), n=-1.42952, Ea=(58.7793,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5517143818398826, var=16.398565006535087, Tref=1000.0, N=6, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]C(F)OOC(F)CF(10610)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(144.923,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F2(78)', 'C=COOC=C(F)F(4698)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(6.07059,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'C=C(F)OOC=C(F)F(1750)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(52.9173,'m^3/(mol*s)'), n=1.09798, Ea=(209.515,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd_Ext-3Cd-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd_Ext-3Cd-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'FC=COOC=C(F)F(2505)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(217.288,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', 'FC#COOC(F)CF(10611)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(269.181,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CF2(43)', '[C-]=[O+]OC(F)CF(10726)'],
    products = ['FCC(F)OOC=C(F)F(5256)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.66533e+06,'cm^3/(mol*s)'), n=1.53, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for CF2
Exact match found for rate rule [CF2]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #3524',
    isomers = [
        'FCC(F)OOC=C(F)F(5256)',
    ],
    reactants = [
        ('O=C(F)CF(476)', 'O=CC(F)F(990)'),
        ('CF2CH(73)', '[O]OC(F)CF(495)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3524',
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

