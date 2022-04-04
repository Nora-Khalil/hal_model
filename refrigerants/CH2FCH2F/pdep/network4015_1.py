species(
    label = '[O]OC(F)=C(F)OOCF(12362)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {2,S} {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-534.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,548,1085,1183,1302,1466,1520,3060,3119,262,390,483,597,572,732,631,807,1275,1439,180],'cm^-1')),
        HinderedRotor(inertia=(0.260258,'amu*angstrom^2'), symmetry=1, barrier=(5.98384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261228,'amu*angstrom^2'), symmetry=1, barrier=(6.00614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261895,'amu*angstrom^2'), symmetry=1, barrier=(6.02148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260222,'amu*angstrom^2'), symmetry=1, barrier=(5.98301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101812,0.103807,-0.000203801,2.04373e-07,-7.66069e-11,-64144.6,34.5034], Tmin=(100,'K'), Tmax=(843.38,'K')), NASAPolynomial(coeffs=[4.03564,0.0440605,-2.44581e-05,4.841e-09,-3.36536e-13,-63346.8,24.8585], Tmin=(843.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-534.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = 'CHFO(47)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)C(=O)F(5430)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-582.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.519618,'amu*angstrom^2'), symmetry=1, barrier=(11.947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54034,'amu*angstrom^2'), symmetry=1, barrier=(35.4154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3788.12,'J/mol'), sigma=(5.7366,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.70 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66212,0.0558676,-8.74785e-05,7.3629e-08,-2.48671e-11,-69953.7,20.9885], Tmin=(100,'K'), Tmax=(757.894,'K')), NASAPolynomial(coeffs=[8.40295,0.0184183,-9.65364e-06,1.91161e-09,-1.34842e-13,-70921.6,-9.30864], Tmin=(757.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=C(F)OO(12063)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {3,S} {9,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {3,S} {8,D}
8 C u0 p0 c0 {2,S} {4,S} {7,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-319.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439],'cm^-1')),
        HinderedRotor(inertia=(0.225342,'amu*angstrom^2'), symmetry=1, barrier=(5.18105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226784,'amu*angstrom^2'), symmetry=1, barrier=(5.21421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225602,'amu*angstrom^2'), symmetry=1, barrier=(5.18704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01189,0.0821807,-0.000172845,1.75767e-07,-6.54215e-11,-38387.2,27.3796], Tmin=(100,'K'), Tmax=(862.081,'K')), NASAPolynomial(coeffs=[3.35234,0.0322777,-1.80806e-05,3.54828e-09,-2.43747e-13,-37339.9,24.8503], Tmin=(862.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=C(F)OOF(12064)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {9,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,S} {9,S}
6 O u0 p2 c0 {3,S} {4,S}
7 O u1 p2 c0 {5,S}
8 C u0 p0 c0 {2,S} {4,S} {9,D}
9 C u0 p0 c0 {1,S} {5,S} {8,D}
"""),
    E0 = (-224.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,277,555,632,262,390,483,597,572,732,631,807,1275,1439,180,2414.25],'cm^-1')),
        HinderedRotor(inertia=(0.199336,'amu*angstrom^2'), symmetry=1, barrier=(4.58314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199175,'amu*angstrom^2'), symmetry=1, barrier=(4.57942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198769,'amu*angstrom^2'), symmetry=1, barrier=(4.57008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648634,0.0918949,-0.000195796,1.98318e-07,-7.35179e-11,-26913.8,30.0644], Tmin=(100,'K'), Tmax=(861.101,'K')), NASAPolynomial(coeffs=[4.02078,0.033372,-1.91937e-05,3.79163e-09,-2.61054e-13,-25905.5,23.5265], Tmin=(861.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(OCF)C(=O)F(12372)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-984.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.237902,0.104389,-0.000183673,1.65902e-07,-5.79482e-11,-118247,31.4905], Tmin=(100,'K'), Tmax=(827.034,'K')), NASAPolynomial(coeffs=[10.6749,0.0321585,-1.73914e-05,3.42755e-09,-2.38648e-13,-119387,-15.0638], Tmin=(827.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-984.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHHO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = 'O=C(F)[C](F)OOCF(12373)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {9,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-809.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,280,501,1494,1531,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.9671,'amu*angstrom^2'), symmetry=1, barrier=(45.2274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30019,'amu*angstrom^2'), symmetry=1, barrier=(69.9538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97556,'amu*angstrom^2'), symmetry=1, barrier=(45.422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41397,0.0672354,-6.82973e-05,5.22959e-09,3.13971e-11,-97276.3,25.3162], Tmin=(100,'K'), Tmax=(500.139,'K')), NASAPolynomial(coeffs=[7.54364,0.0344291,-1.85446e-05,3.74465e-09,-2.68382e-13,-98092.3,-2.0366], Tmin=(500.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-809.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHHO) + group(COCsFO) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]CF(180)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-220.7,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(49.009,'amu')),
        NonlinearRotor(inertia=([8.89593,46.7012,52.3875],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([558.515,760.247,1040.95,1161.07,1162.34,1315.05,1378.4,2858.27,2872.23],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (49.0244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3043.69,'J/mol'), sigma=(5.07813,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.42 K, Pc=52.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08121,-0.00607739,5.67497e-05,-8.01653e-08,3.67029e-11,-26544.2,7.61089], Tmin=(10,'K'), Tmax=(682.778,'K')), NASAPolynomial(coeffs=[1.79187,0.0157666,-9.76404e-06,2.86626e-09,-3.21991e-13,-26428.1,16.3427], Tmin=(682.778,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-220.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=C(F)OOO1(12067)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u0 p0 c0 {1,S} {4,S} {6,D}
"""),
    E0 = (-242.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11876,0.00894386,3.87231e-05,-5.40048e-08,1.99787e-11,-29142.3,19.1532], Tmin=(100,'K'), Tmax=(1010.02,'K')), NASAPolynomial(coeffs=[8.91172,0.0124165,-5.66298e-06,1.18562e-09,-9.10784e-14,-31659.8,-15.522], Tmin=(1010.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CdCFO) + group(CdCFO) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + ring(Cyclopentane)"""),
)

species(
    label = 'FCOOC1(F)OO[C]1F(12374)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {5,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-595.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281081,0.0898456,-0.00013467,1.1096e-07,-3.73562e-11,-71553.3,31.0932], Tmin=(100,'K'), Tmax=(722.76,'K')), NASAPolynomial(coeffs=[10.4028,0.0338299,-1.84195e-05,3.73377e-09,-2.68053e-13,-73016.5,-14.4522], Tmin=(722.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + group(CsFHHO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCOO[C](F)C1(F)OO1(12375)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-606.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.256884,0.090357,-0.000133303,1.07511e-07,-3.54046e-11,-72839.6,30.5083], Tmin=(100,'K'), Tmax=(738.626,'K')), NASAPolynomial(coeffs=[10.71,0.0337409,-1.83123e-05,3.70866e-09,-2.66278e-13,-74383.6,-16.7536], Tmin=(738.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHHO) + ring(O2s-O2s-Cs(C-F)) + radical(CsCsF1sO2s)"""),
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
    label = '[CH2]OOC(F)=C(F)O[O](12158)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {2,S} {5,S} {7,D}
9  C u1 p0 c0 {4,S} {10,S} {11,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-126.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.202864,'amu*angstrom^2'), symmetry=1, barrier=(4.66425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201722,'amu*angstrom^2'), symmetry=1, barrier=(4.63799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201901,'amu*angstrom^2'), symmetry=1, barrier=(4.64211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202176,'amu*angstrom^2'), symmetry=1, barrier=(4.64842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440055,0.0954177,-0.000189341,1.90221e-07,-7.10866e-11,-15092.1,31.3824], Tmin=(100,'K'), Tmax=(849.516,'K')), NASAPolynomial(coeffs=[3.78324,0.0402154,-2.21932e-05,4.37083e-09,-3.02539e-13,-14236.3,24.1794], Tmin=(849.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsJOOC)"""),
)

species(
    label = '[O]O[C]=C(F)OOCF(12376)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-107.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,548,1085,1183,1302,1466,1520,3060,3119,293,496,537,1218,1685,370,180,1350.29],'cm^-1')),
        HinderedRotor(inertia=(0.251283,'amu*angstrom^2'), symmetry=1, barrier=(5.77749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251152,'amu*angstrom^2'), symmetry=1, barrier=(5.77447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250958,'amu*angstrom^2'), symmetry=1, barrier=(5.77003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250603,'amu*angstrom^2'), symmetry=1, barrier=(5.76185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715978,0.0901732,-0.000182713,1.88598e-07,-7.18525e-11,-12831.6,34.7088], Tmin=(100,'K'), Tmax=(848.828,'K')), NASAPolynomial(coeffs=[1.82377,0.0425104,-2.34842e-05,4.63335e-09,-3.21196e-13,-11490.6,38.5525], Tmin=(848.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]OC(F)=[C]OOCF(6265)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-107.527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,548,1085,1183,1302,1466,1520,3060,3119,293,496,537,1218,1685,370,180,1350.29],'cm^-1')),
        HinderedRotor(inertia=(0.251283,'amu*angstrom^2'), symmetry=1, barrier=(5.77749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251152,'amu*angstrom^2'), symmetry=1, barrier=(5.77447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250958,'amu*angstrom^2'), symmetry=1, barrier=(5.77003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250603,'amu*angstrom^2'), symmetry=1, barrier=(5.76185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715978,0.0901732,-0.000182713,1.88598e-07,-7.18525e-11,-12831.6,34.7088], Tmin=(100,'K'), Tmax=(848.828,'K')), NASAPolynomial(coeffs=[1.82377,0.0425104,-2.34842e-05,4.63335e-09,-3.21196e-13,-11490.6,38.5525], Tmin=(848.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]O[C](F)C(=O)F(11195)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
"""),
    E0 = (-443.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,280,501,1494,1531,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(4.168,'amu*angstrom^2'), symmetry=1, barrier=(95.8305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.16864,'amu*angstrom^2'), symmetry=1, barrier=(95.8453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93924,0.0519598,-9.19914e-05,8.61317e-08,-3.09905e-11,-53220.8,20.2085], Tmin=(100,'K'), Tmax=(837.266,'K')), NASAPolynomial(coeffs=[5.80391,0.0200067,-1.05786e-05,2.06413e-09,-1.42884e-13,-53395.1,5.0737], Tmin=(837.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
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
    label = '[O]OC(F)=C(F)O[O](12072)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u1 p2 c0 {3,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {2,S} {3,S} {8,D}
8 C u0 p0 c0 {1,S} {4,S} {7,D}
"""),
    E0 = (-167.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,262,390,483,597,572,732,631,807,1275,1439],'cm^-1')),
        HinderedRotor(inertia=(0.108211,'amu*angstrom^2'), symmetry=1, barrier=(2.48799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108864,'amu*angstrom^2'), symmetry=1, barrier=(2.50299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27926,0.0785338,-0.000177258,1.84473e-07,-6.87794e-11,-20116.9,26.4865], Tmin=(100,'K'), Tmax=(877.777,'K')), NASAPolynomial(coeffs=[1.70078,0.0307284,-1.71552e-05,3.32335e-09,-2.24948e-13,-18423.2,34.577], Tmin=(877.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OCF(414)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-199.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,548,1085,1183,1302,1466,1520,3060,3119],'cm^-1')),
        HinderedRotor(inertia=(0.499266,'amu*angstrom^2'), symmetry=1, barrier=(11.4791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0238,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3197.75,'J/mol'), sigma=(5.27922,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=499.48 K, Pc=49.32 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91825,0.0104144,6.14512e-06,-1.3041e-08,5.04468e-12,-24048.7,9.83513], Tmin=(10,'K'), Tmax=(997.565,'K')), NASAPolynomial(coeffs=[4.63802,0.0134503,-7.32448e-06,1.91158e-09,-1.93959e-13,-24487,4.88755], Tmin=(997.565,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-199.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]OCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)=[C]F(7257)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u1 p0 c0 {2,S} {5,D}
"""),
    E0 = (0.671019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,293,496,537,1218,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.147305,'amu*angstrom^2'), symmetry=1, barrier=(3.38683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94088,0.054583,-0.000112996,1.10271e-07,-3.96438e-11,145.961,18.1124], Tmin=(100,'K'), Tmax=(865.502,'K')), NASAPolynomial(coeffs=[5.76969,0.0156862,-8.83981e-06,1.73988e-09,-1.19298e-13,277.292,4.78112], Tmin=(865.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.671019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C]=C(F)OOCF(10369)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-339.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,293,496,537,1218,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.212349,'amu*angstrom^2'), symmetry=1, barrier=(4.88231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213151,'amu*angstrom^2'), symmetry=1, barrier=(4.90077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0554862,'amu*angstrom^2'), symmetry=1, barrier=(51.6616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08344,0.0745967,-0.000133358,1.30166e-07,-4.91309e-11,-40688.9,26.8017], Tmin=(100,'K'), Tmax=(811.207,'K')), NASAPolynomial(coeffs=[5.09713,0.0336661,-1.85848e-05,3.71893e-09,-2.62131e-13,-40644.5,12.5648], Tmin=(811.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    label = '[O]OC(F)=C(F)OO[CH]F(12070)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {5,S} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-335.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.258054,'amu*angstrom^2'), symmetry=1, barrier=(5.93316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257854,'amu*angstrom^2'), symmetry=1, barrier=(5.92858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258064,'amu*angstrom^2'), symmetry=1, barrier=(5.93341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257685,'amu*angstrom^2'), symmetry=1, barrier=(5.92469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0190581,0.109512,-0.000228652,2.31393e-07,-8.61047e-11,-40201.3,35.1758], Tmin=(100,'K'), Tmax=(856.969,'K')), NASAPolynomial(coeffs=[3.85172,0.041976,-2.3852e-05,4.71174e-09,-3.25212e-13,-39048.3,27.6973], Tmin=(856.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsF1sHO2s)"""),
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
    label = '[O]OC#COOCF(6266)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {4,S} {7,T}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (102.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,548,1085,1183,1302,1466,1520,3060,3119,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(1.40982,'amu*angstrom^2'), symmetry=1, barrier=(32.4146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41034,'amu*angstrom^2'), symmetry=1, barrier=(32.4265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41,'amu*angstrom^2'), symmetry=1, barrier=(32.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4099,'amu*angstrom^2'), symmetry=1, barrier=(32.4165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862625,0.075594,-0.000121268,1.05201e-07,-3.65748e-11,12382.7,28.2577], Tmin=(100,'K'), Tmax=(766.87,'K')), NASAPolynomial(coeffs=[9.30125,0.0261766,-1.40425e-05,2.80055e-09,-1.98109e-13,11247.2,-9.17826], Tmin=(766.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(CsFHHO) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)=C(F)OO[CH]F(12377)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {5,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-487.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,262,390,483,597,572,732,631,807,1275,1439,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.28187,'amu*angstrom^2'), symmetry=1, barrier=(6.48076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280784,'amu*angstrom^2'), symmetry=1, barrier=(6.45577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282103,'amu*angstrom^2'), symmetry=1, barrier=(6.48609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280833,'amu*angstrom^2'), symmetry=1, barrier=(6.4569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24776,'amu*angstrom^2'), symmetry=1, barrier=(51.6804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.280259,0.113093,-0.000224061,2.22577e-07,-8.27811e-11,-58471.9,35.3534], Tmin=(100,'K'), Tmax=(841.091,'K')), NASAPolynomial(coeffs=[5.44233,0.0436314,-2.48396e-05,4.95155e-09,-3.45257e-13,-57940.2,17.6188], Tmin=(841.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sHO2s)"""),
)

species(
    label = '[CH2]OOC(F)=C(F)OOF(12378)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {3,S} {6,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {5,S} {11,S} {12,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-183.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,262,390,483,597,572,732,631,807,1275,1439,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.11223,'amu*angstrom^2'), symmetry=1, barrier=(48.5643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412881,'amu*angstrom^2'), symmetry=1, barrier=(9.49294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.181295,0.108679,-0.000207604,2.03888e-07,-7.58705e-11,-21889.4,34.2335], Tmin=(100,'K'), Tmax=(830.651,'K')), NASAPolynomial(coeffs=[6.01434,0.0430141,-2.43228e-05,4.86091e-09,-3.40472e-13,-21682.5,12.9333], Tmin=(830.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cs-OsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsJOOC)"""),
)

species(
    label = 'FCOO[C]=C(F)OOF(12379)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {3,S} {5,S}
8  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-164.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,548,1085,1183,1302,1466,1520,3060,3119,293,496,537,1218,1685,370],'cm^-1')),
        HinderedRotor(inertia=(2.04724,'amu*angstrom^2'), symmetry=1, barrier=(47.0701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04196,'amu*angstrom^2'), symmetry=1, barrier=(46.9487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340018,'amu*angstrom^2'), symmetry=1, barrier=(7.81767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337123,'amu*angstrom^2'), symmetry=1, barrier=(7.75111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0933241,0.103452,-0.000201045,2.02368e-07,-7.66874e-11,-19628.7,37.5645], Tmin=(100,'K'), Tmax=(831.251,'K')), NASAPolynomial(coeffs=[4.0573,0.0453046,-2.56111e-05,5.12279e-09,-3.59074e-13,-18937.9,27.293], Tmin=(831.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'FCOOC(F)=[C]OOF(12380)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {3,S} {6,S}
8  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-164.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,548,1085,1183,1302,1466,1520,3060,3119,293,496,537,1218,1685,370],'cm^-1')),
        HinderedRotor(inertia=(2.04724,'amu*angstrom^2'), symmetry=1, barrier=(47.0701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04196,'amu*angstrom^2'), symmetry=1, barrier=(46.9487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340018,'amu*angstrom^2'), symmetry=1, barrier=(7.81767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337123,'amu*angstrom^2'), symmetry=1, barrier=(7.75111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0933241,0.103452,-0.000201045,2.02368e-07,-7.66874e-11,-19628.7,37.5645], Tmin=(100,'K'), Tmax=(831.251,'K')), NASAPolynomial(coeffs=[4.0573,0.0453046,-2.56111e-05,5.12279e-09,-3.59074e-13,-18937.9,27.293], Tmin=(831.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-333.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (9.81288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (392.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-276.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-344.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-273.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-218.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-298.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (136.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (155.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (155.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-344.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-20.8494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-9.62431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-158.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (66.1644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (282.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-235.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (73.3021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (87.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (110.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)=C(F)OOCF(12362)'],
    products = ['CHFO(47)', '[O]OC(F)C(=O)F(5430)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(10.9484,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[O]OC(F)=C(F)OO(12063)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.67448e-06,'m^3/(mol*s)'), n=3.15288, Ea=(1.34274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node OH_2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', '[O]OC(F)=C(F)OOF(12064)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(8.52905,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)=C(F)OOCF(12362)'],
    products = ['[O]OC(F)(OCF)C(=O)F(12372)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(68.3004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[C](F)OOCF(12373)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(32.1045,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 32.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)=C(F)OOCF(12362)'],
    products = ['[O]CF(180)', 'FC1=C(F)OOO1(12067)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(71.0072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 67.6 to 71.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(F)=C(F)OOCF(12362)'],
    products = ['FCOOC1(F)OO[C]1F(12374)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)=C(F)OOCF(12362)'],
    products = ['FCOO[C](F)C1(F)OO1(12375)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[CH2]OOC(F)=C(F)O[O](12158)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]O[C]=C(F)OOCF(12376)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]OC(F)=[C]OOCF(6265)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]CF(180)', '[O]O[C](F)C(=O)F(11195)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(129.417,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Ea raised from 127.2 to 129.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2F(46)', '[O]OC(F)=C(F)O[O](12072)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.0381e+07,'m^3/(mol*s)'), n=-0.198371, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.054859309462099146, var=1.178102591472186, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OCF(414)', '[O]OC(F)=[C]F(7257)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', 'F[C]=C(F)OOCF(10369)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(F)=C(F)OO[CH]F(12070)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(78)', '[O]OC#COOCF(6266)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OOC(F)=C(F)OO[CH]F(12377)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2000,'s^-1'), n=1.9, Ea=(62.3416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H_OOCs4;C_rad_out_single;XH_out] for rate rule [R7H_OOCs4;C_rad_out_1H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OOC(F)=C(F)OOF(12378)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(66.7707,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FCOO[C]=C(F)OOF(12379)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.3079,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FCOOC(F)=[C]OOF(12380)'],
    products = ['[O]OC(F)=C(F)OOCF(12362)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(84.7403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #4015',
    isomers = [
        '[O]OC(F)=C(F)OOCF(12362)',
    ],
    reactants = [
        ('CHFO(47)', '[O]OC(F)C(=O)F(5430)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4015',
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

