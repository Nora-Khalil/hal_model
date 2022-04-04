species(
    label = 'FCC(=C(F)F)C(F)(CF)OF(11293)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {11,S} {15,S} {16,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {4,S} {5,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1081.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20949,'amu*angstrom^2'), symmetry=1, barrier=(27.8085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2087,'amu*angstrom^2'), symmetry=1, barrier=(27.7905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20817,'amu*angstrom^2'), symmetry=1, barrier=(27.7782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20827,'amu*angstrom^2'), symmetry=1, barrier=(27.7806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27707,0.126457,-0.000189029,1.50249e-07,-4.81006e-11,-129835,34.8691], Tmin=(100,'K'), Tmax=(762.728,'K')), NASAPolynomial(coeffs=[14.8933,0.0416532,-2.22516e-05,4.47519e-09,-3.19817e-13,-132301,-38.7637], Tmin=(762.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1081.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF)"""),
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
    label = 'FC=C=C(F)F(5101)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CF(126)',
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
    label = 'FCC(OF)=C(F)F(5237)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-598.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,1498.71],'cm^-1')),
        HinderedRotor(inertia=(0.301279,'amu*angstrom^2'), symmetry=1, barrier=(6.92699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0801501,'amu*angstrom^2'), symmetry=1, barrier=(41.2917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24256,0.0803686,-0.000294279,6.75169e-07,-5.89056e-10,-71970.2,12.281], Tmin=(10,'K'), Tmax=(349.554,'K')), NASAPolynomial(coeffs=[5.37599,0.0369442,-2.63563e-05,8.59739e-09,-1.04869e-12,-72003.2,5.89225], Tmin=(349.554,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-598.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""FCC(OF)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FCOF(367)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-309.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,548,1085,1183,1302,1466,1520,3060,3119,393.083],'cm^-1')),
        HinderedRotor(inertia=(0.39906,'amu*angstrom^2'), symmetry=1, barrier=(43.7624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0228,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98149,0.000984005,5.7382e-05,-1.00169e-07,5.47413e-11,-37258.4,9.41718], Tmin=(10,'K'), Tmax=(553.18,'K')), NASAPolynomial(coeffs=[1.76479,0.0223639,-1.51021e-05,4.67299e-09,-5.4307e-13,-37095.1,18.059], Tmin=(553.18,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-309.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""FCOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C(CF)=C(F)F(10927)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {3,S} {6,D}
8  C u0 p1 c0 {4,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-438.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,315,622,1128,231.143],'cm^-1')),
        HinderedRotor(inertia=(0.00665875,'amu*angstrom^2'), symmetry=1, barrier=(12.1547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837676,'amu*angstrom^2'), symmetry=1, barrier=(31.756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45201,0.0618751,-7.93874e-05,5.58699e-08,-1.63123e-11,-52632.6,22.2372], Tmin=(100,'K'), Tmax=(823.188,'K')), NASAPolynomial(coeffs=[8.84755,0.0259383,-1.39029e-05,2.83567e-09,-2.05605e-13,-53850.1,-12.0029], Tmin=(823.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(CJ2_singlet-FC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212708,-2.78619e-07,3.40262e-10,-7.7602e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.09,'K')), NASAPolynomial(coeffs=[2.78813,0.000587684,1.5899e-07,-5.52698e-11,4.34281e-15,-596.125,0.112931], Tmin=(1959.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C]C(F)(OF)C(CF)=C(F)F(11384)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  F u0 p3 c0 {12,S}
7  O u0 p2 c0 {5,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {10,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 C u0 p1 c0 {6,S} {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-729.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,617,898,1187,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02301,'amu*angstrom^2'), symmetry=1, barrier=(23.5209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2724,'amu*angstrom^2'), symmetry=1, barrier=(29.2549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02363,'amu*angstrom^2'), symmetry=1, barrier=(23.5352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02283,'amu*angstrom^2'), symmetry=1, barrier=(23.5168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (192.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35885,0.129428,-0.000220259,1.92529e-07,-6.60789e-11,-87519.7,36.0826], Tmin=(100,'K'), Tmax=(793.731,'K')), NASAPolynomial(coeffs=[14.5632,0.0363476,-2.00852e-05,4.01645e-09,-2.83071e-13,-89642.7,-34.5056], Tmin=(793.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + group(CJ2_singlet-FCs)"""),
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
    label = 'FCC(=C(F)F)C(F)OF(11385)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-852.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,236,527,855,1015,1182,1348,3236,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0241548,'amu*angstrom^2'), symmetry=1, barrier=(14.2089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.326597,'amu*angstrom^2'), symmetry=1, barrier=(7.50911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74563,'amu*angstrom^2'), symmetry=1, barrier=(40.1355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0168041,0.0964715,-0.000143552,1.15525e-07,-3.76777e-11,-102369,28.9069], Tmin=(100,'K'), Tmax=(747.579,'K')), NASAPolynomial(coeffs=[11.6697,0.0339407,-1.80836e-05,3.63557e-09,-2.59814e-13,-104116,-24.0739], Tmin=(747.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-852.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF)"""),
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
    label = 'FCC(=C(F)F)C(F)(F)OF(10923)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {10,S} {12,S} {13,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1086.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.133994,'amu*angstrom^2'), symmetry=1, barrier=(3.08079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134366,'amu*angstrom^2'), symmetry=1, barrier=(3.08935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26655,'amu*angstrom^2'), symmetry=1, barrier=(29.1205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3182.86,'J/mol'), sigma=(5.20665,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=497.15 K, Pc=51.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.657965,0.111888,-0.000184223,1.58552e-07,-5.41756e-11,-130576,31.4027], Tmin=(100,'K'), Tmax=(769.13,'K')), NASAPolynomial(coeffs=[13.2222,0.0332159,-1.81442e-05,3.6353e-09,-2.57372e-13,-132519,-30.6708], Tmin=(769.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1086.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF)"""),
)

species(
    label = 'F[C]C(=C(F)F)C(F)(CF)OF(11386)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  F u0 p3 c0 {12,S}
7  O u0 p2 c0 {5,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {11,D} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 C u0 p1 c0 {6,S} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-742.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,182,240,577,636,1210,1413,315,622,1128,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33146,'amu*angstrom^2'), symmetry=1, barrier=(30.6129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33165,'amu*angstrom^2'), symmetry=1, barrier=(30.6172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33166,'amu*angstrom^2'), symmetry=1, barrier=(30.6176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33144,'amu*angstrom^2'), symmetry=1, barrier=(30.6124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (192.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07803,0.121776,-0.000184074,1.42258e-07,-4.38316e-11,-89138.7,34.1505], Tmin=(100,'K'), Tmax=(793.968,'K')), NASAPolynomial(coeffs=[16.2248,0.0346024,-1.93787e-05,3.96613e-09,-2.86243e-13,-91886.2,-45.3333], Tmin=(793.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsCs) + group(CdCFF) + group(CJ2_singlet-FC)"""),
)

species(
    label = 'FCC(F)(C=C(F)F)OF(4577)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-866.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37868,'amu*angstrom^2'), symmetry=1, barrier=(31.6986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37728,'amu*angstrom^2'), symmetry=1, barrier=(31.6664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37342,'amu*angstrom^2'), symmetry=1, barrier=(31.5776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0258184,0.0946838,-0.000127922,9.01424e-08,-2.55964e-11,-104073,27.7929], Tmin=(100,'K'), Tmax=(855.899,'K')), NASAPolynomial(coeffs=[13.6323,0.0310939,-1.64761e-05,3.33553e-09,-2.40611e-13,-106402,-35.7331], Tmin=(855.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-866.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF)"""),
)

species(
    label = 'FCC(F)(OF)C(F)=C(F)F(11387)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {8,S} {11,D}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1025.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,323,467,575,827,1418,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35942,'amu*angstrom^2'), symmetry=1, barrier=(31.2558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35748,'amu*angstrom^2'), symmetry=1, barrier=(31.2113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35972,'amu*angstrom^2'), symmetry=1, barrier=(31.2626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.803174,0.1144,-0.000181654,1.48005e-07,-4.77997e-11,-123145,30.814], Tmin=(100,'K'), Tmax=(760.388,'K')), NASAPolynomial(coeffs=[15.1014,0.03074,-1.66297e-05,3.32925e-09,-2.36176e-13,-125563,-41.5606], Tmin=(760.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1025.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F))"""),
)

species(
    label = 'FCC(=C(F)F)C(F)(F)COF(11388)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {11,S} {15,S} {16,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {4,S} {5,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1076.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,274,345,380,539,705,1166,1213,2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.120335,'amu*angstrom^2'), symmetry=1, barrier=(2.76673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120783,'amu*angstrom^2'), symmetry=1, barrier=(2.77704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121154,'amu*angstrom^2'), symmetry=1, barrier=(2.78557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873955,'amu*angstrom^2'), symmetry=1, barrier=(20.094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4202,0.13136,-0.000217328,1.91415e-07,-6.66435e-11,-129326,36.4422], Tmin=(100,'K'), Tmax=(794.806,'K')), NASAPolynomial(coeffs=[12.9732,0.043734,-2.32899e-05,4.61565e-09,-3.24432e-13,-131134,-26.6742], Tmin=(794.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1076.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFF) + group(Cs-CsOsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF)"""),
)

species(
    label = 'OF(482)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C(F)C(CF)=C(F)F(9874)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {4,S} {5,S} {7,D}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-873.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,280,518,736,852,873,182,240,577,636,1210,1413,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0272862,'amu*angstrom^2'), symmetry=1, barrier=(3.83818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807993,'amu*angstrom^2'), symmetry=1, barrier=(18.5774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273798,0.0888236,-0.000116838,8.17774e-08,-2.32417e-11,-104984,26.8436], Tmin=(100,'K'), Tmax=(852.802,'K')), NASAPolynomial(coeffs=[12.439,0.0317614,-1.64666e-05,3.31028e-09,-2.38025e-13,-107059,-29.9089], Tmin=(852.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-873.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FOF(488)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {1,S} {2,S}
"""),
    E0 = (16.0257,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(53.9917,'amu')),
        NonlinearRotor(inertia=([8.34135,45.8552,54.1965],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([486.983,915.713,1034.58],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (53.9962,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03449,-0.00336126,4.20188e-05,-7.82672e-08,4.57369e-11,1928.28,6.37844], Tmin=(10,'K'), Tmax=(574.663,'K')), NASAPolynomial(coeffs=[3.91346,0.00598694,-4.58407e-06,1.5533e-09,-1.93094e-13,1801.75,5.67331], Tmin=(574.663,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(16.0257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""FOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(F)C(CF)=C(F)F(11389)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {4,S} {6,D}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-707.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,280,518,736,852,873,182,240,577,636,1210,1413,2950,3100,1380,975,1025,1650,313.353,315.725],'cm^-1')),
        HinderedRotor(inertia=(0.109726,'amu*angstrom^2'), symmetry=1, barrier=(7.83488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306599,'amu*angstrom^2'), symmetry=1, barrier=(22.4195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.779928,0.074984,-8.14427e-05,4.65316e-08,-1.08884e-11,-84929.9,24.2479], Tmin=(100,'K'), Tmax=(1020.99,'K')), NASAPolynomial(coeffs=[12.3824,0.0295276,-1.46586e-05,2.92339e-09,-2.10278e-13,-87299,-31.9679], Tmin=(1020.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-707.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(CdCFF) + group(Cds-CdsHH)"""),
)

species(
    label = 'FCC(OF)=C(CF)C(F)(F)F(11390)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {3,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
10 C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
11 C u0 p0 c0 {8,S} {9,S} {12,D}
12 C u0 p0 c0 {7,S} {10,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
"""),
    E0 = (-1099.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11034,0.123732,-0.000199087,1.75193e-07,-6.1624e-11,-132081,35.7473], Tmin=(100,'K'), Tmax=(779.488,'K')), NASAPolynomial(coeffs=[11.7305,0.0452846,-2.39702e-05,4.76091e-09,-3.35919e-13,-133701,-20.5571], Tmin=(779.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1099.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(CsCdFFF) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs)"""),
)

species(
    label = 'FCC(F)=C(CF)C(F)(F)OF(11391)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
11 C u0 p0 c0 {8,S} {9,S} {12,D}
12 C u0 p0 c0 {5,S} {10,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
"""),
    E0 = (-1103.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39761,0.129997,-0.000210763,1.82778e-07,-6.31257e-11,-132581,36.0915], Tmin=(100,'K'), Tmax=(778.242,'K')), NASAPolynomial(coeffs=[13.4973,0.0429762,-2.28683e-05,4.54457e-09,-3.20529e-13,-134583,-29.9974], Tmin=(778.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1103.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CsCFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsCs) + group(CdCsCdF)"""),
)

species(
    label = 'FCC(=C(F)OF)C(F)(F)CF(11392)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {11,S} {15,S} {16,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {5,S} {7,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1063.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,274,345,380,539,705,1166,1213,528,1116,1182,1331,1402,1494,3075,3110,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,326,540,652,719,1357,180,180,180,1347.2],'cm^-1')),
        HinderedRotor(inertia=(0.340613,'amu*angstrom^2'), symmetry=1, barrier=(7.83135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340619,'amu*angstrom^2'), symmetry=1, barrier=(7.83151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340622,'amu*angstrom^2'), symmetry=1, barrier=(7.83158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69922,'amu*angstrom^2'), symmetry=1, barrier=(39.0683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34494,0.131334,-0.000222544,2.01874e-07,-7.19181e-11,-127766,36.5755], Tmin=(100,'K'), Tmax=(804.568,'K')), NASAPolynomial(coeffs=[11.2056,0.0472647,-2.54034e-05,5.04247e-09,-3.54282e-13,-129084,-16.8853], Tmin=(804.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1063.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFO)"""),
)

species(
    label = 'C=C(C(F)(F)F)C(F)(CF)OF(11393)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {11,D} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1155.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08142,0.120729,-0.00016494,1.1714e-07,-3.34137e-11,-138780,33.0354], Tmin=(100,'K'), Tmax=(853.482,'K')), NASAPolynomial(coeffs=[16.4886,0.038379,-2.0201e-05,4.07562e-09,-2.93321e-13,-141779,-48.9453], Tmin=(853.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1155.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = 'F[CH]C([C](F)F)C(F)(CF)OF(7432)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {9,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
10 C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
11 C u1 p0 c0 {3,S} {8,S} {16,S}
12 C u1 p0 c0 {4,S} {5,S} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-840.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,333,384,448,608,1254,1480,528,1116,1182,1331,1402,1494,3075,3110,334,575,1197,1424,3202,190,488,555,1236,1407,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39682,0.131571,-0.000197318,1.47547e-07,-4.04546e-11,-100893,38.0668], Tmin=(100,'K'), Tmax=(632.152,'K')), NASAPolynomial(coeffs=[15.5228,0.0411359,-2.218e-05,4.45e-09,-3.16593e-13,-103364,-38.428], Tmin=(632.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-840.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(F)(OF)C(CF)[C](F)F(11394)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 C u1 p0 c0 {3,S} {9,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-833.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,316,385,515,654,689,1295,528,1116,1182,1331,1402,1494,3075,3110,190,488,555,1236,1407,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09289,0.149056,-0.000261497,2.33302e-07,-8.05687e-11,-100066,39.4959], Tmin=(100,'K'), Tmax=(825.087,'K')), NASAPolynomial(coeffs=[14.8157,0.0425003,-2.30872e-05,4.55652e-09,-3.1744e-13,-102019,-33.7564], Tmin=(825.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-833.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Csj(Cs-F1sO2sCs)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)(OF)[C](CF)C(F)F(11395)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {13,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 C u1 p0 c0 {5,S} {8,S} {16,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-880.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,316,385,515,654,689,1295,551,1088,1226,1380,1420,1481,3057,3119,522,611,926,1093,1137,1374,1416,3112,360,370,350,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85205,0.141158,-0.000235215,2.04903e-07,-7.04848e-11,-105654,37.3675], Tmin=(100,'K'), Tmax=(788.763,'K')), NASAPolynomial(coeffs=[14.9648,0.0427674,-2.31746e-05,4.61542e-09,-3.25136e-13,-107899,-37.1886], Tmin=(788.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-880.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCJ(C)CO) + radical(Csj(Cs-F1sO2sCs)(F1s)(H))"""),
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
    label = 'FCC(OF)=C(CF)[C](F)F(11396)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {1,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {11,S}
10 C u0 p0 c0 {6,S} {8,S} {9,D}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-679.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,166,258,854,1008,1066,1142,1198,1304,1253,1397,1428,1520,3056,3148,3059,3251,325,375,415,465,420,450,1700,1750,161,297,490,584,780,1358,199.698,220.886,1388.93],'cm^-1')),
        HinderedRotor(inertia=(0.145108,'amu*angstrom^2'), symmetry=1, barrier=(6.14318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00446397,'amu*angstrom^2'), symmetry=1, barrier=(6.1354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983002,'amu*angstrom^2'), symmetry=1, barrier=(35.6255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15892,'amu*angstrom^2'), symmetry=1, barrier=(35.6129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445229,0.0947672,-8.6767e-05,-4.07408e-08,9.38918e-11,-81583.7,32.8195], Tmin=(100,'K'), Tmax=(475.596,'K')), NASAPolynomial(coeffs=[10.015,0.0445385,-2.37817e-05,4.74254e-09,-3.35418e-13,-82836.1,-9.83485], Tmin=(475.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(CsCFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(OF)C(CF)=C(F)F(11397)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {7,S} {14,S} {15,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-695.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.857292,'amu*angstrom^2'), symmetry=1, barrier=(19.7108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858608,'amu*angstrom^2'), symmetry=1, barrier=(19.7411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21331,'amu*angstrom^2'), symmetry=1, barrier=(27.8964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21735,'amu*angstrom^2'), symmetry=1, barrier=(27.9892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.850099,0.118891,-0.000174879,1.25168e-07,-3.06516e-11,-83456.8,32.6998], Tmin=(100,'K'), Tmax=(615.805,'K')), NASAPolynomial(coeffs=[14.3377,0.0378102,-2.01821e-05,4.03251e-09,-2.86149e-13,-85660.5,-35.9146], Tmin=(615.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(=C(F)F)C(F)(CF)OF(11398)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 C u1 p0 c0 {9,S} {14,S} {15,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-754.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42007,'amu*angstrom^2'), symmetry=1, barrier=(32.6501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4251,'amu*angstrom^2'), symmetry=1, barrier=(32.7658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42601,'amu*angstrom^2'), symmetry=1, barrier=(32.7869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42208,'amu*angstrom^2'), symmetry=1, barrier=(32.6964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729002,0.111124,-0.000146085,9.86845e-08,-2.66758e-11,-90522,31.5405], Tmin=(100,'K'), Tmax=(900.274,'K')), NASAPolynomial(coeffs=[16.63,0.033999,-1.75887e-05,3.53453e-09,-2.54325e-13,-93647.7,-50.3838], Tmin=(900.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-754.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Allyl_P)"""),
)

species(
    label = 'F[C]=C(CF)C(F)(CF)OF(11399)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-629.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,167,640,1190,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17894,'amu*angstrom^2'), symmetry=1, barrier=(27.1062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17717,'amu*angstrom^2'), symmetry=1, barrier=(27.0654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17324,'amu*angstrom^2'), symmetry=1, barrier=(26.9752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18103,'amu*angstrom^2'), symmetry=1, barrier=(27.1542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.832395,0.118614,-0.000173912,1.24299e-07,-3.03225e-11,-75552.5,32.7941], Tmin=(100,'K'), Tmax=(613.044,'K')), NASAPolynomial(coeffs=[14.1303,0.0384738,-2.06158e-05,4.12786e-09,-2.93355e-13,-77715.6,-34.751], Tmin=(613.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-629.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = '[O]C(F)(CF)C(CF)=C(F)F(11400)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-985.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2948,'amu*angstrom^2'), symmetry=1, barrier=(29.7701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29462,'amu*angstrom^2'), symmetry=1, barrier=(29.7658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29504,'amu*angstrom^2'), symmetry=1, barrier=(29.7755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.545785,0.108925,-0.000151838,1.13582e-07,-3.44878e-11,-118331,31.8871], Tmin=(100,'K'), Tmax=(800.416,'K')), NASAPolynomial(coeffs=[13.4911,0.0387792,-2.03881e-05,4.10196e-09,-2.94295e-13,-120578,-32.7085], Tmin=(800.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-985.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]F(128)',
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
    label = 'FC[C](F)C(CF)=C(F)F(11401)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-860.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,551,1088,1226,1380,1420,1481,3057,3119,350,440,435,1725,346,659,817,1284,182,240,577,636,1210,1413,247.238,247.254,247.256],'cm^-1')),
        HinderedRotor(inertia=(0.520088,'amu*angstrom^2'), symmetry=1, barrier=(22.5677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229299,'amu*angstrom^2'), symmetry=1, barrier=(9.94535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.986137,'amu*angstrom^2'), symmetry=1, barrier=(42.7666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24561,0.0901763,-0.000114559,7.98813e-08,-2.29879e-11,-103335,29.0172], Tmin=(100,'K'), Tmax=(837.077,'K')), NASAPolynomial(coeffs=[11.3969,0.0368898,-1.90727e-05,3.83406e-09,-2.75825e-13,-105202,-22.7981], Tmin=(837.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-860.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCdCsF1s)"""),
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
    label = 'FCC([C](F)F)=C(F)OF(10932)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-666.902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,326,540,652,719,1357,161,297,490,584,780,1358,180,2027.51],'cm^-1')),
        HinderedRotor(inertia=(0.337459,'amu*angstrom^2'), symmetry=1, barrier=(7.75884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75226,'amu*angstrom^2'), symmetry=1, barrier=(40.288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0138073,'amu*angstrom^2'), symmetry=1, barrier=(40.3529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125828,0.0961427,-0.000165027,1.53249e-07,-5.58818e-11,-80080.9,31.984], Tmin=(100,'K'), Tmax=(799.506,'K')), NASAPolynomial(coeffs=[8.20427,0.0374193,-2.05073e-05,4.10219e-09,-2.89554e-13,-80787.5,-1.52319], Tmin=(799.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-666.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFO) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = 'FC[C](F)OF(3020)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-334.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.263408,'amu*angstrom^2'), symmetry=1, barrier=(6.05627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4082,'amu*angstrom^2'), symmetry=1, barrier=(32.3773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91888,0.0521859,-8.94187e-05,8.43998e-08,-3.10903e-11,-40169.1,18.2725], Tmin=(100,'K'), Tmax=(810.437,'K')), NASAPolynomial(coeffs=[5.57801,0.0221788,-1.17676e-05,2.33426e-09,-1.64035e-13,-40369.8,3.80899], Tmin=(810.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s)"""),
)

species(
    label = 'FC[C]=C(F)F(9672)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-311.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,562,600,623,1070,1265,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.13897,'amu*angstrom^2'), symmetry=1, barrier=(3.19519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0431,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90019,0.026583,-1.69828e-05,4.78788e-09,-4.4109e-13,-37467.3,12.0998], Tmin=(10,'K'), Tmax=(1562.95,'K')), NASAPolynomial(coeffs=[13.6343,0.00613775,-1.648e-06,7.55221e-11,2.01777e-14,-41055.7,-40.9543], Tmin=(1562.95,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-311.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""FC[C]DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[CH]C(F)(OF)C(CF)=C(F)F(11402)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {10,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {3,S} {8,S} {15,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-885.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,334,575,1197,1424,3202,182,240,577,636,1210,1413,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13921,'amu*angstrom^2'), symmetry=1, barrier=(26.1927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15812,'amu*angstrom^2'), symmetry=1, barrier=(26.6275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1514,'amu*angstrom^2'), symmetry=1, barrier=(26.473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04718,'amu*angstrom^2'), symmetry=1, barrier=(24.0766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (193.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65788,0.136791,-0.000232254,2.03386e-07,-6.99032e-11,-106341,36.894], Tmin=(100,'K'), Tmax=(797.755,'K')), NASAPolynomial(coeffs=[14.8516,0.0393697,-2.15443e-05,4.29375e-09,-3.02048e-13,-108509,-36.1041], Tmin=(797.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-885.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCsF1sH)"""),
)

species(
    label = 'FCC(F)([C]=C(F)F)OF(11403)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-628.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33536,'amu*angstrom^2'), symmetry=1, barrier=(30.7025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33089,'amu*angstrom^2'), symmetry=1, barrier=(30.5998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33549,'amu*angstrom^2'), symmetry=1, barrier=(30.7056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0337222,0.099129,-0.000149011,1.07369e-07,-2.62303e-11,-75469,28.504], Tmin=(100,'K'), Tmax=(612.155,'K')), NASAPolynomial(coeffs=[12.8924,0.0298182,-1.63024e-05,3.27734e-09,-2.33015e-13,-77335.4,-29.832], Tmin=(612.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-628.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)(CF)OF(11404)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {11,S} {12,D}
11 C u1 p0 c0 {3,S} {10,S} {15,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-935.238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40807,'amu*angstrom^2'), symmetry=1, barrier=(32.3742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40739,'amu*angstrom^2'), symmetry=1, barrier=(32.3586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40797,'amu*angstrom^2'), symmetry=1, barrier=(32.3719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40717,'amu*angstrom^2'), symmetry=1, barrier=(32.3536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (193.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07007,0.121083,-0.000174574,1.30791e-07,-3.92799e-11,-112309,35.2852], Tmin=(100,'K'), Tmax=(812.583,'K')), NASAPolynomial(coeffs=[15.9272,0.037412,-2.01189e-05,4.07048e-09,-2.92568e-13,-115071,-43.1891], Tmin=(812.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-935.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C]C(F)(CF)C(F)(CF)OF(11405)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
11 C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
12 C u0 p1 c0 {5,S} {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-843.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,1380,1390,370,380,2900,435,442,614,1089,1143,1126,1238,1254,1408,1339,1465,1463,1525,3023,3127,3088,3132,617,898,1187,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (194.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47534,0.13043,-0.000191095,1.44633e-07,-4.37016e-11,-101238,36.3733], Tmin=(100,'K'), Tmax=(809.174,'K')), NASAPolynomial(coeffs=[17.1305,0.0384524,-2.05886e-05,4.1517e-09,-2.97577e-13,-104249,-49.4491], Tmin=(809.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-843.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCCF) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C(CF)C(CF)=C(F)F(11406)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u0 p0 c0 {5,D} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-902.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,244,1102,1254,1379,1468,1475,3037,3083,350,440,435,1725,375,552.5,462.5,1710,182,240,577,636,1210,1413,189.276,216.987,1761.99],'cm^-1')),
        HinderedRotor(inertia=(1.01574,'amu*angstrom^2'), symmetry=1, barrier=(29.5778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00588763,'amu*angstrom^2'), symmetry=1, barrier=(0.153116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197979,'amu*angstrom^2'), symmetry=1, barrier=(4.55194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273644,0.0915881,-0.00014284,1.31715e-07,-4.93654e-11,-108364,30.5083], Tmin=(100,'K'), Tmax=(764.49,'K')), NASAPolynomial(coeffs=[6.78153,0.0441426,-2.3466e-05,4.69754e-09,-3.33884e-13,-108967,3.41923], Tmin=(764.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-902.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCFF)"""),
)

species(
    label = 'C=C(OF)C(CF)=C(F)F(11407)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u0 p0 c0 {5,S} {7,S} {10,D}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 C u0 p0 c0 {8,D} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-568.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,182,240,577,636,1210,1413,2950,3100,1380,975,1025,1650,374.617,374.618,374.618],'cm^-1')),
        HinderedRotor(inertia=(0.0995561,'amu*angstrom^2'), symmetry=1, barrier=(9.91452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995562,'amu*angstrom^2'), symmetry=1, barrier=(9.91452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279246,'amu*angstrom^2'), symmetry=1, barrier=(27.8093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232945,0.0997482,-0.000140627,1.04867e-07,-3.1314e-11,-68177.7,29.0242], Tmin=(100,'K'), Tmax=(818.388,'K')), NASAPolynomial(coeffs=[13.7174,0.0315645,-1.56566e-05,3.06609e-09,-2.16389e-13,-70461.1,-35.4825], Tmin=(818.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-568.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFF) + group(Cds-CdsHH)"""),
)

species(
    label = 'FC=C(OF)C(CF)=C(F)F(11408)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {6,S} {8,S} {11,D}
10 C u0 p0 c0 {2,S} {3,S} {8,D}
11 C u0 p0 c0 {4,S} {9,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-733.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,325,375,415,465,420,450,1700,1750,182,240,577,636,1210,1413,194,682,905,1196,1383,3221,265.086,265.739,1524.22],'cm^-1')),
        HinderedRotor(inertia=(0.158799,'amu*angstrom^2'), symmetry=1, barrier=(7.9214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158412,'amu*angstrom^2'), symmetry=1, barrier=(7.91214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847092,'amu*angstrom^2'), symmetry=1, barrier=(42.3753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714835,0.114461,-0.000190298,1.67813e-07,-5.82127e-11,-88009.3,32.1976], Tmin=(100,'K'), Tmax=(807.456,'K')), NASAPolynomial(coeffs=[11.6772,0.0379181,-1.99536e-05,3.92597e-09,-2.74478e-13,-89516.4,-21.8776], Tmin=(807.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'C=C=C(F)F(10940)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {5,D} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-181.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92399,0.00482044,8.42754e-05,-1.86975e-07,1.24244e-10,-21859.1,8.85451], Tmin=(10,'K'), Tmax=(505.441,'K')), NASAPolynomial(coeffs=[3.60214,0.0228479,-1.51662e-05,4.78207e-09,-5.74586e-13,-22024.3,8.23155], Tmin=(505.441,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-181.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC[C]C(F)(CF)OF(11409)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
9  C u0 p1 c0 {6,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-347.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,249,734,1109,1255,1358,2983,3011,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0593179,0.0968454,-0.000135956,9.61699e-08,-2.36351e-11,-41615.9,36.9756], Tmin=(100,'K'), Tmax=(608.328,'K')), NASAPolynomial(coeffs=[11.2365,0.0364661,-1.94141e-05,3.89246e-09,-2.77574e-13,-43218.4,-13.3871], Tmin=(608.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-551.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-67.4919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-90.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-192.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-228.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-43.9425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-202.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-182.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-2.68427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-314.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-274.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-24.8149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-318.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-398.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-447.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-345.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-354.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-306.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-382.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-143.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-159.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-217.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-93.3188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-448.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-294.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-246.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-182.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-210.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-207.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-254.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-261.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-375.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-105.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-320.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-299.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-87.5081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    products = ['HF(38)', 'O=C(F)CF(879)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(65.8388,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]CF(126)', 'FCC(OF)=C(F)F(5237)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(173.439,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FCOF(367)', 'F[C]C(CF)=C(F)F(10927)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.14358e+07,'m^3/(mol*s)'), n=-0.641018, Ea=(193.845,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17570294811609, var=29.477358162589667, Tref=1000.0, N=2, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2(8)', 'F[C]C(F)(OF)C(CF)=C(F)F(11384)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHF(40)', 'FCC(=C(F)F)C(F)OF(11385)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.79957e-05,'m^3/(mol*s)'), n=3.10993, Ea=(22.0603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(S)(25)', 'FCC(=C(F)F)C(F)(F)OF(10923)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(160.638,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H2(8)', 'F[C]C(=C(F)F)C(F)(CF)OF(11386)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(85.3813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CHF(40)', 'FCC(F)(C=C(F)F)OF(4577)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.44768e+28,'m^3/(mol*s)'), n=-6.4458, Ea=(81.9127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.6591608693425184, var=5.4995193120720405, Tref=1000.0, N=19, data_mean=0.0, correlation='CH',), comment="""Estimated from node CH"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(S)(25)', 'FCC(F)(OF)C(F)=C(F)F(11387)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(140.159,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FCC(=C(F)F)C(F)(F)COF(11388)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OF(482)', 'FC=C(F)C(CF)=C(F)F(9874)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FOF(488)', 'C=C(F)C(CF)=C(F)F(11389)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2500,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    products = ['FCC(OF)=C(CF)C(F)(F)F(11390)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(298.965,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    products = ['FCC(F)=C(CF)C(F)(F)OF(11391)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(219.006,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FCC(=C(F)OF)C(F)(F)CF(11392)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(152.495,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    products = ['C=C(C(F)(F)F)C(F)(CF)OF(11393)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(272.428,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH]C([C](F)F)C(F)(CF)OF(7432)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH]C(F)(OF)C(CF)[C](F)F(11394)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(F)(OF)[C](CF)C(F)F(11395)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F(37)', 'FCC(OF)=C(CF)[C](F)F(11396)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F(37)', '[CH2]C(F)(OF)C(CF)=C(F)F(11397)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F(37)', '[CH2]C(=C(F)F)C(F)(CF)OF(11398)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F(37)', 'F[C]=C(CF)C(F)(CF)OF(11399)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F(37)', '[O]C(F)(CF)C(CF)=C(F)F(11400)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]F(128)', 'FC[C](F)C(CF)=C(F)F(11401)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2F(46)', 'FCC([C](F)F)=C(F)OF(10932)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -8.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['FC[C](F)OF(3020)', 'FC[C]=C(F)F(9672)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -26.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(5)', 'F[CH]C(F)(OF)C(CF)=C(F)F(11402)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2F(46)', 'FCC(F)([C]=C(F)F)OF(11403)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -19.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(5)', 'F[CH]C(=C(F)F)C(F)(CF)OF(11404)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(5.77755,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]C(F)(CF)C(F)(CF)OF(11405)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(118.897,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F2(78)', 'O=C(CF)C(CF)=C(F)F(11406)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(72.4901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F2(78)', 'C=C(OF)C(CF)=C(F)F(11407)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(8.22857,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction34',
    reactants = ['HF(38)', 'FC=C(OF)C(CF)=C(F)F(11408)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(230.74,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction35',
    reactants = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    products = ['F2(78)', 'O=C(F)CF(879)', 'C=C=C(F)F(10940)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(318.445,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CF2(43)', 'FC[C]C(F)(CF)OF(11409)'],
    products = ['FCC(=C(F)F)C(F)(CF)OF(11293)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.66533e+06,'cm^3/(mol*s)'), n=1.53, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for CF2
Exact match found for rate rule [CF2]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #2945',
    isomers = [
        'FCC(=C(F)F)C(F)(CF)OF(11293)',
    ],
    reactants = [
        ('HF(38)', 'O=C(F)CF(879)', 'FC=C=C(F)F(5101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2945',
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

