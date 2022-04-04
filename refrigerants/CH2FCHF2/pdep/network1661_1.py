species(
    label = 'O=C(F)OC(O)(F)CF(5377)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {7,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1238.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,482,664,788,1296,1923,306.078,306.089,306.103,306.123],'cm^-1')),
        HinderedRotor(inertia=(0.277992,'amu*angstrom^2'), symmetry=1, barrier=(18.4856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46166,'amu*angstrom^2'), symmetry=1, barrier=(30.6905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157277,'amu*angstrom^2'), symmetry=1, barrier=(10.4597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46158,'amu*angstrom^2'), symmetry=1, barrier=(30.691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660057,0.0774759,-9.15758e-05,5.48374e-08,-1.31836e-11,-148846,26.4392], Tmin=(100,'K'), Tmax=(1004.87,'K')), NASAPolynomial(coeffs=[13.9095,0.0247351,-1.2848e-05,2.60664e-09,-1.89192e-13,-151509,-37.546], Tmin=(1004.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1238.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COFOO)"""),
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
    label = 'CO(13)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-118.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2193.04],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(F)(CF)OF(5192)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-731.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19839,'amu*angstrom^2'), symmetry=1, barrier=(27.5534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19769,'amu*angstrom^2'), symmetry=1, barrier=(27.5373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19543,'amu*angstrom^2'), symmetry=1, barrier=(27.4852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74088,0.0168844,0.000203555,-5.93451e-07,4.85287e-10,-88000.4,11.3494], Tmin=(10,'K'), Tmax=(439.969,'K')), NASAPolynomial(coeffs=[7.09612,0.0354355,-2.69383e-05,9.22664e-09,-1.16633e-12,-88770.5,-7.47863], Tmin=(439.969,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-731.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)(CF)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)OC(O)(F)[C]F(6006)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {6,D}
9  C u0 p1 c0 {3,S} {7,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-886.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,482,664,788,1296,1923,617,898,1187,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.990305,'amu*angstrom^2'), symmetry=1, barrier=(22.7691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70374,'amu*angstrom^2'), symmetry=1, barrier=(39.1724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983441,'amu*angstrom^2'), symmetry=1, barrier=(22.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99064,'amu*angstrom^2'), symmetry=1, barrier=(22.7768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684465,0.0788566,-0.000115577,8.5261e-08,-2.49418e-11,-106535,27.2934], Tmin=(100,'K'), Tmax=(836.683,'K')), NASAPolynomial(coeffs=[12.8243,0.0208191,-1.15292e-05,2.35647e-09,-1.7028e-13,-108566,-29.1099], Tmin=(836.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-886.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + group(COFOO) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C(F)OC(O)F(5365)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-1020.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,482,664,788,1296,1923,222.971,223.169,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.265891,'amu*angstrom^2'), symmetry=1, barrier=(9.39111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817569,'amu*angstrom^2'), symmetry=1, barrier=(28.8891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608996,'amu*angstrom^2'), symmetry=1, barrier=(21.5108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3648.6,'J/mol'), sigma=(5.56031,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.90 K, Pc=48.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71583,0.0500006,-5.30385e-05,2.7775e-08,-5.77561e-12,-122651,21.0141], Tmin=(100,'K'), Tmax=(1161.05,'K')), NASAPolynomial(coeffs=[12.0233,0.0144896,-7.16047e-06,1.43206e-09,-1.03373e-13,-125045,-30.2527], Tmin=(1161.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1020.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO)"""),
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
    label = 'O=C(F)OC(O)(F)F(5371)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {7,S} {8,S}
5 O u0 p2 c0 {7,S} {9,S}
6 O u0 p2 c0 {8,D}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8 C u0 p0 c0 {3,S} {4,S} {6,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-1252.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,482,664,788,1296,1923,180,180,1673.6],'cm^-1')),
        HinderedRotor(inertia=(0.585359,'amu*angstrom^2'), symmetry=1, barrier=(13.4586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42777,'amu*angstrom^2'), symmetry=1, barrier=(32.8271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04463,'amu*angstrom^2'), symmetry=1, barrier=(24.0181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3648.54,'J/mol'), sigma=(5.76422,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.89 K, Pc=43.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49949,0.0562362,-6.36982e-05,3.52927e-08,-7.75565e-12,-150562,22.4514], Tmin=(100,'K'), Tmax=(1101,'K')), NASAPolynomial(coeffs=[12.7157,0.0154869,-8.1812e-06,1.6764e-09,-1.22484e-13,-153032,-32.7393], Tmin=(1101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1252.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFFOO) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OCC(O)(F)F(6007)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1250.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782498,0.0756341,-9.1622e-05,5.75721e-08,-1.4635e-11,-150284,26.8607], Tmin=(100,'K'), Tmax=(949.759,'K')), NASAPolynomial(coeffs=[12.5321,0.026151,-1.34729e-05,2.71806e-09,-1.96455e-13,-152516,-29.2187], Tmin=(949.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1250.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFFO) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OC(F)(F)CO(6008)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1250.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782506,0.0756341,-9.16217e-05,5.75718e-08,-1.46349e-11,-150284,26.8607], Tmin=(100,'K'), Tmax=(949.775,'K')), NASAPolynomial(coeffs=[12.5321,0.026151,-1.34729e-05,2.71806e-09,-1.96454e-13,-152516,-29.2188], Tmin=(949.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1250.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFFO) + group(Cs-CsOsHH) + group(COFOO)"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787932,2.90876e-06,-1.47517e-09,2.12837e-13,-30281.6,-0.311363], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81029e-07,1.07243e-10,-5.79388e-15,-29958.6,5.91041], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(F)OC(F)=CF(4076)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {2,S} {6,D} {9,S}
8 C u0 p0 c0 {3,S} {4,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-868.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,194,682,905,1196,1383,3221,482,664,788,1296,1923,180,180,266.032],'cm^-1')),
        HinderedRotor(inertia=(0.100201,'amu*angstrom^2'), symmetry=1, barrier=(2.30381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0242,'amu*angstrom^2'), symmetry=1, barrier=(23.5484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29258,0.0654495,-9.66411e-05,7.51369e-08,-2.3604e-11,-104410,21.7264], Tmin=(100,'K'), Tmax=(775.321,'K')), NASAPolynomial(coeffs=[9.75599,0.0217847,-1.21622e-05,2.49577e-09,-1.80659e-13,-105722,-16.9508], Tmin=(775.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-868.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COFOO)"""),
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
    label = 'C=C(F)OC(=O)F(6009)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {5,D} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-702.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,2950,3100,1380,975,1025,1650,482,664,788,1296,1923,180,180,776.271],'cm^-1')),
        HinderedRotor(inertia=(0.509238,'amu*angstrom^2'), symmetry=1, barrier=(11.7084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03417,'amu*angstrom^2'), symmetry=1, barrier=(23.7776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88769,0.0504497,-5.66512e-05,3.31093e-08,-7.95626e-12,-84359.2,18.8189], Tmin=(100,'K'), Tmax=(991.707,'K')), NASAPolynomial(coeffs=[9.43022,0.0200277,-1.06372e-05,2.17713e-09,-1.5866e-13,-85855.3,-17.5067], Tmin=(991.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-702.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(CdCFO) + group(Cds-CdsHH) + group(COFOO)"""),
)

species(
    label = '[O]C(F)OC([O])(F)CF(6010)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-757.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,451,553,637,1069,1180,1265,1301,3056,242.898,245.23,247.578,249.335],'cm^-1')),
        HinderedRotor(inertia=(0.341552,'amu*angstrom^2'), symmetry=1, barrier=(14.4197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166562,'amu*angstrom^2'), symmetry=1, barrier=(7.14781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966625,'amu*angstrom^2'), symmetry=1, barrier=(44.6313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579368,0.0816392,-0.000110953,8.07876e-08,-2.38723e-11,-91017.5,28.4557], Tmin=(100,'K'), Tmax=(821.804,'K')), NASAPolynomial(coeffs=[11.3717,0.0291126,-1.50843e-05,3.02128e-09,-2.1644e-13,-92791.4,-21.4935], Tmin=(821.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-757.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sCs)) + radical(O2sj(Cs-F1sO2sH))"""),
)

species(
    label = '[O]C(F)OC(O)(F)[CH]F(6011)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {12,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
9  C u1 p0 c0 {3,S} {7,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-815.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,451,553,637,1069,1180,1265,1301,3056,334,575,1197,1424,3202,257.021,257.056,257.062,257.063],'cm^-1')),
        HinderedRotor(inertia=(0.13762,'amu*angstrom^2'), symmetry=1, barrier=(6.45325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253569,'amu*angstrom^2'), symmetry=1, barrier=(11.8905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894172,'amu*angstrom^2'), symmetry=1, barrier=(41.9462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137584,'amu*angstrom^2'), symmetry=1, barrier=(6.45292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117223,0.0943734,-0.00014189,1.02663e-07,-2.58744e-11,-97939,29.874], Tmin=(100,'K'), Tmax=(634.034,'K')), NASAPolynomial(coeffs=[13.1426,0.0261961,-1.37101e-05,2.70478e-09,-1.90113e-13,-99872,-29.2495], Tmin=(634.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-815.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = '[O]C(F)(CF)O[C](O)F(6012)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-808.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,482,586,761,1411,286.197,286.199,286.204,286.207],'cm^-1')),
        HinderedRotor(inertia=(0.202295,'amu*angstrom^2'), symmetry=1, barrier=(11.7579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.487124,'amu*angstrom^2'), symmetry=1, barrier=(28.318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202286,'amu*angstrom^2'), symmetry=1, barrier=(11.7578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662183,'amu*angstrom^2'), symmetry=1, barrier=(38.489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.382892,0.0852134,-0.000114154,7.84246e-08,-2.15309e-11,-97107.4,29.4969], Tmin=(100,'K'), Tmax=(887.336,'K')), NASAPolynomial(coeffs=[13.6128,0.0255767,-1.33446e-05,2.68803e-09,-1.93533e-13,-99455.3,-32.749], Tmin=(887.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-808.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = 'O[C](F)OC(O)(F)[CH]F(6013)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {9,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u1 p0 c0 {2,S} {7,S} {10,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-866.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,364,474,579,619,655,1345,334,575,1197,1424,3202,482,586,761,1411,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.280571,0.100843,-0.000158438,1.24361e-07,-3.82714e-11,-104020,31.6087], Tmin=(100,'K'), Tmax=(800.405,'K')), NASAPolynomial(coeffs=[15.3531,0.0227077,-1.19958e-05,2.37715e-09,-1.67649e-13,-106523,-40.333], Tmin=(800.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-866.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHOO) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
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
    label = 'O=C(F)O[C](O)CF(6014)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {4,S} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-814.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,360,370,350,482,664,788,1296,1923,316.138,318.002,319.587,999.353],'cm^-1')),
        HinderedRotor(inertia=(0.00163697,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115736,'amu*angstrom^2'), symmetry=1, barrier=(8.29769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346244,'amu*angstrom^2'), symmetry=1, barrier=(24.4639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344809,'amu*angstrom^2'), symmetry=1, barrier=(24.4617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30689,0.0630789,-7.08256e-05,4.10012e-08,-9.65799e-12,-97812.7,26.6787], Tmin=(100,'K'), Tmax=(1017.04,'K')), NASAPolynomial(coeffs=[11.4509,0.023182,-1.19817e-05,2.42846e-09,-1.76163e-13,-99876,-22.4314], Tmin=(1017.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(COFOO) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C(O)(F)OC(=O)F(6015)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {6,S} {11,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-849.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,3000,3100,440,815,1455,1000,482,664,788,1296,1923,253.076,253.1,253.112],'cm^-1')),
        HinderedRotor(inertia=(0.184194,'amu*angstrom^2'), symmetry=1, barrier=(8.37347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696132,'amu*angstrom^2'), symmetry=1, barrier=(31.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00263162,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696066,'amu*angstrom^2'), symmetry=1, barrier=(31.6402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.840005,0.0747886,-9.97039e-05,6.85926e-08,-1.89313e-11,-102100,24.9359], Tmin=(100,'K'), Tmax=(881.294,'K')), NASAPolynomial(coeffs=[12.1514,0.0234483,-1.23201e-05,2.48981e-09,-1.79593e-13,-104093,-28.2059], Tmin=(881.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-849.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsHHH) + group(COFOO) + radical(Csj(Cs-F1sO2sO2s)(H)(H))"""),
)

species(
    label = 'O=[C]OC(O)(F)CF(6016)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {6,S} {11,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {3,S} {5,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-819.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,1855,455,950,228.066,228.091,228.124],'cm^-1')),
        HinderedRotor(inertia=(1.30706,'amu*angstrom^2'), symmetry=1, barrier=(48.239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.546248,'amu*angstrom^2'), symmetry=1, barrier=(20.1732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30704,'amu*angstrom^2'), symmetry=1, barrier=(48.2391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601909,'amu*angstrom^2'), symmetry=1, barrier=(22.2218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.942648,0.0720422,-9.08026e-05,5.96772e-08,-1.58402e-11,-98456.5,23.367], Tmin=(100,'K'), Tmax=(912.496,'K')), NASAPolynomial(coeffs=[11.7978,0.0244569,-1.25782e-05,2.52556e-09,-1.81804e-13,-100438,-28.0084], Tmin=(912.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-819.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(=O)F(341)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)CF(2940)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {8,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-471.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180],'cm^-1')),
        HinderedRotor(inertia=(0.625922,'amu*angstrom^2'), symmetry=1, barrier=(14.3912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628429,'amu*angstrom^2'), symmetry=1, barrier=(14.4488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78656,0.0175844,0.000121424,-4.80809e-07,5.06984e-10,-56694.5,10.4186], Tmin=(10,'K'), Tmax=(351.977,'K')), NASAPolynomial(coeffs=[6.00134,0.0208372,-1.35648e-05,4.29025e-09,-5.21054e-13,-57026.5,-0.454778], Tmin=(351.977,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-471.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""O[C](F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CFO(51)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)(F)CF(3973)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6 C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-628.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,180],'cm^-1')),
        HinderedRotor(inertia=(1.00515,'amu*angstrom^2'), symmetry=1, barrier=(23.1104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00489,'amu*angstrom^2'), symmetry=1, barrier=(23.1043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78076,0.0150815,0.000141095,-3.93332e-07,3.04195e-10,-75641.9,11.4432], Tmin=(10,'K'), Tmax=(466.103,'K')), NASAPolynomial(coeffs=[6.42151,0.0289801,-2.12929e-05,7.21919e-09,-9.08945e-13,-76285.2,-3.54139], Tmin=(466.103,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-628.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(O)(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(F)O[C](F)CF(3266)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-835.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,482,664,788,1296,1923,329.436,329.437,329.437,329.437],'cm^-1')),
        HinderedRotor(inertia=(0.505038,'amu*angstrom^2'), symmetry=1, barrier=(38.895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113551,'amu*angstrom^2'), symmetry=1, barrier=(8.74503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377194,'amu*angstrom^2'), symmetry=1, barrier=(29.0486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4682,0.0608494,-7.31499e-05,4.67426e-08,-1.22986e-11,-100405,23.601], Tmin=(100,'K'), Tmax=(910.997,'K')), NASAPolynomial(coeffs=[9.8081,0.0242308,-1.28558e-05,2.61957e-09,-1.90196e-13,-101925,-15.8568], Tmin=(910.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-835.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COFOO) + radical(CsCsF1sO2s)"""),
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
    label = '[O]C(F)(CF)OC(=O)F(6017)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-979.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,482,664,788,1296,1923,323.497,323.497,323.497,323.497],'cm^-1')),
        HinderedRotor(inertia=(0.147601,'amu*angstrom^2'), symmetry=1, barrier=(10.9612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458382,'amu*angstrom^2'), symmetry=1, barrier=(34.0405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458382,'amu*angstrom^2'), symmetry=1, barrier=(34.0405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01492,0.0713049,-8.65506e-05,5.48583e-08,-1.41915e-11,-117749,25.9882], Tmin=(100,'K'), Tmax=(929.44,'K')), NASAPolynomial(coeffs=[11.4493,0.0263989,-1.40777e-05,2.87498e-09,-2.09028e-13,-119688,-23.588], Tmin=(929.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-979.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COFOO) + radical(O2sj(Cs-F1sO2sCs))"""),
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
    label = 'O=C(F)O[C](O)F(5802)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-814.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,586,761,1411,482,664,788,1296,1923,189.176,1054.98,1905.38],'cm^-1')),
        HinderedRotor(inertia=(0.583046,'amu*angstrom^2'), symmetry=1, barrier=(14.7893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08391,'amu*angstrom^2'), symmetry=1, barrier=(27.4997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58374,'amu*angstrom^2'), symmetry=1, barrier=(40.1069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00526,0.0434888,-4.46044e-05,2.23764e-08,-4.46436e-12,-97863,22.1567], Tmin=(100,'K'), Tmax=(1205,'K')), NASAPolynomial(coeffs=[11.2161,0.012913,-6.54267e-06,1.31844e-09,-9.54185e-14,-100083,-23.9976], Tmin=(1205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'O=C(F)OC(O)(F)[CH]F(6018)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8  C u1 p0 c0 {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1037.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,482,664,788,1296,1923,320.218,320.218,320.218,320.218],'cm^-1')),
        HinderedRotor(inertia=(0.10086,'amu*angstrom^2'), symmetry=1, barrier=(7.33901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10086,'amu*angstrom^2'), symmetry=1, barrier=(7.33901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45515,'amu*angstrom^2'), symmetry=1, barrier=(33.1188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45515,'amu*angstrom^2'), symmetry=1, barrier=(33.1188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346992,0.0869297,-0.000130541,1.00029e-07,-3.04389e-11,-124661,28.1196], Tmin=(100,'K'), Tmax=(805.125,'K')), NASAPolynomial(coeffs=[13.0584,0.0237764,-1.28814e-05,2.602e-09,-1.86407e-13,-126708,-30.4502], Tmin=(805.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1037.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COFOO) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'O=C(O)F(343)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)=CF(878)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {2,S} {4,D} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-508.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.537415,'amu*angstrom^2'), symmetry=1, barrier=(12.3562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3292.1,'J/mol'), sigma=(5.16209,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.22 K, Pc=54.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86684,0.00827899,9.49292e-05,-2.30798e-07,1.56874e-10,-61152.3,9.15162], Tmin=(10,'K'), Tmax=(532.988,'K')), NASAPolynomial(coeffs=[6.38193,0.0179589,-1.26771e-05,4.31723e-09,-5.56992e-13,-61826,-5.20451], Tmin=(532.988,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-508.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)OC(=O)CF(6019)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-927.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,482,664,788,1296,1923,180,180,180,180,1600,1849.78,2665.76,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163803,'amu*angstrom^2'), symmetry=1, barrier=(3.76615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163803,'amu*angstrom^2'), symmetry=1, barrier=(3.76615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163803,'amu*angstrom^2'), symmetry=1, barrier=(3.76615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41372,0.061843,-8.65607e-05,6.51757e-08,-1.989e-11,-111407,26.2882], Tmin=(100,'K'), Tmax=(797.245,'K')), NASAPolynomial(coeffs=[9.38016,0.0218691,-1.13433e-05,2.2717e-09,-1.62548e-13,-112677,-10.3393], Tmin=(797.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(COFOO)"""),
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
    label = 'C=C(O)OC(=O)F(6020)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,D}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-673.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,482,664,788,1296,1923,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.811387,'amu*angstrom^2'), symmetry=1, barrier=(18.6554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810385,'amu*angstrom^2'), symmetry=1, barrier=(18.6323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80981,'amu*angstrom^2'), symmetry=1, barrier=(18.6191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928027,0.0657785,-8.1053e-05,4.78981e-08,-1.09205e-11,-80943.9,20.6301], Tmin=(100,'K'), Tmax=(1081.52,'K')), NASAPolynomial(coeffs=[15.7761,0.010863,-4.88863e-06,9.49113e-10,-6.79352e-14,-84155.6,-52.1667], Tmin=(1081.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-673.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OC(O)=CF(6021)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-838.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,194,682,905,1196,1383,3221,482,664,788,1296,1923,180,180,594.532],'cm^-1')),
        HinderedRotor(inertia=(0.586396,'amu*angstrom^2'), symmetry=1, barrier=(13.4824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753661,'amu*angstrom^2'), symmetry=1, barrier=(17.3282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753688,'amu*angstrom^2'), symmetry=1, barrier=(17.3288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.792115,0.0759801,-0.000112817,8.43377e-08,-2.49177e-11,-100790,22.59], Tmin=(100,'K'), Tmax=(830.001,'K')), NASAPolynomial(coeffs=[12.5458,0.0193363,-1.04493e-05,2.11487e-09,-1.51864e-13,-102741,-31.9246], Tmin=(830.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COFOO)"""),
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
    E0 = (-612.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-219.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-274.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-318.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-130.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-399.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-399.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-349.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-54.7239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-193.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-250.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-243.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-301.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-201.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-237.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-206.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-313.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-279.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-267.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-228.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-317.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-286.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-544.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-496.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-514.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-138.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-373.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)OC(O)(F)CF(5377)'],
    products = ['HF(38)', 'CO2(14)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.29418e+10,'s^-1'), n=0.745914, Ea=(86.0745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'OC(F)(CF)OF(5192)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(91.6069,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H2(8)', 'O=C(F)OC(O)(F)[C]F(6006)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(40)', 'O=C(F)OC(O)F(5365)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.79957e-05,'m^3/(mol*s)'), n=3.10993, Ea=(23.7398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(S)(25)', 'O=C(F)OC(O)(F)F(5371)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(163.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)OC(O)(F)CF(5377)'],
    products = ['O=C(F)OCC(O)(F)F(6007)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)OC(O)(F)CF(5377)'],
    products = ['O=C(F)OC(F)(F)CO(6008)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H2O(3)', 'O=C(F)OC(F)=CF(4076)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OF(482)', 'C=C(F)OC(=O)F(6009)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1250,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd/unsub_Cd/disub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)OC([O])(F)CF(6010)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)OC(O)(F)[CH]F(6011)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)(CF)O[C](O)F(6012)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O[C](F)OC(O)(F)[CH]F(6013)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C(F)O[C](O)CF(6014)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[CH2]C(O)(F)OC(=O)F(6015)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=[C]OC(O)(F)CF(6016)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(=O)F(341)', 'O[C](F)CF(2940)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(51)', '[O]C(O)(F)CF(3973)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(7)', 'O=C(F)O[C](F)CF(3266)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(F)(CF)OC(=O)F(6017)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2F(46)', 'O=C(F)O[C](O)F(5802)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'O=C(F)OC(O)(F)[CH]F(6018)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -1.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(F)OC(O)(F)CF(5377)'],
    products = ['O=C(O)F(343)', 'O=C(F)CF(879)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.82035e+06,'s^-1'), n=1.67734, Ea=(154.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(F)OC(O)(F)CF(5377)'],
    products = ['O=C(O)F(343)', 'OC(F)=CF(878)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.91252e+12,'s^-1'), n=0, Ea=(202.327,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_Ext-2R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_Ext-2R!H-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'O=C(F)OC(=O)CF(6019)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(154.107,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F2(78)', 'C=C(O)OC(=O)F(6020)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(4.65313,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', 'O=C(F)OC(O)=CF(6021)'],
    products = ['O=C(F)OC(O)(F)CF(5377)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(207.011,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #1661',
    isomers = [
        'O=C(F)OC(O)(F)CF(5377)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1661',
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

