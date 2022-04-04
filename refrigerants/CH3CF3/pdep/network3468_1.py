species(
    label = 'FOC(F)(F)C1C(F)C1(F)F(11638)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {8,S} {9,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-1033.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,180,180,180,775.902,775.914,775.923],'cm^-1')),
        HinderedRotor(inertia=(0.0123734,'amu*angstrom^2'), symmetry=1, barrier=(5.28645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111055,'amu*angstrom^2'), symmetry=1, barrier=(47.4418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11578,0.0886657,-0.000101708,5.85191e-08,-1.34627e-11,-124223,28.2933], Tmin=(100,'K'), Tmax=(1050.2,'K')), NASAPolynomial(coeffs=[16.2446,0.0272355,-1.39693e-05,2.82376e-09,-2.04715e-13,-127611,-50.3095], Tmin=(1050.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1033.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs)"""),
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
    label = 'CF2O(49)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CC1(F)F(974)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-329.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,323,467,575,827,1418,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.45,'J/mol'), sigma=(5.0155,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.04 K, Pc=52.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90533,0.00560985,9.08669e-05,-1.89954e-07,1.15706e-10,-39564.3,10.133], Tmin=(10,'K'), Tmax=(565.045,'K')), NASAPolynomial(coeffs=[4.14369,0.0253568,-1.84551e-05,6.16339e-09,-7.67829e-13,-39933.5,6.09131], Tmin=(565.045,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-329.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC=CC(F)(F)OF(354)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {9,S}
8  C u0 p0 c0 {3,S} {7,D} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-651.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,242.553,246.912],'cm^-1')),
        HinderedRotor(inertia=(0.236378,'amu*angstrom^2'), symmetry=1, barrier=(10.0214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763594,'amu*angstrom^2'), symmetry=1, barrier=(31.0347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73368,0.0181252,0.000172678,-4.71502e-07,3.57348e-10,-78385.6,13.6687], Tmin=(10,'K'), Tmax=(474.15,'K')), NASAPolynomial(coeffs=[6.82303,0.0367988,-2.79207e-05,9.53211e-09,-1.19983e-12,-79181.5,-4.23317], Tmin=(474.15,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-651.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""FCDCC(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)C=C(F)F(1384)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-858.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,267.751,268.886],'cm^-1')),
        HinderedRotor(inertia=(0.3004,'amu*angstrom^2'), symmetry=1, barrier=(15.2411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636275,'amu*angstrom^2'), symmetry=1, barrier=(32.4498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.6881,0.0211397,0.000187466,-5.21693e-07,3.96244e-10,-103252,14.5945], Tmin=(10,'K'), Tmax=(478.832,'K')), NASAPolynomial(coeffs=[7.96478,0.0379937,-3.00443e-05,1.0468e-08,-1.33235e-12,-104265,-9.18287], Tmin=(478.832,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-858.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""FOC(F)(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC1C(F)C1(F)F(11812)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-535.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,250,417,511,1155,1315,1456,3119,222,329,445,522,589,1214,1475,528.834,528.882,528.884,528.909,528.947],'cm^-1')),
        HinderedRotor(inertia=(0.000602744,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4883,0.0581569,-6.22119e-05,3.4446e-08,-7.7835e-12,-64317.8,19.9702], Tmin=(100,'K'), Tmax=(1056.31,'K')), NASAPolynomial(coeffs=[10.9778,0.0222228,-1.11846e-05,2.24153e-09,-1.61642e-13,-66322.6,-26.3309], Tmin=(1056.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-535.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs(O2))"""),
)

species(
    label = 'FOF(328)',
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
    label = 'F[C]C1C(F)C1(F)F-2(9908)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p1 c0 {4,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-388.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,250,417,511,1155,1315,1456,3119,222,329,445,522,589,1214,1475,617,898,1187,180,953.423,954.327,955.077],'cm^-1')),
        HinderedRotor(inertia=(0.109187,'amu*angstrom^2'), symmetry=1, barrier=(2.51042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72724,0.0519603,-4.92345e-05,2.3956e-08,-4.79624e-12,-46677.1,21.2052], Tmin=(100,'K'), Tmax=(1174.08,'K')), NASAPolynomial(coeffs=[10.5464,0.0219142,-1.08476e-05,2.15906e-09,-1.54962e-13,-48748,-22.7573], Tmin=(1174.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-388.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CJ2_singlet-FCs) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + ring(Cs(F)-Cs-Cs)"""),
)

species(
    label = 'OF(195)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2889.11,'J/mol'), sigma=(4.75593,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=451.27 K, Pc=60.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=C1C(F)C1(F)F(3033)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-758.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,284,328,853,1146,1135,1297,3239,182,240,577,636,1210,1413,180,652.645,1539.85,1540.11],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57343,0.0440084,-2.36286e-05,-7.86367e-09,8.39801e-12,-91251.3,13.6077], Tmin=(10,'K'), Tmax=(860.316,'K')), NASAPolynomial(coeffs=[10.0944,0.0269715,-1.70818e-05,5.00836e-09,-5.57198e-13,-92864.9,-19.7278], Tmin=(860.316,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-758.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC(F)DC1C(F)C1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C([C](F)F)C(F)(F)OF(3640)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u1 p0 c0 {3,S} {8,S} {13,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-846.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,190,488,555,1236,1407,262.908,262.908,262.908,2141.51],'cm^-1')),
        HinderedRotor(inertia=(0.196949,'amu*angstrom^2'), symmetry=1, barrier=(9.66026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196948,'amu*angstrom^2'), symmetry=1, barrier=(9.66024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.630458,'amu*angstrom^2'), symmetry=1, barrier=(30.9238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.630458,'amu*angstrom^2'), symmetry=1, barrier=(30.9238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.947429,0.119568,-0.000204968,1.79281e-07,-6.13582e-11,-101627,35.1782], Tmin=(100,'K'), Tmax=(800.848,'K')), NASAPolynomial(coeffs=[13.8847,0.0326209,-1.80179e-05,3.5954e-09,-2.52812e-13,-103590,-30.5084], Tmin=(800.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-846.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCFFO) + group(CsCsFHH) + group(CsCsFFH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]C(F)(F)OF(11813)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 C u1 p0 c0 {5,S} {8,S} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-844.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,334,575,1197,1424,3202,232.812,233.083,233.454,1576.1],'cm^-1')),
        HinderedRotor(inertia=(0.217403,'amu*angstrom^2'), symmetry=1, barrier=(8.39771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218218,'amu*angstrom^2'), symmetry=1, barrier=(8.40088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216855,'amu*angstrom^2'), symmetry=1, barrier=(8.4116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809457,'amu*angstrom^2'), symmetry=1, barrier=(31.3434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07597,0.122927,-0.000214265,1.8791e-07,-6.40419e-11,-101348,35.7651], Tmin=(100,'K'), Tmax=(813.339,'K')), NASAPolynomial(coeffs=[14.4012,0.0315131,-1.74637e-05,3.47471e-09,-2.43331e-13,-103360,-32.5953], Tmin=(813.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-844.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + group(CsCsFHH) + radical(CCJCO) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FOC(F)(F)[CH]C(F)[C](F)F(11814)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {13,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-831.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,190,488,555,1236,1407,208.913,210.92,211.197,4000],'cm^-1')),
        HinderedRotor(inertia=(0.439996,'amu*angstrom^2'), symmetry=1, barrier=(13.7054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224514,'amu*angstrom^2'), symmetry=1, barrier=(6.94519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.066,'amu*angstrom^2'), symmetry=1, barrier=(34.5034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08285,'amu*angstrom^2'), symmetry=1, barrier=(34.5009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.839455,0.115853,-0.000192842,1.64543e-07,-5.53913e-11,-99840.7,36.1163], Tmin=(100,'K'), Tmax=(776.481,'K')), NASAPolynomial(coeffs=[14.4582,0.03093,-1.69683e-05,3.39402e-09,-2.39669e-13,-102032,-32.6282], Tmin=(776.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-831.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(CCJCO) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = 'FOC(F)(F)C1[C](F)C1F(11815)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-611.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,212,367,445,1450,180,180,829.195,831.086,831.11,2412.27],'cm^-1')),
        HinderedRotor(inertia=(1.44519,'amu*angstrom^2'), symmetry=1, barrier=(33.2278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211323,'amu*angstrom^2'), symmetry=1, barrier=(33.2307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.535884,0.0795243,-9.03343e-05,5.16315e-08,-1.18439e-11,-73420.7,26.829], Tmin=(100,'K'), Tmax=(1050.66,'K')), NASAPolynomial(coeffs=[14.7326,0.0254758,-1.31711e-05,2.67005e-09,-1.93813e-13,-76403.9,-42.3637], Tmin=(1050.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFFO) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(F1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FOC(F)(F)C1[CH]C1(F)F(11816)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
10 C u1 p0 c0 {7,S} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-620.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,455.673,459.065,461.028,461.273,461.566,461.838,462.007,462.274],'cm^-1')),
        HinderedRotor(inertia=(0.000791674,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000790785,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364666,0.0822794,-9.50865e-05,5.4609e-08,-1.24716e-11,-74499.1,27.5552], Tmin=(100,'K'), Tmax=(1060.77,'K')), NASAPolynomial(coeffs=[15.9307,0.0235818,-1.20836e-05,2.44346e-09,-1.77263e-13,-77801.5,-48.4603], Tmin=(1060.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-620.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(H)_ring)"""),
)

species(
    label = 'FO[C](F)C1C(F)C1(F)F(11817)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u1 p0 c0 {4,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-604.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,250,417,511,1155,1315,1456,3119,222,329,445,522,589,1214,1475,395,473,707,1436,435.95,435.951,435.953,435.955,435.956,435.957],'cm^-1')),
        HinderedRotor(inertia=(0.000887014,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000886979,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475667,0.0835962,-0.000108352,7.3866e-08,-2.03775e-11,-72610.5,27.6047], Tmin=(100,'K'), Tmax=(878.66,'K')), NASAPolynomial(coeffs=[12.4924,0.0288911,-1.49617e-05,3.00726e-09,-2.16355e-13,-74722.2,-28.8144], Tmin=(878.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCFHO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)C1C(F)C1(F)F(11818)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-908.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,250,417,511,1155,1315,1456,3119,222,329,445,522,589,1214,1475,351,323,533,609,664,892,1120,1201,180,551.87,1141.08,1141.43,1141.55],'cm^-1')),
        HinderedRotor(inertia=(0.203453,'amu*angstrom^2'), symmetry=1, barrier=(4.67779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723673,0.0733764,-7.78844e-05,4.17783e-08,-9.00009e-12,-109154,26.5231], Tmin=(100,'K'), Tmax=(1115.7,'K')), NASAPolynomial(coeffs=[14.4036,0.0243303,-1.19435e-05,2.3758e-09,-1.70851e-13,-112207,-40.9722], Tmin=(1115.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-908.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCFFO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(O2sj(Cs-CsF1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]F(127)',
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
    collisionModel = TransportData(shapeIndex=1, epsilon=(2889.11,'J/mol'), sigma=(4.75593,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=451.27 K, Pc=60.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)C1C(F)C1(F)F(9346)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-758.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,190,488,555,1236,1407,630.499,630.536,630.538,630.573,630.729],'cm^-1')),
        HinderedRotor(inertia=(0.000424134,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.051,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2975.88,'J/mol'), sigma=(5.17347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.83 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56124,0.042768,-2.59143e-06,-4.35885e-08,2.70195e-11,-91177.6,14.2639], Tmin=(10,'K'), Tmax=(746.717,'K')), NASAPolynomial(coeffs=[8.52041,0.0344153,-2.23976e-05,6.75748e-09,-7.71863e-13,-92425.9,-11.6127], Tmin=(746.717,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-758.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), label="""F[C](F)C1C(F)C1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FO[C](F)F(236)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
"""),
    E0 = (-317.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.523179,'amu*angstrom^2'), symmetry=1, barrier=(12.0289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66912,0.0316874,-4.80467e-05,3.64259e-08,-1.08808e-11,-38142.7,14.3743], Tmin=(100,'K'), Tmax=(821.789,'K')), NASAPolynomial(coeffs=[7.60337,0.00767068,-4.20998e-06,8.64397e-10,-6.26025e-14,-38953.7,-8.46215], Tmin=(821.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsF1sF1sO2s)"""),
)

species(
    label = 'FC1[CH]C1(F)F(5918)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
6 C u1 p0 c0 {4,S} {5,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-285.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,2950,1000,779.861,1225.17],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0431,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88868,0.00674243,0.000107888,-2.31767e-07,1.45664e-10,-34351,11.6691], Tmin=(10,'K'), Tmax=(546.414,'K')), NASAPolynomial(coeffs=[4.14002,0.0292647,-2.08185e-05,6.86042e-09,-8.47418e-13,-34742.1,7.28048], Tmin=(546.414,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-285.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1[CH]C1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)[C]1C(F)C1(F)F(11819)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-881.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,253,525,597,667,842,1178,1324,210.762,323.019,324.359,1364.66,1366.6],'cm^-1')),
        HinderedRotor(inertia=(1.88255,'amu*angstrom^2'), symmetry=1, barrier=(43.2835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582061,'amu*angstrom^2'), symmetry=1, barrier=(43.2831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371306,0.0822511,-9.13644e-05,5.03064e-08,-1.10693e-11,-105881,28.2613], Tmin=(100,'K'), Tmax=(1095.12,'K')), NASAPolynomial(coeffs=[16.0123,0.0251207,-1.31115e-05,2.66864e-09,-1.94158e-13,-109307,-48.6188], Tmin=(1095.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-881.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCFFO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CCJ(C)CO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FOC(F)(F)C1[C](F)C1(F)F(11820)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-809.352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,212,367,445,1450,414.337,414.342,414.343,414.346,414.348,1343.86],'cm^-1')),
        HinderedRotor(inertia=(0.294828,'amu*angstrom^2'), symmetry=1, barrier=(35.9188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0067558,'amu*angstrom^2'), symmetry=1, barrier=(35.9188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0652457,0.091777,-0.000117655,7.54421e-08,-1.92425e-11,-97205.2,29.0883], Tmin=(100,'K'), Tmax=(953.806,'K')), NASAPolynomial(coeffs=[15.824,0.0256906,-1.37275e-05,2.80325e-09,-2.03832e-13,-100211,-46.1939], Tmin=(953.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-809.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'O=C(F)C1C(F)C1(F)F(11821)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-890.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,250,417,511,1155,1315,1456,3119,222,329,445,522,589,1214,1475,486,617,768,1157,1926,180,180,1137.77,1137.84,1137.84],'cm^-1')),
        HinderedRotor(inertia=(0.108315,'amu*angstrom^2'), symmetry=1, barrier=(2.49037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31089,0.0605806,-6.03111e-05,3.06249e-08,-6.31749e-12,-107053,23.5917], Tmin=(100,'K'), Tmax=(1153.1,'K')), NASAPolynomial(coeffs=[12.1978,0.0228149,-1.11838e-05,2.22191e-09,-1.59519e-13,-109563,-30.4822], Tmin=(1153.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(COCsFO) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FOC(F)(F)C1=C(F)C1F(11822)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-622.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,275,321,533,585,746,850,1103,323,467,575,827,1418,180,180,180,1153.68],'cm^-1')),
        HinderedRotor(inertia=(0.81569,'amu*angstrom^2'), symmetry=1, barrier=(18.7543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07612,'amu*angstrom^2'), symmetry=1, barrier=(47.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.389423,0.0868657,-0.000128866,1.00512e-07,-3.15771e-11,-74772.1,26.6329], Tmin=(100,'K'), Tmax=(776.471,'K')), NASAPolynomial(coeffs=[11.7808,0.0281772,-1.54795e-05,3.15085e-09,-2.26899e-13,-76540.9,-25.441], Tmin=(776.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cs-Cd(C)-Cd(F))"""),
)

species(
    label = 'FOC(F)(F)C1=CC1(F)F(11781)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-680.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,275,321,533,585,746,850,1103,2950,1000,301.482,301.518,301.518,301.52,301.526,301.542,301.557],'cm^-1')),
        HinderedRotor(inertia=(0.632661,'amu*angstrom^2'), symmetry=1, barrier=(40.8446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633594,'amu*angstrom^2'), symmetry=1, barrier=(40.845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499558,0.0810741,-0.000102267,6.43106e-08,-1.60615e-11,-81710.3,26.6008], Tmin=(100,'K'), Tmax=(974.189,'K')), NASAPolynomial(coeffs=[14.8073,0.0223257,-1.18076e-05,2.40517e-09,-1.74786e-13,-84497.9,-42.0511], Tmin=(974.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-680.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFF) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd-Cd(C)-Cs(F)(F))"""),
)

species(
    label = 'FOC(F)=C1C(F)C1(F)F(11823)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-599.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,284,328,853,1146,1135,1297,3239,274,345,380,539,705,1166,1213,326,540,652,719,1357,180,180,1429.65,1430.38,1431.04],'cm^-1')),
        HinderedRotor(inertia=(0.303443,'amu*angstrom^2'), symmetry=1, barrier=(6.97675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.730219,0.0806577,-0.000113526,7.94307e-08,-1.86758e-11,-72011.4,26.709], Tmin=(100,'K'), Tmax=(600.556,'K')), NASAPolynomial(coeffs=[10.0316,0.0301992,-1.62033e-05,3.26038e-09,-2.32925e-13,-73335.9,-15.1475], Tmin=(600.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + ring(Cs-Cs(F)(F)-Cd(Cd))"""),
)

species(
    label = 'FOC(F)(F)C1C=C1F(11824)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-413.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,223,363,546,575,694,1179,1410,323,467,575,827,1418,277.736,277.737,277.737,277.737,277.738,1911.85,1911.86],'cm^-1')),
        HinderedRotor(inertia=(0.00218539,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513747,'amu*angstrom^2'), symmetry=1, barrier=(28.1219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660814,0.079557,-0.000110547,8.00239e-08,-2.32745e-11,-49598.5,22.454], Tmin=(100,'K'), Tmax=(837.471,'K')), NASAPolynomial(coeffs=[12.0119,0.0253389,-1.3433e-05,2.71389e-09,-1.953e-13,-51499.7,-30.2949], Tmin=(837.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'FOC(F)(F)C1C(F)=C1F(11825)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-578.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,223,363,546,575,694,1179,1410,260,386,409,525,515,635,761,893,1354,1482,180,180,180,882.09],'cm^-1')),
        HinderedRotor(inertia=(0.244152,'amu*angstrom^2'), symmetry=1, barrier=(5.61354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242232,'amu*angstrom^2'), symmetry=1, barrier=(5.56939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.196529,0.0933351,-0.000142235,1.02424e-07,-2.47622e-11,-69422.2,24.2416], Tmin=(100,'K'), Tmax=(615.322,'K')), NASAPolynomial(coeffs=[12.8327,0.0260057,-1.42156e-05,2.84847e-09,-2.01809e-13,-71257.7,-32.8632], Tmin=(615.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-578.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCFFO) + group(CdCsCdF) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'FC1=CC1F(8567)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {3,S} {4,D} {7,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-116.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,2950,1000,531.912],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94882,0.00290669,7.19373e-05,-1.30708e-07,7.19211e-11,-13976,9.81734], Tmin=(10,'K'), Tmax=(588.625,'K')), NASAPolynomial(coeffs=[2.23287,0.026471,-1.8446e-05,6.01458e-09,-7.38973e-13,-13980.2,15.4347], Tmin=(588.625,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-116.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1(F)C=C1(3904)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {3,S} {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-162.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,289.626,289.697,289.921,290.327],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2895.87,'J/mol'), sigma=(5.13861,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.33 K, Pc=48.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94756,0.00289645,7.1868e-05,-1.26533e-07,6.69252e-11,-19583.6,7.86706], Tmin=(10,'K'), Tmax=(615.385,'K')), NASAPolynomial(coeffs=[2.03451,0.0280179,-2.02892e-05,6.80462e-09,-8.52836e-13,-19588.4,14.2158], Tmin=(615.385,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-162.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1(F)CDC1""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-565.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-411.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-300.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-141.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (49.8931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-213.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-429.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-427.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-415.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-129.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-138.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-123.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-426.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-246.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-194.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-260.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-188.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-406.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-328.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-325.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-314.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-13.3208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-302.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-297.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-333.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    products = ['HF(38)', 'CF2O(49)', 'FC1=CC1(F)F(974)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(59.8647,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'FC=CC(F)(F)OF(354)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.52676e-06,'m^3/(mol*s)'), n=2.97694, Ea=(35.343,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CF2;mb_db] for rate rule [CF2;mb_db_twocdisub]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'FOC(F)(F)C=C(F)F(1384)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(36.1911,'m^3/(mol*s)'), n=1.36816, Ea=(10.0248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;mb_db] for rate rule [carbene;mb_db_trisub]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(43)', 'FOC1C(F)C1(F)F(11812)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(188.831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['FOF(328)', 'F[C]C1C(F)C1(F)F-2(9908)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(13.8048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OF(195)', 'FC(F)=C1C(F)C1(F)F(3033)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C([C](F)F)C(F)(F)OF(3640)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.16053e+11,'s^-1'), n=0.0244333, Ea=(7.84361,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_single;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(F)(F)[CH]C(F)(F)OF(11813)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FOC(F)(F)[CH]C(F)[C](F)F(11814)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'FOC(F)(F)C1[C](F)C1F(11815)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FOC(F)(F)C1[CH]C1(F)F(11816)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'FO[C](F)C1C(F)C1(F)F(11817)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]C(F)(F)C1C(F)C1(F)F(11818)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(127)', 'F[C](F)C1C(F)C1(F)F(9346)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FO[C](F)F(236)', 'FC1[CH]C1(F)F(5918)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.17885e+10,'m^3/(mol*s)'), n=-0.943326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'FOC(F)(F)[C]1C(F)C1(F)F(11819)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'FOC(F)(F)C1[C](F)C1(F)F(11820)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.03232e+06,'m^3/(mol*s)'), n=-0.00929868, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(78)', 'O=C(F)C1C(F)C1(F)F(11821)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(84.4306,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'FOC(F)(F)C1=C(F)C1F(11822)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(166.565,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'FOC(F)(F)C1=CC1(F)F(11781)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(227.633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'FOC(F)=C1C(F)C1(F)F(11823)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(157.51,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F2(78)', 'FOC(F)(F)C1C=C1F(11824)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'FOC(F)(F)C1C(F)=C1F(11825)'],
    products = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(148.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    products = ['F2(78)', 'CF2O(49)', 'FC1=CC1F(8567)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(327.623,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FOC(F)(F)C1C(F)C1(F)F(11638)'],
    products = ['F2(78)', 'CF2O(49)', 'FC1(F)C=C1(3904)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(292.028,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

network(
    label = 'PDepNetwork #3468',
    isomers = [
        'FOC(F)(F)C1C(F)C1(F)F(11638)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(49)', 'FC1=CC1(F)F(974)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3468',
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

