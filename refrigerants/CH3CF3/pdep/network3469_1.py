species(
    label = 'FOC(F)(F)C1(F)CC1(F)F(11639)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-1030.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,321,327,383,445,494,1500,222,329,445,522,589,1214,1475,2750,3150,900,1100,223,363,546,575,694,1179,1410,278.654,705.711,1231.14,1231.28,1231.3],'cm^-1')),
        HinderedRotor(inertia=(0.115717,'amu*angstrom^2'), symmetry=1, barrier=(6.37028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804358,'amu*angstrom^2'), symmetry=1, barrier=(44.292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459333,0.0923653,-0.000113615,7.12331e-08,-1.7924e-11,-123764,28.7821], Tmin=(100,'K'), Tmax=(962.556,'K')), NASAPolynomial(coeffs=[15.2405,0.0292222,-1.52158e-05,3.08116e-09,-2.23164e-13,-126689,-43.943], Tmin=(962.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1030.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs)"""),
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
    label = 'C=C(F)C(F)(F)OF(339)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-650.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.493364,'amu*angstrom^2'), symmetry=1, barrier=(11.3434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815502,'amu*angstrom^2'), symmetry=1, barrier=(18.75,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.68127,0.0224957,0.000176709,-5.29719e-07,4.28239e-10,-78282.4,12.8424], Tmin=(10,'K'), Tmax=(453.826,'K')), NASAPolynomial(coeffs=[7.96789,0.0342935,-2.61585e-05,9.0192e-09,-1.14616e-12,-79182,-10.0767], Tmin=(453.826,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-650.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)C(F)=C(F)F(2539)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
"""),
    E0 = (-997.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,323,467,575,827,1418,182,240,577,636,1210,1413,235.584,236.82],'cm^-1')),
        HinderedRotor(inertia=(0.218377,'amu*angstrom^2'), symmetry=1, barrier=(8.63839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88098,'amu*angstrom^2'), symmetry=1, barrier=(35.3275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (166.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45466,0.0456506,6.7432e-05,-2.76879e-07,2.16601e-10,-119995,14.3265], Tmin=(10,'K'), Tmax=(516.918,'K')), NASAPolynomial(coeffs=[10.5448,0.0345863,-2.75626e-05,9.55645e-09,-1.20812e-12,-121313,-20.8611], Tmin=(516.918,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-997.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""FOC(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC1(F)CC1(F)F(11767)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-567.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,2750,3150,900,1100,222,329,445,522,589,1214,1475,440.436,441.738,442.537,445.923],'cm^-1')),
        HinderedRotor(inertia=(0.000861747,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44765,0.0583243,-5.90096e-05,3.00187e-08,-6.20263e-12,-68112.1,19.3035], Tmin=(100,'K'), Tmax=(1149.28,'K')), NASAPolynomial(coeffs=[11.984,0.0216531,-1.11479e-05,2.2554e-09,-1.63376e-13,-70533.9,-32.9944], Tmin=(1149.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-567.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs(O2))"""),
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
    label = 'F[C]C1(F)CC1(F)F(11768)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p1 c0 {4,S} {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-395.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,222,329,445,522,589,1214,1475,617,898,1187,180,180,1456.41],'cm^-1')),
        HinderedRotor(inertia=(0.0693165,'amu*angstrom^2'), symmetry=1, barrier=(25.4266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55033,0.0550789,-5.61053e-05,2.90204e-08,-6.07099e-12,-47516.4,21.6455], Tmin=(100,'K'), Tmax=(1141.42,'K')), NASAPolynomial(coeffs=[11.6184,0.0197962,-9.73822e-06,1.93877e-09,-1.39404e-13,-49814.8,-28.2589], Tmin=(1141.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-395.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(CJ2_singlet-FCs) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FOC1(C(F)(F)F)CC1(F)F(11769)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-1065.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.248095,0.0914372,-0.000105275,5.91151e-08,-1.29681e-11,-127959,29.026], Tmin=(100,'K'), Tmax=(1115.11,'K')), NASAPolynomial(coeffs=[19.5965,0.0202531,-9.52102e-06,1.86918e-09,-1.34025e-13,-132385,-68.8748], Tmin=(1115.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1065.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsOs) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + ring(Cs(F)(F)-Cs(C)-Cs)"""),
)

species(
    label = 'FC(F)=C1CC1(F)F(11770)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-594.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,274,345,380,539,705,1166,1213,182,240,577,636,1210,1413,451.098,451.109,451.166,1205.01,1577.4,1577.73,1577.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69374,0.0280997,4.0777e-05,-1.11879e-07,6.68608e-11,-71497.5,12.261], Tmin=(10,'K'), Tmax=(629.181,'K')), NASAPolynomial(coeffs=[5.86687,0.0346011,-2.31592e-05,7.18905e-09,-8.4254e-13,-72173.1,-0.41197], Tmin=(629.181,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-594.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC(F)DC1CC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)([C](F)F)C(F)(F)OF(11771)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u1 p0 c0 {8,S} {12,S} {13,S}
12 H u0 p0 c0 {11,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-829.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,223,363,546,575,694,1179,1410,190,488,555,1236,1407,3000,3100,440,815,1455,1000,254.588,254.594,254.604],'cm^-1')),
        HinderedRotor(inertia=(0.154984,'amu*angstrom^2'), symmetry=1, barrier=(7.13125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155103,'amu*angstrom^2'), symmetry=1, barrier=(7.13147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155152,'amu*angstrom^2'), symmetry=1, barrier=(7.13128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836779,'amu*angstrom^2'), symmetry=1, barrier=(38.4693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26641,0.129018,-0.000232111,2.06951e-07,-7.08627e-11,-99619.8,35.1175], Tmin=(100,'K'), Tmax=(834.657,'K')), NASAPolynomial(coeffs=[14.1892,0.0324048,-1.79691e-05,3.54903e-09,-2.4642e-13,-101415,-31.9495], Tmin=(834.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-829.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + radical(Csj(Cs-F1sCsCs)(F1s)(F1s)) + radical(Csj(Cs-F1sCsCs)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)[C](F)C(F)(F)OF(11772)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
10 C u1 p0 c0 {5,S} {8,S} {9,S}
11 C u1 p0 c0 {8,S} {12,S} {13,S}
12 H u0 p0 c0 {11,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-844.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,253,525,597,667,842,1178,1324,212,367,445,1450,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.236457,'amu*angstrom^2'), symmetry=1, barrier=(5.43662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242289,'amu*angstrom^2'), symmetry=1, barrier=(5.5707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615255,'amu*angstrom^2'), symmetry=1, barrier=(14.1459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64615,'amu*angstrom^2'), symmetry=1, barrier=(37.8482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25506,0.126997,-0.00022208,1.93664e-07,-6.54605e-11,-101361,35.8493], Tmin=(100,'K'), Tmax=(818.48,'K')), NASAPolynomial(coeffs=[15.2674,0.0306734,-1.70049e-05,3.37499e-09,-2.35713e-13,-103544,-37.3653], Tmin=(818.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-844.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'FOC(F)(F)[C](F)C[C](F)F(3657)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {10,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
10 C u1 p0 c0 {3,S} {8,S} {9,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-841.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,253,525,597,667,842,1178,1324,212,367,445,1450,190,488,555,1236,1407,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.203107,'amu*angstrom^2'), symmetry=1, barrier=(4.66982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200208,'amu*angstrom^2'), symmetry=1, barrier=(4.60318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20175,'amu*angstrom^2'), symmetry=1, barrier=(4.63862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22203,'amu*angstrom^2'), symmetry=1, barrier=(28.0968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.990406,0.124202,-0.000226805,2.07865e-07,-7.29562e-11,-101074,36.2264], Tmin=(100,'K'), Tmax=(834.829,'K')), NASAPolynomial(coeffs=[11.7304,0.0362023,-2.00878e-05,3.97595e-09,-2.766e-13,-102256,-17.2022], Tmin=(834.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-841.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + radical(CsCsCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'FOC(F)(F)[C]1CC1(F)F(11773)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-706.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,215,315,519,588,595,1205,1248,253,525,597,667,842,1178,1324,525.086,525.09,525.091,525.095,525.107,525.127,525.144,525.154],'cm^-1')),
        HinderedRotor(inertia=(0.000611385,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000611416,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570995,0.0721568,-6.95309e-05,3.25673e-08,-6.04721e-12,-84830.2,27.3477], Tmin=(100,'K'), Tmax=(1292.81,'K')), NASAPolynomial(coeffs=[17.0517,0.0211657,-1.03686e-05,2.05934e-09,-1.47733e-13,-89091.5,-56.3948], Tmin=(1292.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(CsCFFO) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CCJ(C)CO)"""),
)

species(
    label = 'FOC(F)(F)C1(F)C[C]1F(11774)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-605.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,2750,3150,900,1100,223,363,546,575,694,1179,1410,212,367,445,1450,520.449,521.684,521.773,522.33,522.499],'cm^-1')),
        HinderedRotor(inertia=(0.00061639,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000619208,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.45657,0.0831873,-0.000102039,6.4096e-08,-1.62133e-11,-72677.1,27.3242], Tmin=(100,'K'), Tmax=(956.037,'K')), NASAPolynomial(coeffs=[13.8256,0.0272519,-1.42771e-05,2.89743e-09,-2.10064e-13,-75233.3,-36.5726], Tmin=(956.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-605.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-F1sCsCs)(Cs-CsHH)(F1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FO[C](F)C1(F)CC1(F)F(11775)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u1 p0 c0 {4,S} {6,S} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-613.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,2750,3150,900,1100,222,329,445,522,589,1214,1475,395,473,707,1436,180,180,939.407,941.605,942.938],'cm^-1')),
        HinderedRotor(inertia=(0.0864317,'amu*angstrom^2'), symmetry=1, barrier=(1.98724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0899607,'amu*angstrom^2'), symmetry=1, barrier=(2.06837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193908,0.0911353,-0.000135265,1.07595e-07,-3.45507e-11,-73653.9,28.5565], Tmin=(100,'K'), Tmax=(759.902,'K')), NASAPolynomial(coeffs=[11.6412,0.0308843,-1.63441e-05,3.27429e-09,-2.33506e-13,-75393.9,-23.528], Tmin=(759.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)C1(F)CC1(F)F(11776)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-904.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([321,327,383,445,494,1500,2750,3150,900,1100,222,329,445,522,589,1214,1475,351,323,533,609,664,892,1120,1201,782.796,1425.93,1426.49,1427.03],'cm^-1')),
        HinderedRotor(inertia=(0.281168,'amu*angstrom^2'), symmetry=1, barrier=(6.46459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691813,0.0766242,-8.82114e-05,5.24654e-08,-1.26134e-11,-108697,26.8763], Tmin=(100,'K'), Tmax=(1002.2,'K')), NASAPolynomial(coeffs=[13.2229,0.0266103,-1.3356e-05,2.67189e-09,-1.92474e-13,-111209,-33.6069], Tmin=(1002.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-904.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(O2sj(Cs-CsF1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[C](F)C1(F)CC1(F)F(11777)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-758.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,222,329,445,522,589,1214,1475,190,488,555,1236,1407,821.443,821.452,821.467,821.468],'cm^-1')),
        HinderedRotor(inertia=(0.00892612,'amu*angstrom^2'), symmetry=1, barrier=(4.27423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53588,0.0427112,6.49845e-06,-6.80289e-08,4.40318e-11,-91192.3,14.5094], Tmin=(10,'K'), Tmax=(670.724,'K')), NASAPolynomial(coeffs=[8.00598,0.0362685,-2.43032e-05,7.52294e-09,-8.78078e-13,-92246.6,-8.66065], Tmin=(670.724,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-758.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), label="""F[C](F)C1(F)CC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]1CC1(F)F(7814)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-292.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,3150,900,1100,212,367,445,1450,1382.51,1382.6,1382.67],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0431,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89779,0.00619863,0.000105234,-2.23257e-07,1.39691e-10,-35186.4,11.7365], Tmin=(10,'K'), Tmax=(543.429,'K')), NASAPolynomial(coeffs=[3.74321,0.0296038,-2.08333e-05,6.80265e-09,-8.34139e-13,-35498.4,9.36284], Tmin=(543.429,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-292.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""F[C]1CC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)C1(F)[CH]C1(F)F(11778)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
11 C u1 p0 c0 {8,S} {9,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-795.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,2950,1000,180,180,990.731,992.148],'cm^-1')),
        HinderedRotor(inertia=(2.14596,'amu*angstrom^2'), symmetry=1, barrier=(49.3398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150078,'amu*angstrom^2'), symmetry=1, barrier=(3.45058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.177335,0.0990486,-0.000140077,1.00209e-07,-2.85058e-11,-95553.5,29.9922], Tmin=(100,'K'), Tmax=(858.843,'K')), NASAPolynomial(coeffs=[15.2132,0.0273689,-1.48868e-05,3.03279e-09,-2.19184e-13,-98197.2,-41.9166], Tmin=(858.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-795.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'O=C(F)C1(F)CC1(F)F(8238)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-908.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,222,329,445,522,589,1214,1475,486,617,768,1157,1926,497.345,498.209,498.255,498.539],'cm^-1')),
        HinderedRotor(inertia=(0.000682988,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38496,0.0587656,-5.52928e-05,2.61018e-08,-5.01419e-12,-109169,23.8239], Tmin=(100,'K'), Tmax=(1227.54,'K')), NASAPolynomial(coeffs=[12.3954,0.0228872,-1.14504e-05,2.29101e-09,-1.64839e-13,-111872,-31.5521], Tmin=(1227.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-908.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(COCsFO) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FOC(F)(F)C1=C(F)C1(11779)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-435.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,275,321,533,585,746,850,1103,323,467,575,827,1418,180,180,397.752,1305.66,1305.82,1306.12,3788.66],'cm^-1')),
        HinderedRotor(inertia=(1.71434,'amu*angstrom^2'), symmetry=1, barrier=(39.4162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71439,'amu*angstrom^2'), symmetry=1, barrier=(39.4173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905137,0.0703278,-7.99325e-05,4.58561e-08,-1.05273e-11,-52246.2,25.1031], Tmin=(100,'K'), Tmax=(1052.46,'K')), NASAPolynomial(coeffs=[13.6254,0.0219825,-1.10286e-05,2.20955e-09,-1.59453e-13,-54923.7,-36.9152], Tmin=(1052.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cs-Cd(C)-Cd(F))"""),
)

species(
    label = 'FOC(F)=C1CC1(F)F(11780)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-436.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,3150,900,1100,274,345,380,539,705,1166,1213,326,540,652,719,1357,180,180,257.648,663.717,666.563,1423.98,2042.41,2045.7],'cm^-1')),
        HinderedRotor(inertia=(0.0139838,'amu*angstrom^2'), symmetry=1, barrier=(41.3063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16813,0.0665296,-7.62782e-05,4.63951e-08,-1.15523e-11,-52430.5,25.1345], Tmin=(100,'K'), Tmax=(963.55,'K')), NASAPolynomial(coeffs=[10.9652,0.0258584,-1.29628e-05,2.58748e-09,-1.86016e-13,-54318.5,-21.7668], Tmin=(963.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFO) + ring(Cs-Cs(F)(F)-Cd(Cd))"""),
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
    label = 'FOC(F)(F)C1(F)C=C1F(11782)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {7,S} {10,D} {11,S}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-609.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,2950,1000,323,467,575,827,1418,180,180,540.335],'cm^-1')),
        HinderedRotor(inertia=(0.560953,'amu*angstrom^2'), symmetry=1, barrier=(12.8974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43191,'amu*angstrom^2'), symmetry=1, barrier=(32.9224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184185,0.0915039,-0.000140857,1.11572e-07,-3.52418e-11,-73172.6,25.9884], Tmin=(100,'K'), Tmax=(775.123,'K')), NASAPolynomial(coeffs=[12.8256,0.0262643,-1.45993e-05,2.97388e-09,-2.13787e-13,-75132.2,-31.7784], Tmin=(775.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-609.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'FC1=C(F)C1(8237)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-53.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,260,386,409,525,515,635,761,893,1354,1482,3407.19],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94869,0.00307939,7.47805e-05,-1.45655e-07,8.68807e-11,-6470.71,8.67605], Tmin=(10,'K'), Tmax=(539.608,'K')), NASAPolynomial(coeffs=[2.5013,0.024875,-1.65693e-05,5.21062e-09,-6.23012e-13,-6475.62,13.2731], Tmin=(539.608,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-53.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DC(F)C1""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-565.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-412.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-165.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-161.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (41.4357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-324.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (59.7465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-415.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-429.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-427.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-227.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-125.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-134.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-425.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-249.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-203.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-177.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-418.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-37.6953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-39.1345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-325.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-323.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-249.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    products = ['HF(38)', 'CF2O(49)', 'FC1=CC1(F)F(974)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(58.5091,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'C=C(F)C(F)(F)OF(339)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.52676e-06,'m^3/(mol*s)'), n=2.97694, Ea=(35.343,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CF2;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', 'FOC(F)(F)C(F)=C(F)F(2539)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.56622e+08,'m^3/(mol*s)'), n=-0.441667, Ea=(7.08269,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CH2;mb_db] for rate rule [CH2;mb_db_tetrasub]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(43)', 'FOC1(F)CC1(F)F(11767)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(202.745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['FOF(328)', 'F[C]C1(F)CC1(F)F(11768)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(14.7883,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    products = ['FOC1(C(F)(F)F)CC1(F)F(11769)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(328)', 'FC(F)=C1CC1(F)F(11770)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(F)([C](F)F)C(F)(F)OF(11771)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.16053e+11,'s^-1'), n=0.0244333, Ea=(7.84361,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_noH;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)[C](F)C(F)(F)OF(11772)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.16053e+11,'s^-1'), n=0.0244333, Ea=(7.84361,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_noH;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FOC(F)(F)[C](F)C[C](F)F(3657)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.16053e+11,'s^-1'), n=0.0244333, Ea=(7.84361,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_single;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FOC(F)(F)[C]1CC1(F)F(11773)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'FOC(F)(F)C1(F)C[C]1F(11774)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'FO[C](F)C1(F)CC1(F)F(11775)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]C(F)(F)C1(F)CC1(F)F(11776)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]F(127)', 'F[C](F)C1(F)CC1(F)F(11777)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['FO[C](F)F(236)', 'F[C]1CC1(F)F(7814)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.17885e+10,'m^3/(mol*s)'), n=-0.943326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'FOC(F)(F)C1(F)[CH]C1(F)F(11778)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(78)', 'O=C(F)C1(F)CC1(F)F(8238)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(92.665,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F2(78)', 'FOC(F)(F)C1=C(F)C1(11779)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', 'FOC(F)=C1CC1(F)F(11780)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'FOC(F)(F)C1=CC1(F)F(11781)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(229.61,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'FOC(F)(F)C1(F)C=C1F(11782)'],
    products = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(160.507,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FOC(F)(F)C1(F)CC1(F)F(11639)'],
    products = ['F2(78)', 'CF2O(49)', 'FC1=C(F)C1(8237)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(374.676,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3469',
    isomers = [
        'FOC(F)(F)C1(F)CC1(F)F(11639)',
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
    label = 'PDepNetwork #3469',
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

