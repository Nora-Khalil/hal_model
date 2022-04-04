species(
    label = '[CH2]C(F)CC(F)(F)OF(9725)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {13,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-635.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,3000,3100,440,815,1455,1000,212.774,222.072,226.353],'cm^-1')),
        HinderedRotor(inertia=(0.00359333,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.431364,'amu*angstrom^2'), symmetry=1, barrier=(14.1272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228077,'amu*angstrom^2'), symmetry=1, barrier=(7.36725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01965,'amu*angstrom^2'), symmetry=1, barrier=(32.4832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0532069,0.0966201,-0.000135871,1.02885e-07,-3.15316e-11,-76305,28.7245], Tmin=(100,'K'), Tmax=(794.547,'K')), NASAPolynomial(coeffs=[12.4143,0.0338533,-1.73725e-05,3.45618e-09,-2.4621e-13,-78286.1,-28.5564], Tmin=(794.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-635.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = '[CH2]C(=C)F(1340)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-49.3492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0540617,'amu*angstrom^2'), symmetry=1, barrier=(25.2827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2755.66,'J/mol'), sigma=(4.80499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=430.43 K, Pc=56.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93971,0.00348861,8.08985e-05,-1.51352e-07,8.61066e-11,-5933.05,7.77268], Tmin=(10,'K'), Tmax=(574.141,'K')), NASAPolynomial(coeffs=[2.51346,0.0273693,-1.79226e-05,5.69534e-09,-6.965e-13,-5999.1,11.8607], Tmin=(574.141,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-49.3492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[CH2]C(DC)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(F)C(F)(F)OF(4806)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-587.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.221478,'amu*angstrom^2'), symmetry=1, barrier=(5.09221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221529,'amu*angstrom^2'), symmetry=1, barrier=(5.0934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973973,'amu*angstrom^2'), symmetry=1, barrier=(22.3936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28654,0.06132,-7.32724e-05,4.55384e-08,-1.14604e-11,-70687.1,14.7537], Tmin=(10,'K'), Tmax=(932.421,'K')), NASAPolynomial(coeffs=[11.5409,0.0259095,-1.63068e-05,4.8088e-09,-5.39998e-13,-72226.4,-24.4912], Tmin=(932.421,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-587.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""[CH2]C(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(F)COF(505)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-154.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1997.95],'cm^-1')),
        HinderedRotor(inertia=(0.394396,'amu*angstrom^2'), symmetry=1, barrier=(9.06794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395458,'amu*angstrom^2'), symmetry=1, barrier=(9.09235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0782,'amu*angstrom^2'), symmetry=1, barrier=(24.7899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.72551,0.0339383,0.000309611,-2.5763e-06,5.70343e-09,-18581.5,11.8009], Tmin=(10,'K'), Tmax=(168.905,'K')), NASAPolynomial(coeffs=[5.20448,0.0351569,-2.30788e-05,7.22866e-09,-8.64338e-13,-18683.2,5.76549], Tmin=(168.905,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-154.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""[CH2]C(F)COF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]CCC(F)(F)OF(9717)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
9  C u1 p0 c0 {3,S} {7,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-626.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,351.038,351.107,351.422,351.463],'cm^-1')),
        HinderedRotor(inertia=(0.0899149,'amu*angstrom^2'), symmetry=1, barrier=(7.84776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.089565,'amu*angstrom^2'), symmetry=1, barrier=(7.84765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895907,'amu*angstrom^2'), symmetry=1, barrier=(7.84673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319023,'amu*angstrom^2'), symmetry=1, barrier=(28.0481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3066.47,'J/mol'), sigma=(5.51817,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.98 K, Pc=41.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0211082,0.0958907,-0.000136596,1.05804e-07,-3.32243e-11,-75244.9,29.3822], Tmin=(100,'K'), Tmax=(775.766,'K')), NASAPolynomial(coeffs=[11.9228,0.0343057,-1.7517e-05,3.47214e-09,-2.46592e-13,-77098.1,-25.2076], Tmin=(775.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-626.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'OF(142)',
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
    label = '[CH2]C(F)C=C(F)F(3097)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-409.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.368265,'amu*angstrom^2'), symmetry=1, barrier=(8.46715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368204,'amu*angstrom^2'), symmetry=1, barrier=(8.46575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.64267,0.0314281,0.000119565,-5.01592e-07,5.30426e-10,-49273.2,13.7197], Tmin=(10,'K'), Tmax=(350.362,'K')), NASAPolynomial(coeffs=[5.75055,0.0362338,-2.4613e-05,7.94239e-09,-9.72017e-13,-49598.1,3.23236], Tmin=(350.362,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-409.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[CH2]C(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]CC(F)(F)OF(4813)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-602.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,223.564,223.652,223.744],'cm^-1')),
        HinderedRotor(inertia=(0.00952725,'amu*angstrom^2'), symmetry=1, barrier=(22.2344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390294,'amu*angstrom^2'), symmetry=1, barrier=(13.8516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19314,'amu*angstrom^2'), symmetry=1, barrier=(42.3752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48408,0.0447428,1.57218e-05,-1.18875e-07,9.14758e-11,-72462,14.1075], Tmin=(10,'K'), Tmax=(552.965,'K')), NASAPolynomial(coeffs=[7.93165,0.036248,-2.54641e-05,8.21622e-09,-9.92547e-13,-73315.9,-7.98774], Tmin=(552.965,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-602.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[CH]CC(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]C(F)CC(F)(F)OF-2(9937)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u2 p0 c0 {7,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-394.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0225925,0.0955412,-0.00013582,1.00833e-07,-2.99886e-11,-47261.8,28.2143], Tmin=(100,'K'), Tmax=(820.833,'K')), NASAPolynomial(coeffs=[13.4902,0.0296918,-1.54854e-05,3.09891e-09,-2.21495e-13,-49480.1,-34.3089], Tmin=(820.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
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
    label = 'C=CCC(F)(F)OF(9938)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-502.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,890.007],'cm^-1')),
        HinderedRotor(inertia=(0.709275,'amu*angstrom^2'), symmetry=1, barrier=(16.3076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706152,'amu*angstrom^2'), symmetry=1, barrier=(16.2358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707601,'amu*angstrom^2'), symmetry=1, barrier=(16.2691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909224,0.0695373,-6.98337e-05,3.62057e-08,-7.63244e-12,-60297.6,24.4295], Tmin=(100,'K'), Tmax=(1130.72,'K')), NASAPolynomial(coeffs=[13.0947,0.0264306,-1.2649e-05,2.49001e-09,-1.77979e-13,-63053.3,-35.8552], Tmin=(1130.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-502.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(F)(F)OF(4811)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {5,S} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-374.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.764002,'amu*angstrom^2'), symmetry=1, barrier=(17.5659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7645,'amu*angstrom^2'), symmetry=1, barrier=(17.5774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78203,0.0146305,0.000135739,-3.80847e-07,2.91218e-10,-45081.4,11.4384], Tmin=(10,'K'), Tmax=(479.598,'K')), NASAPolynomial(coeffs=[7.42144,0.0240799,-1.8304e-05,6.32762e-09,-8.08268e-13,-45888.3,-8.21775], Tmin=(479.598,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-374.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""[CH2]C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(F)CC(F)(F)OF(9939)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-710.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,323,467,575,827,1418,2950,3100,1380,975,1025,1650,254.079,258.304,261.643],'cm^-1')),
        HinderedRotor(inertia=(0.371767,'amu*angstrom^2'), symmetry=1, barrier=(18.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109684,'amu*angstrom^2'), symmetry=1, barrier=(5.46122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291209,'amu*angstrom^2'), symmetry=1, barrier=(14.4739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448632,0.0831427,-0.000101817,6.50322e-08,-1.67712e-11,-85360.5,26.1926], Tmin=(100,'K'), Tmax=(938.285,'K')), NASAPolynomial(coeffs=[13.2859,0.0284142,-1.43208e-05,2.86254e-09,-2.05825e-13,-87769.4,-34.9214], Tmin=(938.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-710.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(CdCsCdF) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH]CC(F)(F)OF(9940)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u1 p0 c0 {5,S} {8,S} {11,S}
8  C u1 p0 c0 {7,S} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-230.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00881128,'amu*angstrom^2'), symmetry=1, barrier=(2.71556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119486,'amu*angstrom^2'), symmetry=1, barrier=(2.74723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118784,'amu*angstrom^2'), symmetry=1, barrier=(2.73109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00882581,'amu*angstrom^2'), symmetry=1, barrier=(2.70528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621552,0.0809153,-0.000119333,1.01269e-07,-3.49316e-11,-27589,29.5048], Tmin=(100,'K'), Tmax=(781.256,'K')), NASAPolynomial(coeffs=[8.58315,0.0337935,-1.66515e-05,3.22937e-09,-2.25547e-13,-28639,-5.69805], Tmin=(781.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFO) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(F)C[C](F)OF(9941)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-206.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,395,473,707,1436,3000,3100,440,815,1455,1000,221.201,221.207,221.233],'cm^-1')),
        HinderedRotor(inertia=(0.179769,'amu*angstrom^2'), symmetry=1, barrier=(6.24227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179754,'amu*angstrom^2'), symmetry=1, barrier=(6.2427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958839,'amu*angstrom^2'), symmetry=1, barrier=(33.258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179683,'amu*angstrom^2'), symmetry=1, barrier=(6.24259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125212,0.0940968,-0.000153368,1.35151e-07,-4.69672e-11,-24685,28.6605], Tmin=(100,'K'), Tmax=(811.999,'K')), NASAPolynomial(coeffs=[9.73788,0.0335614,-1.71904e-05,3.35304e-09,-2.33486e-13,-25811.5,-13.0372], Tmin=(811.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFHO) + group(Cs-CsHHH) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2]C(F)CC([O])(F)F(9942)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-510.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,351,323,533,609,664,892,1120,1201,3000,3100,440,815,1455,1000,224.326,224.337],'cm^-1')),
        HinderedRotor(inertia=(0.166214,'amu*angstrom^2'), symmetry=1, barrier=(5.93555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166218,'amu*angstrom^2'), symmetry=1, barrier=(5.93538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166201,'amu*angstrom^2'), symmetry=1, barrier=(5.93555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613937,0.0805976,-0.000109346,8.24847e-08,-2.54509e-11,-61238.4,26.7445], Tmin=(100,'K'), Tmax=(786.796,'K')), NASAPolynomial(coeffs=[10.2818,0.0314446,-1.56332e-05,3.07596e-09,-2.17966e-13,-62759.7,-17.5787], Tmin=(786.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = '[CH2]C(F)C[C](F)F(3047)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-363.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.19429,'amu*angstrom^2'), symmetry=1, barrier=(4.46711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194081,'amu*angstrom^2'), symmetry=1, barrier=(4.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513725,'amu*angstrom^2'), symmetry=1, barrier=(11.8116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2934.83,'J/mol'), sigma=(5.15676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.41 K, Pc=48.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773699,0.0797718,-0.000131731,1.20612e-07,-4.33047e-11,-43654.2,24.5844], Tmin=(100,'K'), Tmax=(816.706,'K')), NASAPolynomial(coeffs=[7.21739,0.032635,-1.65473e-05,3.23486e-09,-2.25695e-13,-44187.2,-2.01748], Tmin=(816.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = 'FO[C](F)F(4688)',
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
    label = '[CH2]C([CH2])F(238)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (84.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.0782826,'amu*angstrom^2'), symmetry=1, barrier=(1.79987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0782693,'amu*angstrom^2'), symmetry=1, barrier=(1.79956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68793,0.0314085,-2.75226e-05,1.56508e-08,-4.08194e-12,10246.7,15.5107], Tmin=(100,'K'), Tmax=(875.654,'K')), NASAPolynomial(coeffs=[5.01308,0.0207872,-9.32841e-06,1.79894e-09,-1.2725e-13,9839.46,4.60193], Tmin=(875.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(H)(H)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2]C(F)[CH]C(F)(F)OF(9943)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u1 p0 c0 {6,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-435.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,269.882,269.883,2313.83],'cm^-1')),
        HinderedRotor(inertia=(0.170205,'amu*angstrom^2'), symmetry=1, barrier=(8.7973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170205,'amu*angstrom^2'), symmetry=1, barrier=(8.7973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524634,'amu*angstrom^2'), symmetry=1, barrier=(27.1163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52463,'amu*angstrom^2'), symmetry=1, barrier=(27.1163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162409,0.0911332,-0.000126147,9.16874e-08,-2.6788e-11,-52269.3,30.4486], Tmin=(100,'K'), Tmax=(833.963,'K')), NASAPolynomial(coeffs=[13.0134,0.029494,-1.52782e-05,3.05819e-09,-2.18927e-13,-54412.7,-29.2166], Tmin=(833.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(CCJCO) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2][C](F)CC(F)(F)OF(9944)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {9,S}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-449.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,212,367,445,1450,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0808398,'amu*angstrom^2'), symmetry=1, barrier=(1.85867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0776766,'amu*angstrom^2'), symmetry=1, barrier=(1.78594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0805481,'amu*angstrom^2'), symmetry=1, barrier=(1.85196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892969,'amu*angstrom^2'), symmetry=1, barrier=(20.5311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.215942,0.103465,-0.000176957,1.58762e-07,-5.5375e-11,-53865.3,31.009], Tmin=(100,'K'), Tmax=(826.484,'K')), NASAPolynomial(coeffs=[10.2586,0.0343199,-1.79777e-05,3.51399e-09,-2.44023e-13,-54966.6,-13.7163], Tmin=(826.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-449.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(Csj(Cs-CsHH)(Cs-HHH)(F1s)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = '[CH2]C(F)CC(=O)F(1124)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-486.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.105794,'amu*angstrom^2'), symmetry=1, barrier=(2.43242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105895,'amu*angstrom^2'), symmetry=1, barrier=(2.43473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461169,'amu*angstrom^2'), symmetry=1, barrier=(10.6032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98244,0.0546943,-2.58538e-05,-7.64831e-08,9.56173e-11,-58468.6,21.3024], Tmin=(100,'K'), Tmax=(453.228,'K')), NASAPolynomial(coeffs=[6.70548,0.0323537,-1.5933e-05,3.08872e-09,-2.15478e-13,-59095.4,0.0623724], Tmin=(453.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-486.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(COCsFO) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'C=C[CH]C(F)(F)OF(9211)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u1 p0 c0 {5,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-385.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42474,'amu*angstrom^2'), symmetry=1, barrier=(32.7576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.429,'amu*angstrom^2'), symmetry=1, barrier=(32.8555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42898,'amu*angstrom^2'), symmetry=1, barrier=(32.8551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64874,0.0739924,-7.93384e-05,4.2323e-08,-8.99017e-12,-46225.3,22.5632], Tmin=(100,'K'), Tmax=(1136.21,'K')), NASAPolynomial(coeffs=[15.3703,0.0221655,-1.09174e-05,2.17721e-09,-1.56884e-13,-49570.6,-50.3397], Tmin=(1136.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C(F)C=C(F)OF(9945)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {6,D}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,174,267,591,721,1107,1278,1348,3273,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.338103,'amu*angstrom^2'), symmetry=1, barrier=(7.77366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70935,'amu*angstrom^2'), symmetry=1, barrier=(39.3014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227445,'amu*angstrom^2'), symmetry=1, barrier=(7.75789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719305,0.0808879,-0.00013137,1.20642e-07,-4.40639e-11,-29477.3,27.0476], Tmin=(100,'K'), Tmax=(794.784,'K')), NASAPolynomial(coeffs=[7.01309,0.0351382,-1.84644e-05,3.65547e-09,-2.57287e-13,-30033.2,0.92588], Tmin=(794.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFO) + radical(Csj(Cs-F1sCdH)(H)(H))"""),
)

species(
    label = '[CH2]C=C(1339)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29613,0.00579239,4.33915e-05,-5.99882e-08,2.33811e-11,18908.2,9.01999], Tmin=(100,'K'), Tmax=(942.184,'K')), NASAPolynomial(coeffs=[8.06867,0.0101836,-2.8479e-06,5.0087e-10,-3.7962e-14,16914.7,-19.5274], Tmin=(942.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[C](F)CC(F)(F)OF(9866)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {3,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-655.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,219.722,219.825,2330.54],'cm^-1')),
        HinderedRotor(inertia=(0.281674,'amu*angstrom^2'), symmetry=1, barrier=(9.66106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281637,'amu*angstrom^2'), symmetry=1, barrier=(9.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281632,'amu*angstrom^2'), symmetry=1, barrier=(9.66095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.748341,'amu*angstrom^2'), symmetry=1, barrier=(25.6629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.104045,0.0989899,-0.000156252,1.35789e-07,-4.71161e-11,-78731.4,29.2035], Tmin=(100,'K'), Tmax=(792.823,'K')), NASAPolynomial(coeffs=[10.1949,0.0366672,-1.87351e-05,3.66889e-09,-2.56799e-13,-80038.8,-16.0382], Tmin=(792.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-655.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(Csj(Cs-CsHH)(Cs-HHH)(F1s))"""),
)

species(
    label = 'CC(F)[CH]C(F)(F)OF(9881)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
8  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
9  C u1 p0 c0 {6,S} {7,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-642.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.351129,0.0857604,-0.000102286,6.44164e-08,-1.65041e-11,-77138.7,28.3667], Tmin=(100,'K'), Tmax=(940.757,'K')), NASAPolynomial(coeffs=[13.1,0.0315537,-1.58555e-05,3.16793e-09,-2.27789e-13,-79537.4,-32.3606], Tmin=(940.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-642.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsHHH) + radical(CCJCO)"""),
)

species(
    label = '[O]C(F)(F)CC(F)CF(9946)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-891.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.417474,0.0848943,-0.00010366,6.86431e-08,-1.8636e-11,-107042,27.3597], Tmin=(100,'K'), Tmax=(887.702,'K')), NASAPolynomial(coeffs=[11.8427,0.033412,-1.66673e-05,3.31191e-09,-2.37113e-13,-109070,-26.3997], Tmin=(887.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-891.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'FC[CH]CC(F)(F)OF(9904)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-624.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,354.453,354.453,354.453,354.453],'cm^-1')),
        HinderedRotor(inertia=(0.0799197,'amu*angstrom^2'), symmetry=1, barrier=(7.12524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799199,'amu*angstrom^2'), symmetry=1, barrier=(7.12524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200384,'amu*angstrom^2'), symmetry=1, barrier=(17.8652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375338,'amu*angstrom^2'), symmetry=1, barrier=(33.4632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291458,0.0884455,-0.000112123,7.70032e-08,-2.1662e-11,-74947.9,29.4172], Tmin=(100,'K'), Tmax=(858.229,'K')), NASAPolynomial(coeffs=[11.908,0.0343046,-1.74983e-05,3.50071e-09,-2.51269e-13,-76941.9,-24.8501], Tmin=(858.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-624.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFHH) + radical(Csj(Cs-CsHH)(Cs-F1sHH)(H))"""),
)

species(
    label = 'FCC(F)C[C](F)OF(9947)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
9  C u1 p0 c0 {3,S} {5,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-587.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,250,417,511,1155,1315,1456,3119,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,395,473,707,1436,238.1,238.105,238.126,238.156],'cm^-1')),
        HinderedRotor(inertia=(0.241935,'amu*angstrom^2'), symmetry=1, barrier=(9.74123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242083,'amu*angstrom^2'), symmetry=1, barrier=(9.74119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40491,'amu*angstrom^2'), symmetry=1, barrier=(16.2994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845864,'amu*angstrom^2'), symmetry=1, barrier=(34.0416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.166707,0.0995726,-0.000152127,1.27753e-07,-4.33314e-11,-70484.5,29.6151], Tmin=(100,'K'), Tmax=(756.07,'K')), NASAPolynomial(coeffs=[11.2609,0.0356111,-1.828e-05,3.60362e-09,-2.53946e-13,-72112.3,-21.6588], Tmin=(756.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsF1sO2s)"""),
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
    E0 = (-357.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (226.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (86.0808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-132.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-24.0134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (27.9563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (66.8065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-122.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-260.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-209.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (91.6551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (115.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-183.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-12.0602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-9.59546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (16.4147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (25.2245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (12.5634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-163.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-138.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-103.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-151.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-222.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-239.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-292.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-226.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-119.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    products = ['HF(38)', 'CF2O(49)', '[CH2]C(=C)F(1340)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(28.5431,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[CH2]C(F)C(F)(F)OF(4806)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(146.212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(43)', '[CH2]C(F)COF(505)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(195.188,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]CCC(F)(F)OF(9717)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(142)', '[CH2]C(F)C=C(F)F(3097)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(18)', 'F[CH]CC(F)(F)OF(4813)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', '[CH]C(F)CC(F)(F)OF-2(9937)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'C=CCC(F)(F)OF(9938)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(57.5373,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)OF(4811)', 'CH2CHF(56)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(18.2609,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C=C(F)CC(F)(F)OF(9939)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.08158e+21,'m^3/(mol*s)'), n=-4.30708, Ea=(40.618,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS_Ext-6R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2][CH]CC(F)(F)OF(9940)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[CH2]C(F)C[C](F)OF(9941)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -42.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[CH2]C(F)CC([O])(F)F(9942)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(4.20435,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(126)', '[CH2]C(F)C[C](F)F(3047)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(F)(F)OF(4811)', '[CH2][CH]F(233)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -10.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['FO[C](F)F(4688)', '[CH2]C([CH2])F(238)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -17.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[CH2]C(F)[CH]C(F)(F)OF(9943)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[CH2][C](F)CC(F)(F)OF(9944)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.68633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F2(78)', '[CH2]C(F)CC(=O)F(1124)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(82.446,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'C=C[CH]C(F)(F)OF(9211)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(278.537,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[CH2]C(F)C=C(F)OF(9945)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(174.34,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    products = ['F2(78)', 'CF2O(49)', '[CH2]C=C(1339)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(235.006,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[C](F)CC(F)(F)OF(9866)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.465e+11,'s^-1'), n=0, Ea=(184.64,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 440 used for R2H_S;C_rad_out_noH;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    products = ['CC(F)[CH]C(F)(F)OF(9881)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    products = ['[O]C(F)(F)CC(F)CF(9946)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(93.8362,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC[CH]CC(F)(F)OF(9904)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(148.448,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['FCC(F)C[C](F)OF(9947)'],
    products = ['[CH2]C(F)CC(F)(F)OF(9725)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(218.536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #3461',
    isomers = [
        '[CH2]C(F)CC(F)(F)OF(9725)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(49)', '[CH2]C(=C)F(1340)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3461',
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

