species(
    label = '[O]C(F)C(F)(F)C(O)(F)F(6593)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1292.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,222,329,445,522,589,1214,1475,223,363,546,575,694,1179,1410,391,562,707,872,1109,1210,1289,3137,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.293442,'amu*angstrom^2'), symmetry=1, barrier=(6.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565772,'amu*angstrom^2'), symmetry=1, barrier=(13.0082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291876,'amu*angstrom^2'), symmetry=1, barrier=(6.71081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204143,0.0991034,-0.000154567,1.2224e-07,-3.80524e-11,-155272,30.3008], Tmin=(100,'K'), Tmax=(790.535,'K')), NASAPolynomial(coeffs=[14.5859,0.0242631,-1.25516e-05,2.46897e-09,-1.73444e-13,-157610,-37.5751], Tmin=(790.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1292.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH))"""),
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
    label = 'O=C[C](F)F(234)',
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
    label = '[O]C(F)C(O)(F)F(5575)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-868.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,223,363,546,575,694,1179,1410,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.640098,'amu*angstrom^2'), symmetry=1, barrier=(14.7171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32096,'amu*angstrom^2'), symmetry=1, barrier=(30.3715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77187,0.0149771,0.00015394,-4.10068e-07,3.02187e-10,-104429,11.981], Tmin=(10,'K'), Tmax=(493.108,'K')), NASAPolynomial(coeffs=[7.50234,0.0286239,-2.11371e-05,7.20407e-09,-9.14616e-13,-105331,-8.79165], Tmin=(493.108,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-868.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""[O]C(F)C(O)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)F(364)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-426.241,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(66.9995,'amu')),
        NonlinearRotor(inertia=([47.2408,48.1768,88.2639],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([480.226,502.831,615.399,942.972,1118.44,1159.86,1274.53,1324.03,2842.12],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)[C]F(2989)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-539.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,617,898,1187],'cm^-1')),
        HinderedRotor(inertia=(1.39023,'amu*angstrom^2'), symmetry=1, barrier=(31.9641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38775,'amu*angstrom^2'), symmetry=1, barrier=(31.907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03035,0.046278,-6.69891e-05,4.89417e-08,-1.41202e-11,-64798.6,17.5213], Tmin=(100,'K'), Tmax=(850.031,'K')), NASAPolynomial(coeffs=[9.35661,0.011803,-6.15345e-06,1.22947e-09,-8.77631e-14,-66044.1,-16.6333], Tmin=(850.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-539.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]C(F)C(O)(F)C(F)(F)F(7349)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
10 C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1310.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0183691,0.0925734,-0.0001285,8.90877e-08,-2.43255e-11,-157490,28.8577], Tmin=(100,'K'), Tmax=(897.726,'K')), NASAPolynomial(coeffs=[15.7109,0.0226544,-1.1677e-05,2.33551e-09,-1.67366e-13,-160308,-45.1571], Tmin=(897.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1310.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(O2sj(Cs-CsF1sH))"""),
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
    label = '[O]C(F)C(F)=C(F)F(3844)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7 C u0 p0 c0 {2,S} {6,S} {8,D}
8 C u0 p0 c0 {3,S} {4,S} {7,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-688.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,323,467,575,827,1418,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.490926,'amu*angstrom^2'), symmetry=1, barrier=(11.2874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67943,0.0236114,0.000161177,-5.32832e-07,4.64216e-10,-82755.4,13.1547], Tmin=(10,'K'), Tmax=(426.943,'K')), NASAPolynomial(coeffs=[7.98521,0.0295583,-2.23418e-05,7.66806e-09,-9.73587e-13,-83544.9,-8.89409], Tmin=(426.943,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""[O]C(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)F(1205)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-401.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,451,553,637,1069,1180,1265,1301,3056],'cm^-1')),
        HinderedRotor(inertia=(0.412376,'amu*angstrom^2'), symmetry=1, barrier=(9.48133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0238,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91831,0.00720244,6.35725e-05,-2.08276e-07,2.06569e-10,-48296.9,9.02236], Tmin=(10,'K'), Tmax=(327.664,'K')), NASAPolynomial(coeffs=[3.62426,0.0171767,-1.13157e-05,3.55881e-09,-4.26934e-13,-48311.9,9.58992], Tmin=(327.664,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-401.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]C(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2CF2(61)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
"""),
    E0 = (-688.535,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(99.9936,'amu')),
        NonlinearRotor(inertia=([91.6969,155.94,247.638],'amu*angstrom^2'), symmetry=4),
        HarmonicOscillator(frequencies=([198.559,203.927,398.069,428.925,524.86,551.48,555.021,803.576,1208.21,1359.53,1361.83,1918.98],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86796,0.00896674,9.20781e-05,-2.50451e-07,1.90637e-10,-82806.8,9.05323], Tmin=(10,'K'), Tmax=(470.327,'K')), NASAPolynomial(coeffs=[5.38595,0.0191802,-1.42425e-05,4.78608e-09,-5.97116e-13,-83205.3,0.155967], Tmin=(470.327,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)C(F)(F)[CH]F(2999)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1137.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,180,2078.43],'cm^-1')),
        HinderedRotor(inertia=(0.483295,'amu*angstrom^2'), symmetry=1, barrier=(11.1119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01516,'amu*angstrom^2'), symmetry=1, barrier=(23.3404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33576,'amu*angstrom^2'), symmetry=1, barrier=(53.7037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29488,0.0648389,-7.73539e-05,4.73206e-08,-1.16919e-11,-136776,14.3559], Tmin=(10,'K'), Tmax=(946.697,'K')), NASAPolynomial(coeffs=[12.1948,0.0272346,-1.77712e-05,5.36206e-09,-6.11615e-13,-138461,-28.0935], Tmin=(946.697,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1137.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=CC(F)(F)C(O)(F)F(7350)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1192.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,223,363,546,575,694,1179,1410,2782.5,750,1395,475,1775,1000,223.109],'cm^-1')),
        HinderedRotor(inertia=(0.202518,'amu*angstrom^2'), symmetry=1, barrier=(7.15363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202519,'amu*angstrom^2'), symmetry=1, barrier=(7.15363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47291,'amu*angstrom^2'), symmetry=1, barrier=(16.7047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637862,0.079841,-0.000117252,8.99187e-08,-2.75455e-11,-143368,26.3863], Tmin=(100,'K'), Tmax=(798.615,'K')), NASAPolynomial(coeffs=[11.769,0.0240872,-1.25288e-05,2.49548e-09,-1.77513e-13,-145146,-24.8114], Tmin=(798.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1192.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-OdCsH)"""),
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
    label = 'OC(F)(F)[C](F)F(2998)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-907.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.391062,'amu*angstrom^2'), symmetry=1, barrier=(8.99127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390952,'amu*angstrom^2'), symmetry=1, barrier=(8.98875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56174,0.0389036,1.94218e-05,-2.22218e-07,2.51043e-10,-109199,12.6539], Tmin=(10,'K'), Tmax=(390.957,'K')), NASAPolynomial(coeffs=[7.58136,0.0231046,-1.7134e-05,5.81753e-09,-7.33941e-13,-109707,-5.43885], Tmin=(390.957,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-907.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OC(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(F)(F)C(O)(F)F(7351)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u0 p0 c0 {5,S} {7,D} {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1446.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,223,363,546,575,694,1179,1410,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.236469,'amu*angstrom^2'), symmetry=1, barrier=(5.43688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806119,'amu*angstrom^2'), symmetry=1, barrier=(18.5343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23649,'amu*angstrom^2'), symmetry=1, barrier=(5.43738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15287,0.0915216,-0.000146222,1.1957e-07,-3.86896e-11,-173900,29.0452], Tmin=(100,'K'), Tmax=(759.519,'K')), NASAPolynomial(coeffs=[12.9667,0.0240372,-1.29442e-05,2.5847e-09,-1.82964e-13,-175846,-29.2494], Tmin=(759.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1446.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(COCsFO)"""),
)

species(
    label = '[O]C(F)[C](F)C(O)(F)F(7352)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-886.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,391,562,707,872,1109,1210,1289,3137,212,367,445,1450,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.309599,'amu*angstrom^2'), symmetry=1, barrier=(7.11828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309842,'amu*angstrom^2'), symmetry=1, barrier=(7.12388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309376,'amu*angstrom^2'), symmetry=1, barrier=(7.11317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675482,0.0832861,-0.000117696,6.65699e-08,-2.45801e-12,-106569,28.9721], Tmin=(100,'K'), Tmax=(569.083,'K')), NASAPolynomial(coeffs=[11.7663,0.0243984,-1.27378e-05,2.49176e-09,-1.73453e-13,-108140,-20.9952], Tmin=(569.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-886.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(F)C(F)(F)[C](O)F(7353)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-875.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.289263,'amu*angstrom^2'), symmetry=1, barrier=(6.65073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289472,'amu*angstrom^2'), symmetry=1, barrier=(6.65554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289416,'amu*angstrom^2'), symmetry=1, barrier=(6.65423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00798713,0.101515,-0.000191916,1.77332e-07,-6.15128e-11,-105185,29.9001], Tmin=(100,'K'), Tmax=(868.953,'K')), NASAPolynomial(coeffs=[9.58997,0.0281454,-1.48791e-05,2.8523e-09,-1.93184e-13,-105751,-8.71549], Tmin=(868.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-875.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = '[O][CH]C(F)(F)C(O)(F)F(7354)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9  C u1 p0 c0 {6,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-909.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,3025,407.5,1350,352.5,180,1095.06,1095.85],'cm^-1')),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144401,'amu*angstrom^2'), symmetry=1, barrier=(3.32006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144745,'amu*angstrom^2'), symmetry=1, barrier=(3.32797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.55843,0.0867624,-0.000154124,1.45503e-07,-5.30314e-11,-109306,29.816], Tmin=(100,'K'), Tmax=(823.083,'K')), NASAPolynomial(coeffs=[6.99983,0.0338552,-1.83346e-05,3.62966e-09,-2.5382e-13,-109634,4.43945], Tmin=(823.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-909.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsOsHH) + radical(CCOJ) + radical(Csj(Cs-F1sF1sCs)(O2s-H)(H))"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(F)C(F)(F)[C](F)F(7355)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-878.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.421783,'amu*angstrom^2'), symmetry=1, barrier=(9.69762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.422796,'amu*angstrom^2'), symmetry=1, barrier=(9.72092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.113144,0.0973034,-0.000182433,1.67363e-07,-5.81069e-11,-105549,27.2412], Tmin=(100,'K'), Tmax=(851.609,'K')), NASAPolynomial(coeffs=[10.4127,0.025248,-1.38104e-05,2.71022e-09,-1.86782e-13,-106445,-15.7528], Tmin=(851.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-878.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)C(F)(F)C([O])(F)F(7356)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {9,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-1032.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,351,323,533,609,664,892,1120,1201,391,562,707,872,1109,1210,1289,3137,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.277005,'amu*angstrom^2'), symmetry=1, barrier=(6.36889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275696,'amu*angstrom^2'), symmetry=1, barrier=(6.3388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0498565,0.0959163,-0.000156633,1.29579e-07,-4.21077e-11,-124011,30.33], Tmin=(100,'K'), Tmax=(784.884,'K')), NASAPolynomial(coeffs=[13.6952,0.0231008,-1.21868e-05,2.39757e-09,-1.67663e-13,-126083,-32.1096], Tmin=(784.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1032.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFHO) + radical(O2sj(Cs-CsF1sF1s)) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O[C](F)F(2676)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-470.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.307072,'amu*angstrom^2'), symmetry=1, barrier=(7.06019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94331,0.003562,5.16413e-05,-1.18327e-07,7.89485e-11,-56559.5,10.5484], Tmin=(10,'K'), Tmax=(519.044,'K')), NASAPolynomial(coeffs=[4.2538,0.012858,-9.00307e-06,2.95251e-09,-3.63989e-13,-56749.2,7.73732], Tmin=(519.044,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-470.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)[C](F)F(429)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-445.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.813996,'amu*angstrom^2'), symmetry=1, barrier=(18.7154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14839,0.0397111,-4.32916e-05,2.26821e-08,-4.63769e-12,-53565.5,18.7797], Tmin=(100,'K'), Tmax=(1191.22,'K')), NASAPolynomial(coeffs=[11.353,0.00880293,-4.37161e-06,9.00449e-10,-6.63946e-14,-55758.5,-27.2377], Tmin=(1191.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = '[O][CH]F(388)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53592,0.00799825,-5.22271e-07,-4.56612e-09,2.08752e-12,-2348.33,9.31392], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97187,0.0038822,-1.62979e-06,3.36437e-10,-2.54188e-14,-3160.18,-3.94923], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O][C](F)C(F)(F)C(O)(F)F(7357)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1098.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,395,473,707,1436,180,180,3230.87],'cm^-1')),
        HinderedRotor(inertia=(0.364176,'amu*angstrom^2'), symmetry=1, barrier=(8.37312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3639,'amu*angstrom^2'), symmetry=1, barrier=(8.36678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363384,'amu*angstrom^2'), symmetry=1, barrier=(8.35492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.448554,0.1143,-0.000225478,2.10748e-07,-7.31192e-11,-131960,32.0854], Tmin=(100,'K'), Tmax=(878.514,'K')), NASAPolynomial(coeffs=[10.1238,0.0287559,-1.55478e-05,2.97207e-09,-1.9958e-13,-132374,-9.33499], Tmin=(878.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1098.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = '[O]C(F)C(F)(F)C(=O)F(7358)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {6,D} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-997.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,348,427,405,1245,1236,1280,391,562,707,872,1109,1210,1289,3137,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.334445,'amu*angstrom^2'), symmetry=1, barrier=(7.68954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337065,'amu*angstrom^2'), symmetry=1, barrier=(7.74979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88126,0.0742552,-0.000111051,8.70031e-08,-2.72536e-11,-119886,26.7342], Tmin=(100,'K'), Tmax=(781.196,'K')), NASAPolynomial(coeffs=[10.9671,0.0226105,-1.18826e-05,2.37099e-09,-1.68619e-13,-121462,-19.4333], Tmin=(781.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-997.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(COCsFO) + radical(O2sj(Cs-CsF1sH))"""),
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
    label = '[O]C(F)C(F)=C(O)F(7359)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-719.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,361,769,965,1078,1132,1246,3247,323,467,575,827,1418,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.637093,'amu*angstrom^2'), symmetry=1, barrier=(14.648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63795,'amu*angstrom^2'), symmetry=1, barrier=(14.6677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854226,0.0709516,-9.36217e-05,6.10018e-08,-1.55278e-11,-86417.8,23.8077], Tmin=(100,'K'), Tmax=(965.295,'K')), NASAPolynomial(coeffs=[14.2494,0.0154467,-7.37427e-06,1.43849e-09,-1.02131e-13,-89004,-40.3434], Tmin=(965.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sCdH))"""),
)

species(
    label = '[O]C=C(F)C(O)(F)F(7360)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-846.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,323,467,575,827,1418,3010,987.5,1337.5,450,1655,211.801,213.282],'cm^-1')),
        HinderedRotor(inertia=(0.301435,'amu*angstrom^2'), symmetry=1, barrier=(9.98816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.665904,'amu*angstrom^2'), symmetry=1, barrier=(21.507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922551,0.0698932,-9.00426e-05,5.73455e-08,-1.43421e-11,-101704,23.3073], Tmin=(100,'K'), Tmax=(978.993,'K')), NASAPolynomial(coeffs=[13.99,0.0165034,-8.24208e-06,1.64335e-09,-1.18172e-13,-104262,-39.4587], Tmin=(978.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-846.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'O=C(F)[C](F)C(O)(F)F(7361)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {9,S}
9  C u0 p0 c0 {4,S} {6,D} {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1087.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,234,347,1316,1464,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.497223,'amu*angstrom^2'), symmetry=1, barrier=(11.4321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58301,'amu*angstrom^2'), symmetry=1, barrier=(36.3966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497553,'amu*angstrom^2'), symmetry=1, barrier=(11.4397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817648,0.0749843,-0.000106464,7.6546e-08,-2.1818e-11,-130629,26.2821], Tmin=(100,'K'), Tmax=(858.636,'K')), NASAPolynomial(coeffs=[12.5923,0.0201305,-1.06349e-05,2.14018e-09,-1.53572e-13,-132651,-28.7289], Tmin=(858.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1087.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(COCsFO) + radical(CsCOCsF1s)"""),
)

species(
    label = 'O[C](F)C(F)(F)C(O)(F)F(7362)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1324.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.690507,0.118174,-0.000224607,2.06082e-07,-7.10011e-11,-159201,32.0507], Tmin=(100,'K'), Tmax=(868.042,'K')), NASAPolynomial(coeffs=[11.5029,0.0297303,-1.60366e-05,3.08944e-09,-2.09525e-13,-160102,-18.05], Tmin=(868.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1324.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFHO) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = '[O]C(F)(F)C(F)(F)C(O)F(7363)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1258.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,222,329,445,522,589,1214,1475,261,493,600,1152,1365,1422,3097,351,323,533,609,664,892,1120,1201,195.926,210.961],'cm^-1')),
        HinderedRotor(inertia=(0.237753,'amu*angstrom^2'), symmetry=1, barrier=(6.74208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250396,'amu*angstrom^2'), symmetry=1, barrier=(6.7191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.650803,'amu*angstrom^2'), symmetry=1, barrier=(19.7497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.140399,0.0978155,-0.000147888,1.13124e-07,-3.41545e-11,-151258,29.764], Tmin=(100,'K'), Tmax=(813.372,'K')), NASAPolynomial(coeffs=[14.7266,0.0247017,-1.3052e-05,2.60653e-09,-1.85376e-13,-153676,-38.8897], Tmin=(813.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1258.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'OOC(F)C(F)(F)[C](F)F(7364)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1034.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,215,315,519,588,595,1205,1248,261,493,600,1152,1365,1422,3097,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.83901,'amu*angstrom^2'), symmetry=1, barrier=(42.2825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619836,'amu*angstrom^2'), symmetry=1, barrier=(14.2512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618711,'amu*angstrom^2'), symmetry=1, barrier=(14.2254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83738,'amu*angstrom^2'), symmetry=1, barrier=(42.245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.624146,0.110931,-0.00018773,1.6132e-07,-5.43743e-11,-124289,30.876], Tmin=(100,'K'), Tmax=(792.709,'K')), NASAPolynomial(coeffs=[14.0201,0.0285348,-1.57282e-05,3.13727e-09,-2.20699e-13,-126343,-34.6873], Tmin=(792.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1034.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'OC(F)(F)[C](F)C(F)OF(7365)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u0 p2 c0 {5,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
10 C u1 p0 c0 {4,S} {8,S} {9,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-979.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,487,638,688,1119,1325,1387,3149,212,367,445,1450,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.72212,'amu*angstrom^2'), symmetry=1, barrier=(39.595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363086,'amu*angstrom^2'), symmetry=1, barrier=(8.34807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817007,'amu*angstrom^2'), symmetry=1, barrier=(18.7846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365796,'amu*angstrom^2'), symmetry=1, barrier=(8.41038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.234734,0.10269,-0.000160717,1.23403e-07,-3.53036e-11,-117608,31.4594], Tmin=(100,'K'), Tmax=(648.81,'K')), NASAPolynomial(coeffs=[13.9901,0.0273331,-1.50293e-05,3.02204e-09,-2.14683e-13,-119714,-33.0152], Tmin=(648.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-979.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsCsF1s)"""),
)

species(
    label = 'O[C](F)C(F)(F)C(F)OF(7366)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-967.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,215,315,519,588,595,1205,1248,261,493,600,1152,1365,1422,3097,395,473,707,1436,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.355253,'amu*angstrom^2'), symmetry=1, barrier=(8.16796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355286,'amu*angstrom^2'), symmetry=1, barrier=(8.16873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355181,'amu*angstrom^2'), symmetry=1, barrier=(8.16631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53336,'amu*angstrom^2'), symmetry=1, barrier=(35.255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.741337,0.117937,-0.000218455,1.98575e-07,-6.8591e-11,-116232,31.7998], Tmin=(100,'K'), Tmax=(849.599,'K')), NASAPolynomial(coeffs=[12.1588,0.030461,-1.67997e-05,3.2925e-09,-2.26787e-13,-117459,-22.6538], Tmin=(849.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-967.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCFHO) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = 'OC(F)(F)C(F)(F)[CH]OF(6729)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {5,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u1 p0 c0 {7,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-994.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,3025,407.5,1350,352.5,294.712,294.714,294.721],'cm^-1')),
        HinderedRotor(inertia=(0.167417,'amu*angstrom^2'), symmetry=1, barrier=(10.3191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406081,'amu*angstrom^2'), symmetry=1, barrier=(25.0299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167406,'amu*angstrom^2'), symmetry=1, barrier=(10.3191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167416,'amu*angstrom^2'), symmetry=1, barrier=(10.319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.666482,0.112157,-0.000181182,1.38547e-07,-3.8972e-11,-119480,30.8143], Tmin=(100,'K'), Tmax=(667.979,'K')), NASAPolynomial(coeffs=[16.7234,0.0224918,-1.2324e-05,2.44862e-09,-1.71834e-13,-122126,-48.4812], Tmin=(667.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-994.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsOsHH) + radical(CCsJO)"""),
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
    E0 = (-515.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-134.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-180.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-400.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (40.9794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-168.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-301.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-511.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-660.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-614.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-221.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-210.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-244.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-257.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-227.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-323.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-335.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-294.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-523.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-131.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-248.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-509.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-574.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-633.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-332.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-285.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-301.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-320.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    products = ['HF(38)', 'CF2O(49)', 'O=C[C](F)F(234)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.98542e+10,'s^-1'), n=0.978818, Ea=(183.685,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', '[O]C(F)C(O)(F)F(5575)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.22791e-07,'m^3/(mol*s)'), n=3.59921, Ea=(345.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node CC_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)F(364)', 'OC(F)(F)[C]F(2989)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(192.828,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    products = ['[O]C(F)C(O)(F)C(F)(F)F(7349)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(482)', '[O]C(F)C(F)=C(F)F(3844)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(O)F(1205)', 'CF2CF2(61)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.58e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/disub_Cd/disub;R_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', 'OC(F)(F)C(F)(F)[CH]F(2999)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'O=CC(F)(F)C(O)(F)F(7350)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(15.5688,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CHFO(47)', 'OC(F)(F)[C](F)F(2998)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.22525e-05,'m^3/(mol*s)'), n=2.9005, Ea=(48.8951,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_N-2R!H->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'O=C(F)C(F)(F)C(O)(F)F(7351)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(28.1096,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]C(F)[C](F)C(O)(F)F(7352)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(F)C(F)(F)[C](O)F(7353)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O][CH]C(F)(F)C(O)(F)F(7354)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OH(7)', '[O]C(F)C(F)(F)[C](F)F(7355)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[O]C(F)C(F)(F)C([O])(F)F(7356)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O[C](F)F(2676)', '[O]C(F)[C](F)F(429)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][CH]F(388)', 'OC(F)(F)[C](F)F(2998)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[O][C](F)C(F)(F)C(O)(F)F(7357)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[O]C(F)C(F)(F)C(=O)F(7358)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(162.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', '[O]C(F)C(F)=C(O)F(7359)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(4.29331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F2(78)', '[O]C=C(F)C(O)(F)F(7360)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(14.7727,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'O=C(F)[C](F)C(O)(F)F(7361)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(266.362,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    products = ['O[C](F)C(F)(F)C(O)(F)F(7362)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(F)C(F)(F)C(O)F(7363)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(214655,'s^-1'), n=1.70206, Ea=(32.737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;XH_out] for rate rule [R5H_CCC;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['OOC(F)C(F)(F)[C](F)F(7364)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.79e+10,'s^-1'), n=0, Ea=(110.039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3OOH_SS;C_rad_out_noH]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OC(F)(F)[C](F)C(F)OF(7365)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(101.176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O[C](F)C(F)(F)C(F)OF(7366)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(73.4342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OC(F)(F)C(F)(F)[CH]OF(6729)'],
    products = ['[O]C(F)C(F)(F)C(O)(F)F(6593)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(81.4373,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1946',
    isomers = [
        '[O]C(F)C(F)(F)C(O)(F)F(6593)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(49)', 'O=C[C](F)F(234)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1946',
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

