species(
    label = 'F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {9,S} {13,S}
11 C u0 p0 c0 {5,S} {6,S} {12,S} {14,S}
12 C u1 p0 c0 {9,S} {11,S} {15,S}
13 C u1 p0 c0 {7,S} {8,S} {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1382.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([181,249,286,344,486,552,511,665,553,637,1169,1241,1190,1306,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.206773,'amu*angstrom^2'), symmetry=1, barrier=(4.75412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206774,'amu*angstrom^2'), symmetry=1, barrier=(4.75414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206824,'amu*angstrom^2'), symmetry=1, barrier=(4.75528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567562,'amu*angstrom^2'), symmetry=1, barrier=(13.0494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97389,0.147378,-0.000265638,2.39567e-07,-8.28712e-11,-166055,41.4453], Tmin=(100,'K'), Tmax=(837.357,'K')), NASAPolynomial(coeffs=[14.3424,0.0404432,-2.21437e-05,4.35807e-09,-3.02127e-13,-167771,-28.3058], Tmin=(837.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1382.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2113.54,'J/mol'), sigma=(4.647,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86796,0.00896674,9.20781e-05,-2.50451e-07,1.90637e-10,-82806.8,9.05323], Tmin=(10,'K'), Tmax=(470.327,'K')), NASAPolynomial(coeffs=[5.38595,0.0191802,-1.42425e-05,4.78608e-09,-5.97116e-13,-83205.3,0.155967], Tmin=(470.327,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=CC(F)F(2834)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {3,S} {4,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-792.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.142102,'amu*angstrom^2'), symmetry=1, barrier=(3.2672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2858.21,'J/mol'), sigma=(4.57471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=446.45 K, Pc=67.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74532,0.0233589,0.000133412,-5.9431e-07,7.11464e-10,-95378,11.8568], Tmin=(10,'K'), Tmax=(302.26,'K')), NASAPolynomial(coeffs=[4.97982,0.0307765,-2.12832e-05,6.89182e-09,-8.42907e-13,-95561.1,5.58312], Tmin=(302.26,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-792.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)C(C(F)F)C(F)(F)[C](F)F(11892)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {12,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {10,S} {11,S} {12,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,S} {15,S}
12 C u1 p0 c0 {7,S} {8,S} {9,S}
13 C u1 p0 c0 {5,S} {6,S} {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1418.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,235,523,627,1123,1142,1372,1406,3097,146,234,414,562,504,606,1176,1296,1354,1460,197.964,198.096,199.445,199.776],'cm^-1')),
        HinderedRotor(inertia=(0.215102,'amu*angstrom^2'), symmetry=1, barrier=(6.20869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220961,'amu*angstrom^2'), symmetry=1, barrier=(6.21933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224252,'amu*angstrom^2'), symmetry=1, barrier=(6.20745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92869,'amu*angstrom^2'), symmetry=1, barrier=(26.3485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3144.24,'J/mol'), sigma=(5.32552,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=491.12 K, Pc=47.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02351,0.147202,-0.000261193,2.32035e-07,-7.94254e-11,-170447,40.9135], Tmin=(100,'K'), Tmax=(832.005,'K')), NASAPolynomial(coeffs=[15.4593,0.0385849,-2.10819e-05,4.15246e-09,-2.88319e-13,-172506,-35.1057], Tmin=(832.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1418.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = '[CH]C(F)F-2(3493)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-68.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,200.056,1350.78,1616.29],'cm^-1')),
        HinderedRotor(inertia=(0.00421206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96601,0.0195049,-1.18279e-05,2.09006e-09,3.19831e-13,-8190.14,12.8831], Tmin=(100,'K'), Tmax=(1211.04,'K')), NASAPolynomial(coeffs=[7.87163,0.00800485,-3.40879e-06,6.6204e-10,-4.73288e-14,-9723.19,-13.1468], Tmin=(1211.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)F(11908)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {9,S}
6 F u0 p3 c0 {9,S}
7 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
9 C u1 p0 c0 {5,S} {6,S} {7,S}
"""),
    E0 = (-917.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.339859,'amu*angstrom^2'), symmetry=1, barrier=(7.81403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647917,'amu*angstrom^2'), symmetry=1, barrier=(14.8969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242254,0.0921591,-0.000169904,1.51448e-07,-5.14291e-11,-110269,27.3579], Tmin=(100,'K'), Tmax=(842.414,'K')), NASAPolynomial(coeffs=[12.0021,0.019625,-1.10216e-05,2.18605e-09,-1.51606e-13,-111658,-23.8436], Tmin=(842.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-917.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[C]F(156)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C](F)C(F)(F)[CH]C(F)F(7253)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-971.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.213841,'amu*angstrom^2'), symmetry=1, barrier=(4.91663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21375,'amu*angstrom^2'), symmetry=1, barrier=(4.91453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213515,'amu*angstrom^2'), symmetry=1, barrier=(4.90913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.264982,0.106897,-0.000195154,1.80464e-07,-6.36546e-11,-116763,33.1133], Tmin=(100,'K'), Tmax=(841.38,'K')), NASAPolynomial(coeffs=[9.80404,0.0332253,-1.78111e-05,3.49776e-09,-2.42372e-13,-117544,-8.29738], Tmin=(841.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-971.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'FC(F)C1C(F)(F)C(F)(F)C1(F)F(11976)',
    structure = adjacencyList("""1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {10,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {9,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
13 C u0 p0 c0 {7,S} {8,S} {9,S} {15,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1687.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723675,0.108335,-0.000130281,7.90131e-08,-1.91034e-11,-202753,31.966], Tmin=(100,'K'), Tmax=(1003.72,'K')), NASAPolynomial(coeffs=[18.427,0.032016,-1.62278e-05,3.25996e-09,-2.3544e-13,-206597,-60.4963], Tmin=(1003.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1687.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsFFH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FC(F)=CC(F)(F)C(F)(F)C(F)F(11893)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {14,S}
12 C u0 p0 c0 {10,S} {13,D} {15,S}
13 C u0 p0 c0 {7,S} {8,S} {12,D}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1630.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48001,0.131346,-0.000215403,1.84357e-07,-6.24528e-11,-195932,36.7323], Tmin=(100,'K'), Tmax=(782.232,'K')), NASAPolynomial(coeffs=[14.9784,0.0384711,-2.05956e-05,4.08896e-09,-2.8781e-13,-198241,-36.9235], Tmin=(782.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1630.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(Cds-CdsCsH) + group(CdCFF)"""),
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
    label = 'F[C](F)C(F)(F)C(F)=CC(F)F(11987)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {11,S} {13,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {14,S}
12 C u1 p0 c0 {6,S} {7,S} {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1240.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.295588,'amu*angstrom^2'), symmetry=1, barrier=(6.79616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295641,'amu*angstrom^2'), symmetry=1, barrier=(6.79737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29599,'amu*angstrom^2'), symmetry=1, barrier=(6.80538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (195.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.987891,0.121551,-0.000207384,1.84332e-07,-6.40506e-11,-148995,36.1688], Tmin=(100,'K'), Tmax=(811.449,'K')), NASAPolynomial(coeffs=[12.403,0.0379729,-2.04086e-05,4.03576e-09,-2.82433e-13,-150590,-22.0726], Tmin=(811.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1240.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sF1sCd)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[C](F)F(4977)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
"""),
    E0 = (-493.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([146,234,414,562,504,606,1176,1296,1354,1460,180],'cm^-1')),
        HinderedRotor(inertia=(0.200092,'amu*angstrom^2'), symmetry=1, barrier=(4.60051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04127,0.0516237,-0.000104941,1.02886e-07,-3.74496e-11,-59293.5,18.0916], Tmin=(100,'K'), Tmax=(852.084,'K')), NASAPolynomial(coeffs=[5.53322,0.0162399,-9.21976e-06,1.83692e-09,-1.2751e-13,-59199.1,5.84951], Tmin=(852.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sF1s) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'FC(F)=C(F)C(F)(F)[CH]C(F)F(7321)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {8,S} {9,S} {14,S}
11 C u0 p0 c0 {5,S} {8,S} {12,D}
12 C u0 p0 c0 {6,S} {7,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1228.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,323,467,575,827,1418,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.243206,'amu*angstrom^2'), symmetry=1, barrier=(5.59177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243142,'amu*angstrom^2'), symmetry=1, barrier=(5.59032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243226,'amu*angstrom^2'), symmetry=1, barrier=(5.59223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (195.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11609,0.126096,-0.000222051,2.01512e-07,-7.09177e-11,-147642,36.3434], Tmin=(100,'K'), Tmax=(818.649,'K')), NASAPolynomial(coeffs=[11.9287,0.039549,-2.168e-05,4.30552e-09,-3.0141e-13,-149014,-19.3115], Tmin=(818.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1228.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Cs_S)"""),
)

species(
    label = 'F[C](F)[CH]C(F)F(4845)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u1 p0 c0 {5,S} {7,S} {9,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-554.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,190,488,555,1236,1407,180,2250.74],'cm^-1')),
        HinderedRotor(inertia=(0.397589,'amu*angstrom^2'), symmetry=1, barrier=(9.14135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399336,'amu*angstrom^2'), symmetry=1, barrier=(9.18152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4542,0.0670926,-0.000128678,1.28357e-07,-4.81479e-11,-66597.3,23.7219], Tmin=(100,'K'), Tmax=(836.732,'K')), NASAPolynomial(coeffs=[4.34924,0.0286339,-1.55993e-05,3.09731e-09,-2.16322e-13,-66219.9,15.4209], Tmin=(836.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FC=CC(F)(F)C(F)(F)[C](F)F(11988)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {5,S} {6,S} {9,S}
12 C u0 p0 c0 {7,S} {10,D} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-1226.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,215,315,519,588,595,1205,1248,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,194,682,905,1196,1383,3221,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.288402,'amu*angstrom^2'), symmetry=1, barrier=(6.63093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289626,'amu*angstrom^2'), symmetry=1, barrier=(6.65907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288879,'amu*angstrom^2'), symmetry=1, barrier=(6.6419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (195.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28379,0.12903,-0.000225377,2.00181e-07,-6.88041e-11,-147387,35.5917], Tmin=(100,'K'), Tmax=(827.896,'K')), NASAPolynomial(coeffs=[13.5274,0.0364685,-1.96229e-05,3.85734e-09,-2.68073e-13,-149120,-28.7198], Tmin=(827.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1226.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
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
    label = 'F[C](F)C(F)(F)C(F)(F)C=C(F)F(11989)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u0 p0 c0 {9,S} {13,D} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 C u0 p0 c0 {7,S} {8,S} {11,D}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1427.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,215,315,519,588,595,1205,1248,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.250033,'amu*angstrom^2'), symmetry=1, barrier=(5.74874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249809,'amu*angstrom^2'), symmetry=1, barrier=(5.74361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249902,'amu*angstrom^2'), symmetry=1, barrier=(5.74573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (213.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6408,0.138466,-0.000247626,2.21795e-07,-7.64466e-11,-171554,37.8794], Tmin=(100,'K'), Tmax=(832.157,'K')), NASAPolynomial(coeffs=[14.2684,0.0372207,-2.04726e-05,4.04141e-09,-2.80905e-13,-173344,-30.7965], Tmin=(832.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1427.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
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
    label = 'F[C](F)C(F)=C(F)[CH]C(F)F(7324)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u1 p0 c0 {7,S} {9,S} {13,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {11,S}
11 C u1 p0 c0 {5,S} {6,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-917.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,206,336,431,607,515,611,528,696,1312,1446,161,297,490,584,780,1358,180,824.13],'cm^-1')),
        HinderedRotor(inertia=(0.200404,'amu*angstrom^2'), symmetry=1, barrier=(4.60768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205439,'amu*angstrom^2'), symmetry=1, barrier=(4.72346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30015,'amu*angstrom^2'), symmetry=1, barrier=(52.885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263592,0.1033,-0.000154016,1.18356e-07,-3.49271e-11,-110189,31.7834], Tmin=(100,'K'), Tmax=(650.62,'K')), NASAPolynomial(coeffs=[12.6958,0.0339964,-1.81461e-05,3.63378e-09,-2.58594e-13,-112094,-26.8546], Tmin=(650.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-917.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Allyl_S) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
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
    label = 'F[C](F)CC(F)(F)C(F)(F)[C](F)F(11990)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
11 C u0 p0 c0 {9,S} {12,S} {14,S} {15,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 C u1 p0 c0 {7,S} {8,S} {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1389.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.14314,0.154178,-0.000288447,2.65656e-07,-9.27222e-11,-166901,40.3173], Tmin=(100,'K'), Tmax=(850.746,'K')), NASAPolynomial(coeffs=[13.387,0.0425158,-2.34337e-05,4.5928e-09,-3.16363e-13,-168145,-23.8774], Tmin=(850.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1389.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[CH]C(F)(F)C(F)(F)C(F)F(11991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {14,S}
12 C u1 p0 c0 {10,S} {13,S} {15,S}
13 C u1 p0 c0 {7,S} {8,S} {12,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1384.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80383,0.143563,-0.000257826,2.33973e-07,-8.15967e-11,-166317,41.6232], Tmin=(100,'K'), Tmax=(833.738,'K')), NASAPolynomial(coeffs=[13.3504,0.0419273,-2.29202e-05,4.51879e-09,-3.13981e-13,-167838,-22.7015], Tmin=(833.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1384.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)C(F)C(F)F(11992)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {11,S} {12,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {9,S} {15,S}
12 C u1 p0 c0 {6,S} {9,S} {10,S}
13 C u1 p0 c0 {7,S} {8,S} {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1378.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,190,488,555,1236,1407,180,180,180,180.821],'cm^-1')),
        HinderedRotor(inertia=(0.278126,'amu*angstrom^2'), symmetry=1, barrier=(6.39466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280742,'amu*angstrom^2'), symmetry=1, barrier=(6.45481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76777,'amu*angstrom^2'), symmetry=1, barrier=(40.6444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274649,'amu*angstrom^2'), symmetry=1, barrier=(6.31472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08581,0.150075,-0.000271209,2.43988e-07,-8.4036e-11,-165545,41.3421], Tmin=(100,'K'), Tmax=(841.362,'K')), NASAPolynomial(coeffs=[14.7893,0.0400292,-2.18583e-05,4.28905e-09,-2.9653e-13,-167329,-30.8838], Tmin=(841.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1378.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[C]([CH]C(F)F)C(F)(F)C(F)(F)F(11993)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
11 C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
12 C u1 p0 c0 {8,S} {9,S} {13,S}
13 C u1 p0 c0 {11,S} {12,S} {15,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1433.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53687,0.136059,-0.00023741,2.14908e-07,-7.58124e-11,-172273,41.4899], Tmin=(100,'K'), Tmax=(810.012,'K')), NASAPolynomial(coeffs=[12.67,0.0433283,-2.38846e-05,4.76299e-09,-3.34655e-13,-173834,-19.4851], Tmin=(810.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1433.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFFH) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)(F)C(F)C(F)F(11994)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {9,S} {15,S}
12 C u1 p0 c0 {6,S} {10,S} {13,S}
13 C u1 p0 c0 {7,S} {8,S} {12,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1383.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04265,0.150213,-0.000275005,2.50648e-07,-8.72172e-11,-166210,41.6727], Tmin=(100,'K'), Tmax=(842.168,'K')), NASAPolynomial(coeffs=[13.8236,0.0417851,-2.29838e-05,4.5208e-09,-3.1287e-13,-167709,-25.185], Tmin=(842.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1383.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](C(F)(F)F)C(F)(F)[CH]C(F)F(11995)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {12,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
11 C u0 p0 c0 {3,S} {6,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {9,S} {11,S}
13 C u1 p0 c0 {9,S} {10,S} {15,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1428.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56457,0.13635,-0.000237329,2.13139e-07,-7.44683e-11,-171565,41.2155], Tmin=(100,'K'), Tmax=(815.481,'K')), NASAPolynomial(coeffs=[13.185,0.0418774,-2.28572e-05,4.53501e-09,-3.17508e-13,-173235,-22.4235], Tmin=(815.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1428.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sF1s)(F1s)) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'F[CH]C(F)C(F)(F)C(F)(F)[C](F)F(11996)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
11 C u0 p0 c0 {5,S} {9,S} {12,S} {14,S}
12 C u1 p0 c0 {6,S} {11,S} {15,S}
13 C u1 p0 c0 {7,S} {8,S} {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1353.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,190,488,555,1236,1407,188.77,188.793,188.826,188.84],'cm^-1')),
        HinderedRotor(inertia=(0.319341,'amu*angstrom^2'), symmetry=1, barrier=(8.08057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319586,'amu*angstrom^2'), symmetry=1, barrier=(8.08067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31948,'amu*angstrom^2'), symmetry=1, barrier=(8.08063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54511,'amu*angstrom^2'), symmetry=1, barrier=(39.0911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37023,0.156943,-0.000286569,2.5669e-07,-8.76732e-11,-162621,40.2462], Tmin=(100,'K'), Tmax=(847.819,'K')), NASAPolynomial(coeffs=[16.0916,0.0383599,-2.10705e-05,4.12473e-09,-2.84054e-13,-164620,-39.1015], Tmin=(847.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1353.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][CH]C(F)(F)C(F)(F)C(F)(F)F(11997)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  F u0 p3 c0 {13,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
12 C u1 p0 c0 {10,S} {13,S} {14,S}
13 C u1 p0 c0 {8,S} {12,S} {15,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1411.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (214.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86152,0.14382,-0.000254881,2.28191e-07,-7.87733e-11,-169513,40.2503], Tmin=(100,'K'), Tmax=(831.163,'K')), NASAPolynomial(coeffs=[14.3318,0.0403125,-2.19238e-05,4.31702e-09,-2.99931e-13,-171321,-29.5641], Tmin=(831.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1411.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_noncyclic(CsF2-CsF2-CsF2) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFHH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(H)) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    E0 = (-526.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-369.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-130.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-82.3695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-518.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-501.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-256.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-418.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-248.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-381.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-247.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-358.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-189.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-57.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-319.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-398.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-493.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-342.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-324.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-299.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-366.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-330.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-342.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['CF2CF2(61)', 'FC(F)=CC(F)F(2834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C](F)C(C(F)F)C(F)(F)[C](F)F(11892)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(F)F-2(3493)', 'F[C](F)C(F)(F)[C](F)F(11908)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]F(156)', 'F[C](F)C(F)(F)[CH]C(F)F(7253)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['FC(F)C1C(F)(F)C(F)(F)C1(F)F(11976)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['FC(F)=CC(F)(F)C(F)(F)C(F)F(11893)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'F[C](F)C(F)(F)C(F)=CC(F)F(11987)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(55.0231,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C](F)[C](F)F(4977)', 'FC(F)=CC(F)F(2834)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(12.1732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'FC(F)=C(F)C(F)(F)[CH]C(F)F(7321)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.7953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2CF2(61)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.00392e-07,'m^3/(mol*s)'), n=3.22539, Ea=(6.002,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6344808260200242, var=1.3100047285217136, Tref=1000.0, N=132, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC=CC(F)(F)C(F)(F)[C](F)F(11988)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.1183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'F[C](F)C(F)(F)C(F)(F)C=C(F)F(11989)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.5227,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C](F)[C](F)F(4977)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(2.12259,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F2(78)', 'F[C](F)C(F)=C(F)[CH]C(F)F(7324)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(12.8028,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF2(43)', 'F[C](F)C(F)(F)[CH]C(F)F(7253)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C](F)CC(F)(F)C(F)(F)[C](F)F(11990)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C](F)[CH]C(F)(F)C(F)(F)C(F)F(11991)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(10500,'s^-1'), n=2.14, Ea=(33.3465,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_noH] for rate rule [R5HJ_3;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C](F)C(F)(F)[C](F)C(F)C(F)F(11992)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(179.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C]([CH]C(F)F)C(F)(F)C(F)(F)F(11993)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(201.97,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C](F)[C](F)C(F)(F)C(F)C(F)F(11994)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(227.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[C](C(F)(F)F)C(F)(F)[CH]C(F)F(11995)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(159.542,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(F)C(F)(F)C(F)(F)[C](F)F(11996)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)'],
    products = ['F[CH][CH]C(F)(F)C(F)(F)C(F)(F)F(11997)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(183.951,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3045',
    isomers = [
        'F[C](F)C(F)(F)C(F)(F)[CH]C(F)F(11891)',
    ],
    reactants = [
        ('CF2CF2(61)', 'FC(F)=CC(F)F(2834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3045',
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

