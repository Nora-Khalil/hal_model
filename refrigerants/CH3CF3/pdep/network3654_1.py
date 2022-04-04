species(
    label = 'F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {14,S}
12 C u1 p0 c0 {7,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1229.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,334,575,1197,1424,3202,180,180,180,2137.8],'cm^-1')),
        HinderedRotor(inertia=(0.223515,'amu*angstrom^2'), symmetry=1, barrier=(5.13905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223501,'amu*angstrom^2'), symmetry=1, barrier=(5.13873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223525,'amu*angstrom^2'), symmetry=1, barrier=(5.13927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51429,'amu*angstrom^2'), symmetry=1, barrier=(34.8164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28168,0.128991,-0.000220207,1.9697e-07,-6.90078e-11,-147674,36.7401], Tmin=(100,'K'), Tmax=(806.347,'K')), NASAPolynomial(coeffs=[12.5179,0.0416346,-2.25431e-05,4.47681e-09,-3.1426e-13,-149285,-23.0545], Tmin=(806.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1229.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFF) + group(CsCsFHH) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'CHFCF2(55)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-511.455,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(82.003,'amu')),
        NonlinearRotor(inertia=([47.283,130.848,178.131],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([226.38,309.801,488.171,585.738,628.102,776.692,950.625,1195.27,1295.73,1386.74,1841.68,3239.84],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93201,0.00413887,7.03233e-05,-1.49827e-07,9.42397e-11,-61509.2,9.50863], Tmin=(10,'K'), Tmax=(540.182,'K')), NASAPolynomial(coeffs=[3.82032,0.0197002,-1.38027e-05,4.49179e-09,-5.49698e-13,-61712.1,7.9889], Tmin=(540.182,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=CC(F)(F)F(1027)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u0 p0 c0 {4,S} {6,D} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-831.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,224.569],'cm^-1')),
        HinderedRotor(inertia=(0.0840844,'amu*angstrom^2'), symmetry=1, barrier=(1.93327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3084.67,'J/mol'), sigma=(4.9,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77934,0.0173102,9.62417e-05,-2.61573e-07,1.94701e-10,-100018,11.5092], Tmin=(10,'K'), Tmax=(474.505,'K')), NASAPolynomial(coeffs=[4.81364,0.0316651,-2.20773e-05,7.14072e-09,-8.67875e-13,-100376,4.55317], Tmin=(474.505,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-831.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FCDCC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(C(F)(F)F)C(F)(F)[CH]F(7029)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {8,S} {14,S}
12 C u1 p0 c0 {7,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1232.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,193,295,551,588,656,1146,1192,1350,262,406,528,622,1148,1246,1368,1480,3164,3240,202.519,202.52,202.52,1839.43],'cm^-1')),
        HinderedRotor(inertia=(0.246914,'amu*angstrom^2'), symmetry=1, barrier=(7.18639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00411022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246915,'amu*angstrom^2'), symmetry=1, barrier=(7.18639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40205,'amu*angstrom^2'), symmetry=1, barrier=(40.8065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3173.38,'J/mol'), sigma=(5.84662,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.67 K, Pc=36.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4343,0.132147,-0.000225733,1.99309e-07,-6.85888e-11,-148004,37.2288], Tmin=(100,'K'), Tmax=(818.736,'K')), NASAPolynomial(coeffs=[13.6081,0.039544,-2.10611e-05,4.14151e-09,-2.88606e-13,-149826,-28.4216], Tmin=(818.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1232.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFHH) + group(CsCsFHH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(F)C(F)[CH]C(F)(F)F(12116)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {8,S} {12,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u1 p0 c0 {6,S} {7,S} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1221.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([207,311,436,622,447,691,1074,1182,1276,1366,1311,1469,3045,3235,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,190,488,555,1236,1407,250.487,250.519,250.556,1690.06],'cm^-1')),
        HinderedRotor(inertia=(0.146723,'amu*angstrom^2'), symmetry=1, barrier=(6.53625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146736,'amu*angstrom^2'), symmetry=1, barrier=(6.53638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146801,'amu*angstrom^2'), symmetry=1, barrier=(6.53623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786952,'amu*angstrom^2'), symmetry=1, barrier=(35.028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3192.38,'J/mol'), sigma=(5.58068,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=498.64 K, Pc=41.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04199,0.121879,-0.000198908,1.74364e-07,-6.09451e-11,-146795,36.9815], Tmin=(100,'K'), Tmax=(774.54,'K')), NASAPolynomial(coeffs=[12.4406,0.0414858,-2.23695e-05,4.46919e-09,-3.16187e-13,-148560,-22.5351], Tmin=(774.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1221.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFF) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[CH]C(F)(F)F(462)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u2 p0 c0 {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-319.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([193,295,551,588,656,1146,1192,1350,180,747.919,1907.62],'cm^-1')),
        HinderedRotor(inertia=(0.00850025,'amu*angstrom^2'), symmetry=1, barrier=(3.32661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44977,0.0297269,-2.96856e-05,1.40699e-08,-2.55624e-12,-38337.5,14.102], Tmin=(100,'K'), Tmax=(1355.02,'K')), NASAPolynomial(coeffs=[10.9803,0.0045449,-1.80947e-06,3.54966e-10,-2.58626e-14,-40649.3,-29.6449], Tmin=(1355.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]F(1129)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 C u1 p0 c0 {4,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-512.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,262,406,528,622,1148,1246,1368,1480,3164,3240,295.621,295.798],'cm^-1')),
        HinderedRotor(inertia=(0.139232,'amu*angstrom^2'), symmetry=1, barrier=(8.73602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139328,'amu*angstrom^2'), symmetry=1, barrier=(8.73931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958965,0.0727702,-0.00012328,1.05184e-07,-3.4846e-11,-61489.6,22.3745], Tmin=(100,'K'), Tmax=(820.009,'K')), NASAPolynomial(coeffs=[10.8923,0.0172825,-8.91369e-06,1.7449e-09,-1.2131e-13,-62882.2,-22.1351], Tmin=(820.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = '[CH]F(804)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'F[C](F)C(F)[CH]C(F)(F)F(1140)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1004.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,738.616],'cm^-1')),
        HinderedRotor(inertia=(0.203364,'amu*angstrom^2'), symmetry=1, barrier=(4.67574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203368,'amu*angstrom^2'), symmetry=1, barrier=(4.67583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203349,'amu*angstrom^2'), symmetry=1, barrier=(4.6754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0384801,0.0960041,-0.000159619,1.40805e-07,-4.90751e-11,-120719,30.6515], Tmin=(100,'K'), Tmax=(788.794,'K')), NASAPolynomial(coeffs=[10.7976,0.0311133,-1.65743e-05,3.30396e-09,-2.33216e-13,-122094,-16.6653], Tmin=(788.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1004.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC1C(C(F)(F)F)C(F)C1(F)F(12267)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {8,S} {11,S} {15,S}
11 C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
12 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1493.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206135,0.0889015,-8.61432e-05,4.14114e-08,-7.92183e-12,-179493,30.4259], Tmin=(100,'K'), Tmax=(1256.16,'K')), NASAPolynomial(coeffs=[19.0739,0.0275069,-1.283e-05,2.50216e-09,-1.78024e-13,-184337,-66.9859], Tmin=(1256.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1493.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsFFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FCC(F)(F)C(F)=CC(F)(F)F(12299)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {12,S}
11 C u0 p0 c0 {7,S} {8,S} {12,D}
12 C u0 p0 c0 {10,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1477.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.49731,0.106319,-0.000139202,9.5527e-08,-2.64278e-11,-177537,32.5961], Tmin=(100,'K'), Tmax=(877.845,'K')), NASAPolynomial(coeffs=[15.041,0.0355167,-1.82197e-05,3.64869e-09,-2.61939e-13,-180265,-40.343], Tmin=(877.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1477.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH)"""),
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
    label = 'F[CH]C(F)(F)C=CC(F)(F)F(12300)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u0 p0 c0 {8,S} {9,D} {13,S}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1090.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,219,296,586,564,718,793,1177,1228,2995,3025,975,1000,1300,1375,400,500,1630,1680,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.142934,'amu*angstrom^2'), symmetry=1, barrier=(3.28633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144763,'amu*angstrom^2'), symmetry=1, barrier=(3.32838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44806,'amu*angstrom^2'), symmetry=1, barrier=(10.3018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00994805,0.0954185,-0.000126124,8.93027e-08,-2.57507e-11,-131034,32.8071], Tmin=(100,'K'), Tmax=(839.976,'K')), NASAPolynomial(coeffs=[12.6977,0.034999,-1.82289e-05,3.66969e-09,-2.63988e-13,-133166,-26.1915], Tmin=(839.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1090.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CsCdFFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
)

species(
    label = 'F[CH][C](F)F(588)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1767.79],'cm^-1')),
        HinderedRotor(inertia=(0.511633,'amu*angstrom^2'), symmetry=1, barrier=(11.7635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54625,0.0360263,-6.18605e-05,5.68976e-08,-2.03344e-11,-34259.7,17.1201], Tmin=(100,'K'), Tmax=(816.909,'K')), NASAPolynomial(coeffs=[5.79758,0.0130515,-6.72073e-06,1.3276e-09,-9.30772e-14,-34555.5,3.53259], Tmin=(816.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(F)(F)C(F)=CC(F)(F)F(12301)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {13,S}
12 C u1 p0 c0 {7,S} {8,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-1276.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,219,296,586,564,718,793,1177,1228,323,467,575,827,1418,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.172897,'amu*angstrom^2'), symmetry=1, barrier=(3.97524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172686,'amu*angstrom^2'), symmetry=1, barrier=(3.9704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172525,'amu*angstrom^2'), symmetry=1, barrier=(3.96668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (195.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.438593,0.10506,-0.000144507,1.03193e-07,-2.95189e-11,-153352,35.4153], Tmin=(100,'K'), Tmax=(852.23,'K')), NASAPolynomial(coeffs=[15.0046,0.0325768,-1.693e-05,3.39419e-09,-2.43319e-13,-155984,-36.62], Tmin=(852.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1276.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
)

species(
    label = 'FC=C(F)C(F)[CH]C(F)(F)F(12302)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1068.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,323,467,575,827,1418,194,682,905,1196,1383,3221,180,1981.3],'cm^-1')),
        HinderedRotor(inertia=(0.26209,'amu*angstrom^2'), symmetry=1, barrier=(6.02597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261988,'amu*angstrom^2'), symmetry=1, barrier=(6.02361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26055,'amu*angstrom^2'), symmetry=1, barrier=(5.99056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.372824,0.10608,-0.000171459,1.51667e-07,-5.35596e-11,-128331,33.8043], Tmin=(100,'K'), Tmax=(781.806,'K')), NASAPolynomial(coeffs=[10.4285,0.0391773,-2.07678e-05,4.12626e-09,-2.91124e-13,-129664,-13.3722], Tmin=(781.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1068.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFF) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_S)"""),
)

species(
    label = 'F[CH][CH]C(F)(F)F(9334)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u1 p0 c0 {5,S} {7,S} {8,S}
7 C u1 p0 c0 {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-581.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,334,575,1197,1424,3202,208.405,1612.5],'cm^-1')),
        HinderedRotor(inertia=(0.236761,'amu*angstrom^2'), symmetry=1, barrier=(6.44327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00345971,'amu*angstrom^2'), symmetry=1, barrier=(6.42129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9911,0.0487013,-6.31668e-05,4.72978e-08,-1.48179e-11,-69825.2,23.0747], Tmin=(100,'K'), Tmax=(768.081,'K')), NASAPolynomial(coeffs=[7.09051,0.0221438,-1.13006e-05,2.27841e-09,-1.64263e-13,-70608.5,-0.181341], Tmin=(768.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)(F)C(F)C=C(F)F(9349)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1027.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,215,315,519,588,595,1205,1248,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.122936,'amu*angstrom^2'), symmetry=1, barrier=(2.82654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125574,'amu*angstrom^2'), symmetry=1, barrier=(2.8872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501522,'amu*angstrom^2'), symmetry=1, barrier=(11.531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744547,0.115809,-0.00019571,1.74597e-07,-6.1043e-11,-123397,33.5575], Tmin=(100,'K'), Tmax=(810.309,'K')), NASAPolynomial(coeffs=[11.4139,0.0384671,-2.04726e-05,4.04103e-09,-2.82783e-13,-124799,-19.0332], Tmin=(810.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1027.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    collisionModel = TransportData(shapeIndex=1, epsilon=(1045.13,'J/mol'), sigma=(3.301,'angstroms'), dipoleMoment=(0,'De'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]C(F)=C[CH]C(F)(F)F(12216)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,D} {12,S}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 C u1 p0 c0 {5,S} {9,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-798.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,271,519,563,612,1379,234,589,736,816,1240,3237,759.156,2511.88],'cm^-1')),
        HinderedRotor(inertia=(0.316752,'amu*angstrom^2'), symmetry=1, barrier=(7.28275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88802,'amu*angstrom^2'), symmetry=1, barrier=(43.4092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88781,'amu*angstrom^2'), symmetry=1, barrier=(43.4044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516089,0.0805238,-9.07743e-05,5.32749e-08,-1.26853e-11,-95868.2,27.126], Tmin=(100,'K'), Tmax=(1009.83,'K')), NASAPolynomial(coeffs=[13.4925,0.0291223,-1.44214e-05,2.86764e-09,-2.05899e-13,-98488.9,-35.6046], Tmin=(1009.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-798.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
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
    label = 'F[CH]C(F)=C(F)[CH]C(F)(F)F(12303)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {7,S} {9,S} {12,S}
9  C u0 p0 c0 {4,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {11,S}
11 C u1 p0 c0 {6,S} {10,S} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-975.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,206,336,431,607,515,611,528,696,1312,1446,234,589,736,816,1240,3237,180,2085.7],'cm^-1')),
        HinderedRotor(inertia=(0.0140777,'amu*angstrom^2'), symmetry=1, barrier=(43.4474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135081,'amu*angstrom^2'), symmetry=1, barrier=(9.93639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88932,'amu*angstrom^2'), symmetry=1, barrier=(43.4393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0643374,0.0922046,-0.000116772,7.62192e-08,-1.9977e-11,-117191,29.6382], Tmin=(100,'K'), Tmax=(925.983,'K')), NASAPolynomial(coeffs=[14.5901,0.0294564,-1.51249e-05,3.03648e-09,-2.18636e-13,-119882,-39.3228], Tmin=(925.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-975.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Allyl_S) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
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
    label = 'F[CH]C(F)(F)[C](F)CC(F)(F)F(12304)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {10,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {8,S} {9,S}
12 C u1 p0 c0 {7,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1241.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48441,0.136786,-0.000246894,2.26372e-07,-7.94157e-11,-149106,38.4383], Tmin=(100,'K'), Tmax=(841.773,'K')), NASAPolynomial(coeffs=[11.7178,0.0425498,-2.28368e-05,4.46768e-09,-3.08853e-13,-150212,-16.3517], Tmin=(841.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1241.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + group(CsCsFHH) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FCC(F)(F)[C](F)[CH]C(F)(F)F(12305)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {12,S}
11 C u1 p0 c0 {7,S} {8,S} {12,S}
12 C u1 p0 c0 {10,S} {11,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1238.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17463,0.126318,-0.000213785,1.91474e-07,-6.74689e-11,-148732,37.2308], Tmin=(100,'K'), Tmax=(798.705,'K')), NASAPolynomial(coeffs=[12.078,0.0423273,-2.29561e-05,4.57202e-09,-3.21866e-13,-150287,-20.2085], Tmin=(798.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1238.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CsCsFFF) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]C(F)C(F)(F)F(12293)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
11 C u1 p0 c0 {8,S} {9,S} {14,S}
12 C u1 p0 c0 {7,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1229.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4194,0.134055,-0.000237565,2.16081e-07,-7.58312e-11,-147663,38.237], Tmin=(100,'K'), Tmax=(829.09,'K')), NASAPolynomial(coeffs=[12.1246,0.0420076,-2.27207e-05,4.47836e-09,-3.11727e-13,-148991,-19.0314], Tmin=(829.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1229.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCsFHH) + radical(Cs_S) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FC(F)C(F)(F)[CH][CH]C(F)(F)F(12306)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {7,S} {12,S}
11 C u1 p0 c0 {8,S} {12,S} {14,S}
12 C u1 p0 c0 {10,S} {11,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1251.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.837934,0.117414,-0.000191208,1.6951e-07,-6.00165e-11,-150330,39.7931], Tmin=(100,'K'), Tmax=(775.884,'K')), NASAPolynomial(coeffs=[11.3007,0.0425985,-2.29118e-05,4.57764e-09,-3.23979e-13,-151845,-13.3152], Tmin=(775.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1251.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFF) + radical(Cs_S) + radical(Csj(Cs-CsHH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = 'F[CH][C](F)C(F)C(F)C(F)(F)F(12307)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {9,S} {12,S}
12 C u1 p0 c0 {7,S} {11,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1196.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,259,529,569,1128,1321,1390,3140,193,295,551,588,656,1146,1192,1350,212,367,445,1450,334,575,1197,1424,3202,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.111237,'amu*angstrom^2'), symmetry=1, barrier=(2.55755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111122,'amu*angstrom^2'), symmetry=1, barrier=(2.55492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110991,'amu*angstrom^2'), symmetry=1, barrier=(2.5519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33424,'amu*angstrom^2'), symmetry=1, barrier=(30.6767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52735,0.136143,-0.000240165,2.16307e-07,-7.51084e-11,-143698,38.411], Tmin=(100,'K'), Tmax=(833.78,'K')), NASAPolynomial(coeffs=[12.8992,0.0408717,-2.1886e-05,4.29014e-09,-2.97471e-13,-145198,-23.1342], Tmin=(833.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1196.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sCsH)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[C](C(F)F)C(F)[CH]C(F)(F)F(12286)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {12,S}
10 C u0 p0 c0 {5,S} {6,S} {11,S} {14,S}
11 C u1 p0 c0 {7,S} {8,S} {10,S}
12 C u1 p0 c0 {8,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1224.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.883416,0.121968,-0.000212218,1.98345e-07,-7.22109e-11,-147160,37.8401], Tmin=(100,'K'), Tmax=(811.959,'K')), NASAPolynomial(coeffs=[8.85886,0.0475099,-2.57771e-05,5.12515e-09,-3.59979e-13,-147870,-1.75951], Tmin=(811.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1224.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCsFFF) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = 'F[CH]C(F)(F)C(F)C(F)[C](F)F(12308)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {8,S} {11,S} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {10,S}
12 C u1 p0 c0 {7,S} {9,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1161.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,334,575,1197,1424,3202,204.108,204.108,204.109,204.109],'cm^-1')),
        HinderedRotor(inertia=(0.190294,'amu*angstrom^2'), symmetry=1, barrier=(5.62568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190296,'amu*angstrom^2'), symmetry=1, barrier=(5.62569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190293,'amu*angstrom^2'), symmetry=1, barrier=(5.62567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41976,'amu*angstrom^2'), symmetry=1, barrier=(41.9724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66948,0.139644,-0.000247883,2.2302e-07,-7.71985e-11,-139460,38.5198], Tmin=(100,'K'), Tmax=(836.375,'K')), NASAPolynomial(coeffs=[13.4216,0.04052,-2.17762e-05,4.26691e-09,-2.95506e-13,-141041,-25.9547], Tmin=(836.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1161.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'F[C](F)[CH]C(F)C(F)(F)C(F)F(12309)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {14,S}
11 C u1 p0 c0 {9,S} {12,S} {15,S}
12 C u1 p0 c0 {6,S} {7,S} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1187.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,3025,407.5,1350,352.5,190,488,555,1236,1407,223.137,223.137,223.14,1302.18],'cm^-1')),
        HinderedRotor(inertia=(0.127867,'amu*angstrom^2'), symmetry=1, barrier=(4.51783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127865,'amu*angstrom^2'), symmetry=1, barrier=(4.51784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127868,'amu*angstrom^2'), symmetry=1, barrier=(4.51785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15288,'amu*angstrom^2'), symmetry=1, barrier=(40.7342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17751,0.126546,-0.000214233,1.92327e-07,-6.79073e-11,-142666,39.4027], Tmin=(100,'K'), Tmax=(800.056,'K')), NASAPolynomial(coeffs=[11.8738,0.0429981,-2.32896e-05,4.63513e-09,-3.26157e-13,-144169,-16.9912], Tmin=(800.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1187.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    E0 = (-480.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-321.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-313.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-83.0042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-41.5349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-472.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-417.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-213.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-359.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-309.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-197.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-338.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-166.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-118.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-41.7707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-267.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-117.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-353.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-325.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-298.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-264.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-236.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-297.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-263.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-261.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['CHFCF2(55)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['F[CH]C(C(F)(F)F)C(F)(F)[CH]F(7029)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C](F)C(F)C(F)[CH]C(F)(F)F(12116)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(F)(F)F(462)', 'F[CH]C(F)(F)[CH]F(1129)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(804)', 'F[C](F)C(F)[CH]C(F)(F)F(1140)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['FC1C(C(F)(F)F)C(F)C1(F)F(12267)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['FCC(F)(F)C(F)=CC(F)(F)F(12299)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'F[CH]C(F)(F)C=CC(F)(F)F(12300)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(56.0705,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][C](F)F(588)', 'FC=CC(F)(F)F(1027)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(8.97052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'F[CH]C(F)(F)C(F)=CC(F)(F)F(12301)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.72227,'m^3/(mol*s)'), n=2.06948, Ea=(6.54133,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.15211669081533294, var=0.5516637502444228, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC=C(F)C(F)[CH]C(F)(F)F(12302)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(49.796,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHFCF2(55)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(5.33236,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'F[CH]C(F)(F)C(F)C=C(F)F(9349)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[CH][C](F)F(588)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -6.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['F2(78)', 'F[CH]C(F)=C[CH]C(F)(F)F(12216)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(16.7808,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'F[CH]C(F)=C(F)[CH]C(F)(F)F(12303)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(241.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHF(40)', 'F[C](F)C(F)[CH]C(F)(F)F(1140)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['F[CH]C(F)(F)[C](F)CC(F)(F)F(12304)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['FCC(F)(F)[C](F)[CH]C(F)(F)F(12305)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(91.367,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['F[CH]C(F)(F)[CH]C(F)C(F)(F)F(12293)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.556,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['FC(F)C(F)(F)[CH][CH]C(F)(F)F(12306)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(216.04,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH][C](F)C(F)C(F)C(F)(F)F(12307)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(211.038,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    products = ['F[C](C(F)F)C(F)[CH]C(F)(F)F(12286)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(183.877,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(F)(F)C(F)C(F)[C](F)F(12308)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(149.382,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C](F)[CH]C(F)C(F)(F)C(F)F(12309)'],
    products = ['F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(178.173,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3654',
    isomers = [
        'F[CH]C(F)(F)C(F)[CH]C(F)(F)F(12115)',
    ],
    reactants = [
        ('CHFCF2(55)', 'FC=CC(F)(F)F(1027)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3654',
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

