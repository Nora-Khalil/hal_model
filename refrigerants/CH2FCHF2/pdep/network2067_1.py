species(
    label = 'FOC(F)(F)C(F)(F)CC(F)F(7684)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {8,S}
8  O u0 p2 c0 {7,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
12 C u0 p0 c0 {5,S} {6,S} {10,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1452.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,269.377,269.377,269.377,269.379],'cm^-1')),
        HinderedRotor(inertia=(0.205354,'amu*angstrom^2'), symmetry=1, barrier=(10.5743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205354,'amu*angstrom^2'), symmetry=1, barrier=(10.5743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205359,'amu*angstrom^2'), symmetry=1, barrier=(10.5744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628644,'amu*angstrom^2'), symmetry=1, barrier=(32.3702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (200.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31086,0.12819,-0.000198276,1.54478e-07,-4.61079e-11,-174545,34.192], Tmin=(100,'K'), Tmax=(663.199,'K')), NASAPolynomial(coeffs=[15.9998,0.0365171,-1.97355e-05,3.95619e-09,-2.81066e-13,-177122,-44.3242], Tmin=(663.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1452.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH)"""),
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
    label = 'FOC(F)(F)CC(F)F(2803)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1014.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,253.948,254.028,254.512],'cm^-1')),
        HinderedRotor(inertia=(0.135369,'amu*angstrom^2'), symmetry=1, barrier=(6.21199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13581,'amu*angstrom^2'), symmetry=1, barrier=(6.22595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445322,'amu*angstrom^2'), symmetry=1, barrier=(20.4881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.047,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.45,'J/mol'), sigma=(5.12178,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.04 K, Pc=49.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28621,0.0637681,-6.80208e-05,3.70983e-08,-8.16866e-12,-121994,14.6088], Tmin=(10,'K'), Tmax=(1055.22,'K')), NASAPolynomial(coeffs=[12.754,0.0278787,-1.70035e-05,4.86637e-09,-5.32308e-13,-123992,-31.5764], Tmin=(1055.22,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1014.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""FOC(F)(F)CC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC(F)(F)F(2825)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
"""),
    E0 = (-757.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,431,527,613,668,835,1199,1245,1322,180],'cm^-1')),
        HinderedRotor(inertia=(0.86978,'amu*angstrom^2'), symmetry=1, barrier=(19.9979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88199,0.0071794,8.93998e-05,-2.06387e-07,1.33773e-10,-91085.4,10.4612], Tmin=(10,'K'), Tmax=(550.265,'K')), NASAPolynomial(coeffs=[5.53551,0.0205748,-1.6396e-05,5.72482e-09,-7.29448e-13,-91652.1,-0.0246061], Tmin=(550.265,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-757.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""FOC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CC(F)F(718)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
6 C u0 p1 c0 {3,S} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-353.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.335793,'amu*angstrom^2'), symmetry=1, barrier=(7.72053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710453,'amu*angstrom^2'), symmetry=1, barrier=(16.3347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1379,0.0437333,-4.65407e-05,2.72468e-08,-6.65411e-12,-42416.9,18.0261], Tmin=(100,'K'), Tmax=(972.998,'K')), NASAPolynomial(coeffs=[7.98243,0.0197063,-9.4999e-06,1.86757e-09,-1.33199e-13,-43554.2,-10.0104], Tmin=(972.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsFFH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'F[C]C(F)(F)OF(2807)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u0 p1 c0 {3,S} {6,S}
"""),
    E0 = (-404.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.52222,'amu*angstrom^2'), symmetry=1, barrier=(34.9989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52258,'amu*angstrom^2'), symmetry=1, barrier=(35.0071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52555,0.0589818,-9.49276e-05,7.55255e-08,-2.36062e-11,-48604.5,19.5039], Tmin=(100,'K'), Tmax=(787.109,'K')), NASAPolynomial(coeffs=[10.5434,0.0131539,-7.59296e-06,1.55453e-09,-1.11631e-13,-50024.1,-21.8434], Tmin=(787.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'CH2FCHF2(1)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-683.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([528,1116,1182,1331,1402,1494,3075,3110,235,523,627,1123,1142,1372,1406,3097,180],'cm^-1')),
        HinderedRotor(inertia=(0.105901,'amu*angstrom^2'), symmetry=1, barrier=(14.2006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0404,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85089,0.0125534,3.70918e-05,-7.49617e-08,4.07443e-11,-82160.9,10.5553], Tmin=(10,'K'), Tmax=(616.246,'K')), NASAPolynomial(coeffs=[3.0356,0.0261873,-1.63995e-05,4.87265e-09,-5.54431e-13,-82218.8,12.8087], Tmin=(616.246,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-683.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""FCC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)C(F)(F)C(F)F(7711)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {8,S}
8  O u0 p2 c0 {7,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-1372.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28449,'amu*angstrom^2'), symmetry=1, barrier=(29.5329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369438,'amu*angstrom^2'), symmetry=1, barrier=(8.49412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369145,'amu*angstrom^2'), symmetry=1, barrier=(8.48736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.028,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.13867,0.0750208,-4.82475e-05,-7.04278e-08,8.47378e-11,-165038,15.782], Tmin=(10,'K'), Tmax=(537.347,'K')), NASAPolynomial(coeffs=[12.2318,0.0375834,-2.81875e-05,9.45405e-09,-1.17103e-12,-166452,-26.5029], Tmin=(537.347,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1372.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""FOC(F)(F)C(F)(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]C(F)(F)CC(F)F(7712)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
9  C u0 p1 c0 {5,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-803.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,347,453,1141,1468,235,523,627,1123,1142,1372,1406,3097,617,898,1187,225.185,225.185,225.187,225.187,225.189],'cm^-1')),
        HinderedRotor(inertia=(0.180947,'amu*angstrom^2'), symmetry=1, barrier=(6.51131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0033244,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540119,'amu*angstrom^2'), symmetry=1, barrier=(19.4359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465024,0.0847435,-0.000123724,9.74835e-08,-3.11139e-11,-96510.5,26.7172], Tmin=(100,'K'), Tmax=(763.538,'K')), NASAPolynomial(coeffs=[10.963,0.029746,-1.56777e-05,3.14404e-09,-2.24593e-13,-98113.6,-21.0972], Tmin=(763.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-803.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'CC(F)(F)C(F)(F)OF(2804)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
9  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1014.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,223,363,546,575,694,1179,1410,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.74927,'amu*angstrom^2'), symmetry=1, barrier=(17.2272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.480944,'amu*angstrom^2'), symmetry=1, barrier=(11.0578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479396,'amu*angstrom^2'), symmetry=1, barrier=(11.0223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.047,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2934.99,'J/mol'), sigma=(5.41947,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.44 K, Pc=41.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.58094,0.0310091,0.000206561,-6.70643e-07,5.77215e-10,-121982,13.7981], Tmin=(10,'K'), Tmax=(429.405,'K')), NASAPolynomial(coeffs=[8.65797,0.0413465,-3.08669e-05,1.05258e-08,-1.33064e-12,-122950,-12.5905], Tmin=(429.405,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1014.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""CC(F)(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CC(F)(F)C(F)(F)OF(7713)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
11 C u0 p1 c0 {5,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-864.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,617,898,1187,245.785,245.791,245.805],'cm^-1')),
        HinderedRotor(inertia=(0.208192,'amu*angstrom^2'), symmetry=1, barrier=(8.93117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871229,'amu*angstrom^2'), symmetry=1, barrier=(37.3634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20854,'amu*angstrom^2'), symmetry=1, barrier=(8.93224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208351,'amu*angstrom^2'), symmetry=1, barrier=(8.93221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.908427,0.117217,-0.000195001,1.65422e-07,-5.51831e-11,-103839,33.2411], Tmin=(100,'K'), Tmax=(787.662,'K')), NASAPolynomial(coeffs=[14.884,0.0302124,-1.63525e-05,3.24692e-09,-2.28109e-13,-106116,-37.8387], Tmin=(787.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-864.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CJ2_singlet-FCs)"""),
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
    label = 'FCC(F)(F)C(F)(F)OF(4504)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1168.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,223,363,546,575,694,1179,1410,528,1116,1182,1331,1402,1494,3075,3110,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.484427,'amu*angstrom^2'), symmetry=1, barrier=(11.1379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485497,'amu*angstrom^2'), symmetry=1, barrier=(11.1625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.241,'amu*angstrom^2'), symmetry=1, barrier=(28.5331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.038,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.05,'J/mol'), sigma=(5.3108,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.67 K, Pc=44.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29587,0.0604691,2.34329e-05,-2.32208e-07,2.15507e-10,-140503,13.9768], Tmin=(10,'K'), Tmax=(472.08,'K')), NASAPolynomial(coeffs=[10.1056,0.039999,-2.98204e-05,1.0051e-08,-1.25463e-12,-141561,-18.1588], Tmin=(472.08,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1168.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""FCC(F)(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC(F)(CC(F)F)C(F)(F)F(7714)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {8,S}
8  O u0 p2 c0 {7,S} {9,S}
9  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
10 C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
12 C u0 p0 c0 {5,S} {6,S} {10,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1465.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (200.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17459,0.122933,-0.000177998,1.31591e-07,-3.87133e-11,-176132,33.1429], Tmin=(100,'K'), Tmax=(831.265,'K')), NASAPolynomial(coeffs=[17.1809,0.0346072,-1.86157e-05,3.76768e-09,-2.7089e-13,-179184,-52.0194], Tmin=(831.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1465.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCsFFH)"""),
)

species(
    label = 'FC(F)=C(F)CC(F)F(3719)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-999.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,323,467,575,827,1418,182,240,577,636,1210,1413,180,1397.7],'cm^-1')),
        HinderedRotor(inertia=(0.323852,'amu*angstrom^2'), symmetry=1, barrier=(7.44599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325281,'amu*angstrom^2'), symmetry=1, barrier=(7.47885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67908,0.0421484,0.00045811,-4.23698e-06,1.05977e-08,-120220,11.889], Tmin=(10,'K'), Tmax=(148.37,'K')), NASAPolynomial(coeffs=[5.21405,0.0457303,-3.26812e-05,1.08146e-08,-1.34136e-12,-120315,5.7446], Tmin=(148.37,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-999.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FC(F)DC(F)CC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)[C](F)CC(F)F(7715)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {10,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {11,S}
10 C u0 p0 c0 {1,S} {2,S} {8,S} {14,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1042.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,253,525,597,667,842,1178,1324,235,523,627,1123,1142,1372,1406,3097,212,367,445,1450,243.309,243.31,243.31,243.31],'cm^-1')),
        HinderedRotor(inertia=(0.196643,'amu*angstrom^2'), symmetry=1, barrier=(8.26084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196645,'amu*angstrom^2'), symmetry=1, barrier=(8.26084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196643,'amu*angstrom^2'), symmetry=1, barrier=(8.26084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755922,'amu*angstrom^2'), symmetry=1, barrier=(31.756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0179,0.121132,-0.000203307,1.77337e-07,-6.0852e-11,-125189,34.9654], Tmin=(100,'K'), Tmax=(796.779,'K')), NASAPolynomial(coeffs=[13.4891,0.0358662,-1.93725e-05,3.84776e-09,-2.70377e-13,-127106,-29.2493], Tmin=(796.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1042.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FO[C](F)C(F)(F)CC(F)F(7716)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {14,S}
11 C u1 p0 c0 {5,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1035.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,215,315,519,588,595,1205,1248,235,523,627,1123,1142,1372,1406,3097,395,473,707,1436,185.419,187.837,202.003,227.4],'cm^-1')),
        HinderedRotor(inertia=(1.33758,'amu*angstrom^2'), symmetry=1, barrier=(30.9407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380503,'amu*angstrom^2'), symmetry=1, barrier=(9.33356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403657,'amu*angstrom^2'), symmetry=1, barrier=(9.33735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403253,'amu*angstrom^2'), symmetry=1, barrier=(9.31267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24108,0.127251,-0.000219661,1.93588e-07,-6.64376e-11,-124361,34.1083], Tmin=(100,'K'), Tmax=(814.371,'K')), NASAPolynomial(coeffs=[13.9301,0.0356252,-1.93841e-05,3.83923e-09,-2.68479e-13,-126265,-32.4856], Tmin=(814.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1035.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[CH]CC(F)(F)C(F)(F)OF(7717)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
11 C u1 p0 c0 {5,S} {9,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1039.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,312.825,312.83,312.83,312.833],'cm^-1')),
        HinderedRotor(inertia=(0.122971,'amu*angstrom^2'), symmetry=1, barrier=(8.5403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122974,'amu*angstrom^2'), symmetry=1, barrier=(8.53956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122973,'amu*angstrom^2'), symmetry=1, barrier=(8.53989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493815,'amu*angstrom^2'), symmetry=1, barrier=(34.2938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19583,0.12505,-0.000211067,1.82656e-07,-6.19287e-11,-124803,34.1018], Tmin=(100,'K'), Tmax=(803.47,'K')), NASAPolynomial(coeffs=[14.5403,0.0343059,-1.84996e-05,3.66245e-09,-2.56539e-13,-126932,-35.8811], Tmin=(803.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1039.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[O]C(F)(F)C(F)(F)CC(F)F(7718)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {11,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1327.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,222,329,445,522,589,1214,1475,235,523,627,1123,1142,1372,1406,3097,351,323,533,609,664,892,1120,1201,210.846,210.869,210.872],'cm^-1')),
        HinderedRotor(inertia=(0.248904,'amu*angstrom^2'), symmetry=1, barrier=(7.85584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248874,'amu*angstrom^2'), symmetry=1, barrier=(7.85584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24888,'amu*angstrom^2'), symmetry=1, barrier=(7.85547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.513881,0.110246,-0.000162604,1.1706e-07,-2.92756e-11,-159484,31.7677], Tmin=(100,'K'), Tmax=(623.57,'K')), NASAPolynomial(coeffs=[13.8359,0.0341723,-1.80377e-05,3.58662e-09,-2.5376e-13,-161585,-33.1747], Tmin=(623.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1327.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'F[C](F)C(F)(F)CC(F)F(7719)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {13,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1180.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,215,315,519,588,595,1205,1248,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.643411,'amu*angstrom^2'), symmetry=1, barrier=(14.7933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363014,'amu*angstrom^2'), symmetry=1, barrier=(8.34641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363482,'amu*angstrom^2'), symmetry=1, barrier=(8.35717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.21026,0.0706899,-4.27245e-06,-2.36175e-07,2.69969e-10,-141936,15.2316], Tmin=(10,'K'), Tmax=(414.47,'K')), NASAPolynomial(coeffs=[9.66701,0.0432692,-3.13135e-05,1.0437e-08,-1.29873e-12,-142771,-13.8471], Tmin=(414.47,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1180.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""F[C](F)C(F)(F)CC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC(F)(F)[C](F)F(2831)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {6,S}
6 O u0 p2 c0 {5,S} {7,S}
7 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
"""),
    E0 = (-744.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.821163,'amu*angstrom^2'), symmetry=1, barrier=(18.8802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82472,'amu*angstrom^2'), symmetry=1, barrier=(18.9619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62228,0.0269041,0.0001746,-6.00025e-07,5.19722e-10,-89530.3,13.0155], Tmin=(10,'K'), Tmax=(443.604,'K')), NASAPolynomial(coeffs=[11.1479,0.0214546,-1.80038e-05,6.57499e-09,-8.70046e-13,-90812.1,-24.0954], Tmin=(443.604,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-744.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""FOC(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF2-CH2(64)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-300.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.294105,'amu*angstrom^2'), symmetry=1, barrier=(6.76206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.042,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95374,0.00280607,6.87135e-05,-1.351e-07,8.18166e-11,-36148.1,9.20378], Tmin=(10,'K'), Tmax=(530.258,'K')), NASAPolynomial(coeffs=[2.67787,0.0223019,-1.43606e-05,4.45233e-09,-5.29912e-13,-36151.6,13.2411], Tmin=(530.258,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-300.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[CH2]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FO[C](F)F(2563)',
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
    label = 'F[C](F)CC(F)F(2837)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-747.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,180,760.636],'cm^-1')),
        HinderedRotor(inertia=(0.184363,'amu*angstrom^2'), symmetry=1, barrier=(4.23888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184627,'amu*angstrom^2'), symmetry=1, barrier=(4.24494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.049,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2779.34,'J/mol'), sigma=(4.72663,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=434.13 K, Pc=59.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37926,0.0669634,-0.000254523,6.24431e-07,-5.64626e-10,-89891,13.2828], Tmin=(10,'K'), Tmax=(355.622,'K')), NASAPolynomial(coeffs=[4.07027,0.0357815,-2.42593e-05,7.6622e-09,-9.14224e-13,-89792.1,12.7447], Tmin=(355.622,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-747.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""F[C](F)CC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF2(82)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)(F)C(F)(F)OF(2832)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
9  C u1 p0 c0 {7,S} {10,S} {11,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-802.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.326112,'amu*angstrom^2'), symmetry=1, barrier=(7.49797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24298,'amu*angstrom^2'), symmetry=1, barrier=(28.5785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.326104,'amu*angstrom^2'), symmetry=1, barrier=(7.49776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39501,0.0506731,6.55634e-05,-3.33561e-07,3.0052e-10,-96536.9,15.1099], Tmin=(10,'K'), Tmax=(455.375,'K')), NASAPolynomial(coeffs=[9.9211,0.036171,-2.77248e-05,9.52105e-09,-1.20394e-12,-97575.2,-16.1162], Tmin=(455.375,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-802.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""[CH2]C(F)(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)C(F)(F)[CH]C(F)F(7720)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {8,S}
8  O u0 p2 c0 {7,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
12 C u1 p0 c0 {9,S} {11,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-1245.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,264.89,264.902,265.074,265.296],'cm^-1')),
        HinderedRotor(inertia=(0.143263,'amu*angstrom^2'), symmetry=1, barrier=(7.13297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143026,'amu*angstrom^2'), symmetry=1, barrier=(7.13275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.730587,'amu*angstrom^2'), symmetry=1, barrier=(36.3824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143328,'amu*angstrom^2'), symmetry=1, barrier=(7.1318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (199.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3154,0.127318,-0.000211766,1.80922e-07,-6.10507e-11,-149575,37.2754], Tmin=(100,'K'), Tmax=(773.426,'K')), NASAPolynomial(coeffs=[15.385,0.0344184,-1.89334e-05,3.79266e-09,-2.68115e-13,-151963,-37.7411], Tmin=(773.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1245.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = 'FOC(F)(F)C(F)(F)C[C](F)F(7721)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {8,S}
8  O u0 p2 c0 {7,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1252.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,223,363,546,575,694,1179,1410,190,488,555,1236,1407,193.665,193.665,193.665,193.665],'cm^-1')),
        HinderedRotor(inertia=(1.2018,'amu*angstrom^2'), symmetry=1, barrier=(31.9861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307564,'amu*angstrom^2'), symmetry=1, barrier=(8.18586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00449468,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307564,'amu*angstrom^2'), symmetry=1, barrier=(8.18586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (199.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52724,0.134659,-0.000236636,2.09887e-07,-7.21898e-11,-150420,36.2976], Tmin=(100,'K'), Tmax=(817.942,'K')), NASAPolynomial(coeffs=[14.5848,0.0362122,-2.00564e-05,3.98681e-09,-2.78913e-13,-152398,-34.177], Tmin=(817.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1252.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'O=C(F)C(F)(F)CC(F)F(7472)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {13,S}
10 C u0 p0 c0 {5,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1292.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,231,348,427,405,1245,1236,1280,235,523,627,1123,1142,1372,1406,3097,486,617,768,1157,1926,237.823,239.334,245.425],'cm^-1')),
        HinderedRotor(inertia=(0.174841,'amu*angstrom^2'), symmetry=1, barrier=(7.00378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172024,'amu*angstrom^2'), symmetry=1, barrier=(6.99403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454758,'amu*angstrom^2'), symmetry=1, barrier=(18.7266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135597,0.0926822,-0.000136042,1.08745e-07,-3.52944e-11,-155348,29.1389], Tmin=(100,'K'), Tmax=(750.792,'K')), NASAPolynomial(coeffs=[11.2672,0.0333718,-1.75374e-05,3.51059e-09,-2.50412e-13,-157019,-21.3731], Tmin=(750.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1292.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(COCsFO)"""),
)

species(
    label = 'FOC(F)=C(F)CC(F)F(7722)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-846.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,323,467,575,827,1418,326,540,652,719,1357,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.198229,'amu*angstrom^2'), symmetry=1, barrier=(4.55767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199031,'amu*angstrom^2'), symmetry=1, barrier=(4.57612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200806,'amu*angstrom^2'), symmetry=1, barrier=(4.61693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.202207,0.102706,-0.000171843,1.54012e-07,-5.43171e-11,-101685,30.6327], Tmin=(100,'K'), Tmax=(805.134,'K')), NASAPolynomial(coeffs=[10.0547,0.0362782,-1.92621e-05,3.80654e-09,-2.66868e-13,-102835,-13.5139], Tmin=(805.134,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-846.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FOC(F)(F)C(F)=CC(F)F(7723)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1101.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,171,205,598,1104,1143,1317,1411,3153,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.215032,'amu*angstrom^2'), symmetry=1, barrier=(4.94401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212997,'amu*angstrom^2'), symmetry=1, barrier=(4.89723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974307,'amu*angstrom^2'), symmetry=1, barrier=(22.4012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.61669,0.110619,-0.000180532,1.54055e-07,-5.22346e-11,-132318,31.6869], Tmin=(100,'K'), Tmax=(766.6,'K')), NASAPolynomial(coeffs=[13.3371,0.0325922,-1.76486e-05,3.52645e-09,-2.49395e-13,-134304,-30.9232], Tmin=(766.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1101.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH)"""),
)

species(
    label = 'FC=CC(F)(F)C(F)(F)OF(7724)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {8,S} {11,D} {12,S}
11 C u0 p0 c0 {5,S} {10,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1089.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,274,345,380,539,705,1166,1213,223,363,546,575,694,1179,1410,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,191.206,191.206,191.209],'cm^-1')),
        HinderedRotor(inertia=(0.232096,'amu*angstrom^2'), symmetry=1, barrier=(6.02148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232088,'amu*angstrom^2'), symmetry=1, barrier=(6.02145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0554,'amu*angstrom^2'), symmetry=1, barrier=(27.3815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.487175,0.107127,-0.000163929,1.29796e-07,-4.09582e-11,-130914,30.9393], Tmin=(100,'K'), Tmax=(776.223,'K')), NASAPolynomial(coeffs=[14.2973,0.0309432,-1.67163e-05,3.36643e-09,-2.40416e-13,-133209,-36.6429], Tmin=(776.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1089.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-CdsCsH) + group(CdCFH)"""),
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
    E0 = (-827.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-300.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-343.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-329.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-238.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-194.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-510.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-566.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-320.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-574.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-172.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-390.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-383.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-386.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-675.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-498.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-465.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-485.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-480.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-454.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-461.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-643.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-276.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-612.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-605.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    products = ['HF(38)', 'CF2O(49)', 'FC(F)=CC(F)F(2834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(46.3145,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'FOC(F)(F)CC(F)F(2803)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.22791e-07,'m^3/(mol*s)'), n=3.59921, Ea=(337.811,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node CC_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FOC(F)(F)F(2825)', 'F[C]CC(F)F(718)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.05141e+54,'m^3/(mol*s)'), n=-13.541, Ea=(187.408,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]C(F)(F)OF(2807)', 'CH2FCHF2(1)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(178.725,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(S)(25)', 'FOC(F)(F)C(F)(F)C(F)F(7711)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000166864,'m^3/(mol*s)'), n=2.82973, Ea=(135.162,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FOF(488)', 'F[C]C(F)(F)CC(F)F(7712)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(13.5225,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CF2(43)', 'CC(F)(F)C(F)(F)OF(2804)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.00405e-10,'m^3/(mol*s)'), n=4.72997, Ea=(127.877,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HF(38)', 'F[C]CC(F)(F)C(F)(F)OF(7713)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.9697e-08,'m^3/(mol*s)'), n=3.78805, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3001410487914693, var=2.9282398409420045, Tref=1000.0, N=3, data_mean=0.0, correlation='HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R',), comment="""Estimated from node HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CHF(40)', 'FCC(F)(F)C(F)(F)OF(4504)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96165e-06,'m^3/(mol*s)'), n=3.30609, Ea=(129.274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s',), comment="""Estimated from node CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    products = ['FOC(F)(CC(F)F)C(F)(F)F(7714)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction11',
    reactants = ['FOF(488)', 'FC(F)=C(F)CC(F)F(3719)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'FOC(F)(F)[C](F)CC(F)F(7715)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'FO[C](F)C(F)(F)CC(F)F(7716)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[CH]CC(F)(F)C(F)(F)OF(7717)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(F)(F)C(F)(F)CC(F)F(7718)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]F(128)', 'F[C](F)C(F)(F)CC(F)F(7719)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FOC(F)(F)[C](F)F(2831)', 'CHF2-CH2(64)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -19.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['FO[C](F)F(2563)', 'F[C](F)CC(F)F(2837)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -13.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CHF2(82)', '[CH2]C(F)(F)C(F)(F)OF(2832)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -15.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'FOC(F)(F)C(F)(F)[CH]C(F)F(7720)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', 'FOC(F)(F)C(F)(F)C[C](F)F(7721)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['F2(78)', 'O=C(F)C(F)(F)CC(F)F(7472)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(78.7711,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F2(78)', 'FOC(F)=C(F)CC(F)F(7722)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'FOC(F)(F)C(F)=CC(F)F(7723)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(191.21,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'FC=CC(F)(F)C(F)(F)OF(7724)'],
    products = ['FOC(F)(F)C(F)(F)CC(F)F(7684)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(186.03,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

network(
    label = 'PDepNetwork #2067',
    isomers = [
        'FOC(F)(F)C(F)(F)CC(F)F(7684)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(49)', 'FC(F)=CC(F)F(2834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2067',
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

