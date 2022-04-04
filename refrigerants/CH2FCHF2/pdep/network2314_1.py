species(
    label = 'C=C([O])C(F)(F)[CH]C(F)F(7442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {9,D} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-739.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.258502,'amu*angstrom^2'), symmetry=1, barrier=(5.94347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25733,'amu*angstrom^2'), symmetry=1, barrier=(5.91652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258186,'amu*angstrom^2'), symmetry=1, barrier=(5.9362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.397239,0.107548,-0.000179588,1.6087e-07,-5.65863e-11,-88797.3,32.7089], Tmin=(100,'K'), Tmax=(811.585,'K')), NASAPolynomial(coeffs=[10.2017,0.0382674,-2.0042e-05,3.93941e-09,-2.75257e-13,-89956.4,-12.7544], Tmin=(811.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-739.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'C=C([O])C([C](F)F)C(F)F(7423)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {8,D} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-743.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,374.813,374.914,374.932,374.943],'cm^-1')),
        HinderedRotor(inertia=(0.118405,'amu*angstrom^2'), symmetry=1, barrier=(11.8111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118439,'amu*angstrom^2'), symmetry=1, barrier=(11.8119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118425,'amu*angstrom^2'), symmetry=1, barrier=(11.8119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3708.67,'J/mol'), sigma=(5.91702,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=579.29 K, Pc=40.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232459,0.0970303,-0.000127125,8.50037e-08,-2.24778e-11,-89256.6,33.2258], Tmin=(100,'K'), Tmax=(925.568,'K')), NASAPolynomial(coeffs=[16.1358,0.0262894,-1.24765e-05,2.42139e-09,-1.71059e-13,-92286.4,-44.4749], Tmin=(925.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-743.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = '[O]C(CF)=C(F)[CH]C(F)F(7441)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
7  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
8  C u1 p0 c0 {6,S} {10,S} {14,S}
9  C u0 p0 c0 {5,S} {7,S} {10,D}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-735.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,212,931,1104,1251,1325,1474,3102,3155,3025,407.5,1350,352.5,350,440,435,1725,271,519,563,612,1379,208.742,209.765,211.317,1489.59],'cm^-1')),
        HinderedRotor(inertia=(0.136137,'amu*angstrom^2'), symmetry=1, barrier=(4.37102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140941,'amu*angstrom^2'), symmetry=1, barrier=(4.37351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44826,'amu*angstrom^2'), symmetry=1, barrier=(45.529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.369146,0.0899309,-0.000115395,7.22817e-08,-1.22469e-11,-88398.4,28.9959], Tmin=(100,'K'), Tmax=(584.882,'K')), NASAPolynomial(coeffs=[9.82856,0.0390182,-2.01639e-05,4.01666e-09,-2.86047e-13,-89740.6,-13.5817], Tmin=(584.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-735.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
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
    label = '[CH2]C(=O)[C](F)F(6453)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {3,D} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {2,S} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-279.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,179,346,818,1406,1524,241.327],'cm^-1')),
        HinderedRotor(inertia=(0.0200138,'amu*angstrom^2'), symmetry=1, barrier=(12.2066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28776,'amu*angstrom^2'), symmetry=1, barrier=(53.2065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00969,0.0471598,-6.16119e-05,4.28745e-08,-1.20388e-11,-33525,18.4741], Tmin=(100,'K'), Tmax=(865.239,'K')), NASAPolynomial(coeffs=[8.69208,0.016267,-8.05515e-06,1.60882e-09,-1.15551e-13,-34681.4,-12.7972], Tmin=(865.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(CO-CsO2d)(F1s)(F1s)) + radical(C2JC=O)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C(F)(F)[CH]C(F)F(8838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {9,D} {12,S} {13,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-441.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.503289,'amu*angstrom^2'), symmetry=1, barrier=(11.5716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247276,'amu*angstrom^2'), symmetry=1, barrier=(5.68537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24896,'amu*angstrom^2'), symmetry=1, barrier=(5.72407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08451,0.0768917,-7.47449e-05,-1.49672e-08,5.80788e-11,-52972.5,28.268], Tmin=(100,'K'), Tmax=(486.377,'K')), NASAPolynomial(coeffs=[8.71239,0.0365567,-1.94242e-05,3.88384e-09,-2.75747e-13,-53979.4,-5.75752], Tmin=(486.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH))"""),
)

species(
    label = 'C=C1OC(C(F)F)C1(F)F(8663)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {10,D}
10 C u0 p0 c0 {9,D} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-929.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.947659,0.0893771,-9.22731e-05,4.54534e-08,-8.43322e-12,-111608,26.0809], Tmin=(100,'K'), Tmax=(1471.85,'K')), NASAPolynomial(coeffs=[24.744,0.0101543,-1.95398e-06,2.0412e-10,-1.03172e-14,-118153,-104.338], Tmin=(1471.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-929.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=C(O)C(F)(F)C=C(F)F(8668)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {14,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,D}
8  C u0 p0 c0 {6,S} {10,D} {11,S}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-927.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.305592,0.101573,-0.000143973,1.0654e-07,-3.14724e-11,-111429,28.285], Tmin=(100,'K'), Tmax=(827.663,'K')), NASAPolynomial(coeffs=[14.3619,0.0306873,-1.5507e-05,3.06424e-09,-2.17458e-13,-113857,-39.703], Tmin=(827.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFF)"""),
)

species(
    label = 'FC(F)[CH]C(F)(F)[C]1CO1(8839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 C u1 p0 c0 {6,S} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-615.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2629,0.0741477,-4.50102e-05,-8.46304e-08,1.16378e-10,-73893.4,29.5583], Tmin=(100,'K'), Tmax=(455.196,'K')), NASAPolynomial(coeffs=[7.91276,0.0413257,-2.1255e-05,4.19215e-09,-2.95251e-13,-74764.2,-0.204568], Tmin=(455.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(C-FF)-Cs-O2s) + radical(C2CsJO) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H))"""),
)

species(
    label = '[O][C]1CC(C(F)F)C1(F)F(8840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {14,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-634.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469005,0.0774559,-7.36217e-05,3.54553e-08,-6.9075e-12,-76226.5,28.5423], Tmin=(100,'K'), Tmax=(1220.13,'K')), NASAPolynomial(coeffs=[15.3561,0.0286513,-1.3623e-05,2.67288e-09,-1.90556e-13,-79859.4,-46.2414], Tmin=(1220.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-634.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsFFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1([O])C(C(F)F)C1(F)F(8822)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
10 C u1 p0 c0 {7,S} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-538.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.112027,0.0952514,-0.000117913,7.56286e-08,-1.9451e-11,-64677,28.2699], Tmin=(100,'K'), Tmax=(944.227,'K')), NASAPolynomial(coeffs=[15.1804,0.0304679,-1.49966e-05,2.96477e-09,-2.11885e-13,-67564.9,-44.6296], Tmin=(944.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(CsCsFFH) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = 'C=C([O])C(F)=CC(F)F(8841)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {3,S} {6,D} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-607.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,280,518,736,852,873,350,440,435,1725,2950,3100,1380,975,1025,1650,294.058,294.09,294.105],'cm^-1')),
        HinderedRotor(inertia=(0.180267,'amu*angstrom^2'), symmetry=1, barrier=(11.0544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180171,'amu*angstrom^2'), symmetry=1, barrier=(11.0536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538916,0.0800497,-9.60988e-05,6.05224e-08,-1.5363e-11,-72913.8,26.3156], Tmin=(100,'K'), Tmax=(954.185,'K')), NASAPolynomial(coeffs=[13.1209,0.0273051,-1.3183e-05,2.59098e-09,-1.84711e-13,-75314.9,-33.7952], Tmin=(954.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-607.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(4373)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C([O])C(F)(F)C=CF(8842)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {6,D} {11,S}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-588.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25517,'amu*angstrom^2'), symmetry=1, barrier=(28.8589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25792,'amu*angstrom^2'), symmetry=1, barrier=(28.9221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190375,0.091117,-0.00014099,1.19587e-07,-4.0553e-11,-70695.3,26.9822], Tmin=(100,'K'), Tmax=(793.493,'K')), NASAPolynomial(coeffs=[10.3998,0.0323665,-1.61585e-05,3.13727e-09,-2.18663e-13,-72086.2,-18.4655], Tmin=(793.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-588.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])C(F)(F)C=C(F)F(8843)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u0 p0 c0 {3,S} {4,S} {7,D}
10 C u0 p0 c0 {8,D} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-789.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,182,240,577,636,1210,1413,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25419,'amu*angstrom^2'), symmetry=1, barrier=(28.8364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26647,'amu*angstrom^2'), symmetry=1, barrier=(29.1186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.176622,0.100682,-0.000163743,1.41933e-07,-4.85417e-11,-94862.1,29.3051], Tmin=(100,'K'), Tmax=(808.053,'K')), NASAPolynomial(coeffs=[11.1702,0.0330655,-1.69762e-05,3.31355e-09,-2.30834e-13,-96322.1,-20.7061], Tmin=(808.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-789.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFF) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])C(F)(F)C[C](F)F(8844)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u0 p0 c0 {8,D} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-733.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.185055,'amu*angstrom^2'), symmetry=1, barrier=(4.25479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185171,'amu*angstrom^2'), symmetry=1, barrier=(4.25746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18494,'amu*angstrom^2'), symmetry=1, barrier=(4.25213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.535247,0.111502,-0.000190337,1.71157e-07,-5.97773e-11,-88071.7,32.01], Tmin=(100,'K'), Tmax=(829.766,'K')), NASAPolynomial(coeffs=[10.4256,0.037874,-1.96542e-05,3.82819e-09,-2.65295e-13,-89174.9,-14.5114], Tmin=(829.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH]=C(O)C(F)(F)[CH]C(F)F(8845)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {9,S} {13,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u1 p0 c0 {9,D} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-630.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.57536,'amu*angstrom^2'), symmetry=1, barrier=(13.2287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574301,'amu*angstrom^2'), symmetry=1, barrier=(13.2043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575176,'amu*angstrom^2'), symmetry=1, barrier=(13.2244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575817,'amu*angstrom^2'), symmetry=1, barrier=(13.2392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.796314,0.114503,-0.000186525,1.57793e-07,-5.278e-11,-75636.4,33.0469], Tmin=(100,'K'), Tmax=(782.018,'K')), NASAPolynomial(coeffs=[14.1025,0.0322194,-1.70392e-05,3.37034e-09,-2.36788e-13,-77780.8,-33.9802], Tmin=(782.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-630.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(F)(F)CC(F)F(8846)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {13,S}
9  C u0 p0 c0 {5,S} {7,S} {10,D}
10 C u1 p0 c0 {9,D} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-686.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,274,345,380,539,705,1166,1213,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,3120,650,792.5,1650,235.002,235.172,236.033],'cm^-1')),
        HinderedRotor(inertia=(0.241625,'amu*angstrom^2'), symmetry=1, barrier=(9.39495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2416,'amu*angstrom^2'), symmetry=1, barrier=(9.40155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243063,'amu*angstrom^2'), symmetry=1, barrier=(9.40112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.667425,0.112318,-0.00018472,1.59423e-07,-5.40514e-11,-82465.6,31.5275], Tmin=(100,'K'), Tmax=(812.73,'K')), NASAPolynomial(coeffs=[12.6297,0.0343472,-1.76947e-05,3.4509e-09,-2.39985e-13,-84213.3,-27.3209], Tmin=(812.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-686.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C(F)(F)[CH][C](F)F(8847)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {14,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,D}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-676.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,248,333,466,604,684,796,1061,1199,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,190,488,555,1236,1407,217.486,217.49],'cm^-1')),
        HinderedRotor(inertia=(0.179825,'amu*angstrom^2'), symmetry=1, barrier=(6.03615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179825,'amu*angstrom^2'), symmetry=1, barrier=(6.03615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179831,'amu*angstrom^2'), symmetry=1, barrier=(6.03617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470195,'amu*angstrom^2'), symmetry=1, barrier=(15.7833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.68329,0.113935,-0.000193116,1.70945e-07,-5.91788e-11,-81241.6,33.5967], Tmin=(100,'K'), Tmax=(812.916,'K')), NASAPolynomial(coeffs=[11.955,0.0356438,-1.8937e-05,3.7326e-09,-2.60821e-13,-82764.3,-21.4851], Tmin=(812.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-676.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(OF)=C(F)[CH]C(F)F(8848)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {6,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u1 p0 c0 {8,S} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-401.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,350,440,435,1725,271,519,563,612,1379,3000,3100,440,815,1455,1000,229.831,229.831,1971.75],'cm^-1')),
        HinderedRotor(inertia=(1.06812,'amu*angstrom^2'), symmetry=1, barrier=(40.0377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00239497,'amu*angstrom^2'), symmetry=1, barrier=(6.60741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176273,'amu*angstrom^2'), symmetry=1, barrier=(6.60741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06814,'amu*angstrom^2'), symmetry=1, barrier=(40.0377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.197909,0.100562,-0.000154404,1.31264e-07,-4.50169e-11,-48123.9,31.4889], Tmin=(100,'K'), Tmax=(769.784,'K')), NASAPolynomial(coeffs=[10.8284,0.0371755,-1.90214e-05,3.73821e-09,-2.62759e-13,-49641,-17.6495], Tmin=(769.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Allyl_S) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])[C](F)C(F)C(F)F(8849)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {12,S}
8  C u1 p0 c0 {4,S} {6,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,D}
10 C u0 p0 c0 {9,D} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-747.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152594,0.0949111,-0.000116414,7.29056e-08,-1.81783e-11,-89764,30.6127], Tmin=(100,'K'), Tmax=(976.107,'K')), NASAPolynomial(coeffs=[16.1775,0.0279931,-1.3582e-05,2.67405e-09,-1.90907e-13,-92952.1,-47.7758], Tmin=(976.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-747.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'C=C(OF)C(F)(F)[CH][CH]F(8850)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,D}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  C u0 p0 c0 {7,D} {13,S} {14,S}
10 C u1 p0 c0 {3,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-305.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,248,333,466,604,684,796,1061,1199,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,334,575,1197,1424,3202,180,1921.98,1923.06],'cm^-1')),
        HinderedRotor(inertia=(0.336745,'amu*angstrom^2'), symmetry=1, barrier=(7.74243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59499,'amu*angstrom^2'), symmetry=1, barrier=(36.6719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336154,'amu*angstrom^2'), symmetry=1, barrier=(7.72885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336244,'amu*angstrom^2'), symmetry=1, barrier=(7.73091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.562425,0.117145,-0.000217166,2.08842e-07,-7.63522e-11,-36645.5,34.6185], Tmin=(100,'K'), Tmax=(837.608,'K')), NASAPolynomial(coeffs=[7.07233,0.045129,-2.45242e-05,4.83586e-09,-3.3625e-13,-36677.2,6.58373], Tmin=(837.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-305.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'C=C([O])C(F)(F)C(F)[CH]F(8851)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 C u0 p0 c0 {8,D} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-698.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,259,529,569,1128,1321,1390,3140,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,228.385,229.339,229.43,229.684],'cm^-1')),
        HinderedRotor(inertia=(0.61651,'amu*angstrom^2'), symmetry=1, barrier=(22.8653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173541,'amu*angstrom^2'), symmetry=1, barrier=(6.34491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170397,'amu*angstrom^2'), symmetry=1, barrier=(6.33877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750925,0.114121,-0.000187901,1.61399e-07,-5.43619e-11,-83792.1,31.8986], Tmin=(100,'K'), Tmax=(816.616,'K')), NASAPolynomial(coeffs=[13.094,0.0337828,-1.73296e-05,3.36949e-09,-2.33778e-13,-85635.7,-29.5337], Tmin=(816.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'O=[C]CC(F)(F)[CH]C(F)F(8660)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-737.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.24131,'amu*angstrom^2'), symmetry=1, barrier=(5.5482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241262,'amu*angstrom^2'), symmetry=1, barrier=(5.54709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241348,'amu*angstrom^2'), symmetry=1, barrier=(5.54907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241364,'amu*angstrom^2'), symmetry=1, barrier=(5.54943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3459.05,'J/mol'), sigma=(5.90821,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=540.30 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.385676,0.107693,-0.000182207,1.63833e-07,-5.74391e-11,-88602.2,34.9153], Tmin=(100,'K'), Tmax=(823.858,'K')), NASAPolynomial(coeffs=[10.0725,0.0376396,-1.95645e-05,3.82149e-09,-2.65579e-13,-89671.3,-9.54284], Tmin=(823.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(CCCJ=O)"""),
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
    label = 'O=[C]C(F)(F)[CH]C(F)F(7611)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u1 p0 c0 {5,D} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-679.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([169,350,401,508,725,1243,1223,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1855,455,950,282.999,1970.27],'cm^-1')),
        HinderedRotor(inertia=(0.161933,'amu*angstrom^2'), symmetry=1, barrier=(9.19949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162933,'amu*angstrom^2'), symmetry=1, barrier=(9.20114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415269,'amu*angstrom^2'), symmetry=1, barrier=(23.5268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962636,0.0729463,-9.82022e-05,7.04273e-08,-2.05498e-11,-81656.4,29.0737], Tmin=(100,'K'), Tmax=(830.514,'K')), NASAPolynomial(coeffs=[10.642,0.026328,-1.40048e-05,2.84085e-09,-2.05118e-13,-83264.2,-15.8258], Tmin=(830.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(COj(Cs-F1sF1sCs)(O2d))"""),
)

species(
    label = 'O=C1CC(C(F)F)C1(F)F(8763)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {14,S}
10 C u0 p0 c0 {5,D} {7,S} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-953.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587388,0.0698829,-5.75393e-05,2.33892e-08,-3.81526e-12,-114507,28.2193], Tmin=(100,'K'), Tmax=(1445.97,'K')), NASAPolynomial(coeffs=[16.5882,0.0256194,-1.16215e-05,2.21862e-09,-1.54955e-13,-119134,-54.8762], Tmin=(1445.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-953.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-OdCsCs) + longDistanceInteraction_cyclic(Cs(F)2-CO) + ring(Cyclobutanone)"""),
)

species(
    label = 'CC(=O)C(F)(F)C=C(F)F(8669)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {5,D} {6,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-949.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244291,0.0899179,-0.000125988,1.00157e-07,-3.27822e-11,-114063,29.3574], Tmin=(100,'K'), Tmax=(741.566,'K')), NASAPolynomial(coeffs=[10.0765,0.0368801,-1.87001e-05,3.70033e-09,-2.62629e-13,-115521,-15.1372], Tmin=(741.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-949.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(CdCFF)"""),
)

species(
    label = '[CH2][C]1OC(C(F)F)C1(F)F(8852)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 C u1 p0 c0 {9,S} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-612.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00879514,0.0951289,-0.000133257,1.02238e-07,-3.18924e-11,-73558.7,29.1057], Tmin=(100,'K'), Tmax=(780.085,'K')), NASAPolynomial(coeffs=[11.7229,0.0350584,-1.77397e-05,3.50849e-09,-2.49048e-13,-75386.1,-24.4977], Tmin=(780.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-612.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCsFFH) + group(Cs-CsHHH) + ring(Cs-O2s-Cs-Cs(F)(F)) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = 'CC(=O)C(F)(F)[CH][C](F)F(8853)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {5,D} {6,S} {7,S}
9  C u1 p0 c0 {6,S} {10,S} {14,S}
10 C u1 p0 c0 {3,S} {4,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-693.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,190,488,555,1236,1407,337.076,341.313,2145.06],'cm^-1')),
        HinderedRotor(inertia=(0.00148405,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0917235,'amu*angstrom^2'), symmetry=1, barrier=(7.56905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0926173,'amu*angstrom^2'), symmetry=1, barrier=(7.55831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0940549,'amu*angstrom^2'), symmetry=1, barrier=(7.59121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0350972,0.0970503,-0.000150164,1.29393e-07,-4.49856e-11,-83231.2,34.264], Tmin=(100,'K'), Tmax=(770.013,'K')), NASAPolynomial(coeffs=[10.16,0.0369121,-1.90313e-05,3.75524e-09,-2.64554e-13,-84588.5,-10.8752], Tmin=(770.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-693.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(CF)C(F)(F)[CH][CH]F(8854)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {5,D} {6,S} {7,S}
9  C u1 p0 c0 {6,S} {10,S} {13,S}
10 C u1 p0 c0 {4,S} {9,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-643.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,244,1102,1254,1379,1468,1475,3037,3083,375,552.5,462.5,1710,3025,407.5,1350,352.5,334,575,1197,1424,3202,314.875,314.884,1700.7,1700.7],'cm^-1')),
        HinderedRotor(inertia=(0.150762,'amu*angstrom^2'), symmetry=1, barrier=(10.6073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(10.6074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150763,'amu*angstrom^2'), symmetry=1, barrier=(10.6074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364187,'amu*angstrom^2'), symmetry=1, barrier=(25.6242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2935,0.103363,-0.000164716,1.42552e-07,-4.92481e-11,-77304.7,35.0477], Tmin=(100,'K'), Tmax=(783.748,'K')), NASAPolynomial(coeffs=[11.0973,0.0361783,-1.88116e-05,3.70982e-09,-2.607e-13,-78812.2,-15.3575], Tmin=(783.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    E0 = (-257.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-100.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (53.1086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (133.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (283.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-249.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-232.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-26.6096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-132.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-57.2465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (5.01654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-117.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (18.0179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-94.5562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-130.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (87.4976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-71.3275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (43.5944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-160.963,'kJ/mol'),
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
    E0 = (174.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-79.9336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (223.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-54.6397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-100.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (183.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-249.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-232.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-131.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-155.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-8.23922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['CH2CO(28)', 'FC(F)=CC(F)F(2834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['C=C([O])C([C](F)F)C(F)F(7423)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(CF)=C(F)[CH]C(F)F(7441)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(307.396,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(F)F-2(3493)', '[CH2]C(=O)[C](F)F(6453)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'C=[C]C(F)(F)[CH]C(F)F(8838)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['C=C1OC(C(F)F)C1(F)F(8663)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['C=C(O)C(F)(F)C=C(F)F(8668)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['FC(F)[CH]C(F)(F)[C]1CO1(8839)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['[O][C]1CC(C(F)F)C1(F)F(8840)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['[CH2]C1([O])C(C(F)F)C1(F)F(8822)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(200.579,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 200.0 to 200.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C=C([O])C(F)=CC(F)F(8841)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(57.6677,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(4373)', 'FC(F)=CC(F)F(2834)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1256.55,'m^3/(mol*s)'), n=0.85398, Ea=(33.1913,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.03938838252749272, var=0.8044325964401412, Tref=1000.0, N=80, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_Ext-2C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C=C([O])C(F)(F)C=CF(8842)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(52.2966,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'C=C([O])C(F)(F)C=C(F)F(8843)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.84061,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2CO(28)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.13223e-11,'m^3/(mol*s)'), n=4.48095, Ea=(3.08987,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=O(4373)', 'F[C](F)[CH]C(F)F(4845)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C(F)(F)C[C](F)F(8844)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.76836e+09,'s^-1'), n=1.1815, Ea=(180.499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(O)C(F)(F)[CH]C(F)F(8845)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C([O])C(F)(F)CC(F)F(8846)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C(F)(F)[CH][C](F)F(8847)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(57774.9,'s^-1'), n=1.76812, Ea=(63.6725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;XH_out] for rate rule [R5HJ_1;C_rad_out_noH;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(OF)=C(F)[CH]C(F)F(8848)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(94.0779,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['C=C([O])[C](F)C(F)C(F)F(8849)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(177.892,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(OF)C(F)(F)[CH][CH]F(8850)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(47.904,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C(F)(F)C(F)[CH]F(8851)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(161.686,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['O=[C]CC(F)(F)[CH]C(F)F(8660)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(18)', 'O=[C]C(F)(F)[CH]C(F)F(7611)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['O=C1CC(C(F)F)C1(F)F(8763)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['CC(=O)C(F)(F)C=C(F)F(8669)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    products = ['[CH2][C]1OC(C(F)F)C1(F)F(8852)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(126.791,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 125.7 to 126.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC(=O)C(F)(F)[CH][C](F)F(8853)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(102311,'s^-1'), n=1.78229, Ea=(55.9436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=C(CF)C(F)(F)[CH][CH]F(8854)'],
    products = ['C=C([O])C(F)(F)[CH]C(F)F(7442)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(154.016,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #2314',
    isomers = [
        'C=C([O])C(F)(F)[CH]C(F)F(7442)',
    ],
    reactants = [
        ('CH2CO(28)', 'FC(F)=CC(F)F(2834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2314',
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

