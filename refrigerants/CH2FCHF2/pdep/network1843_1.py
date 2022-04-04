species(
    label = 'O=CC(F)OOC(F)=CF(6221)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-727.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,326,540,652,719,1357,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.262554,'amu*angstrom^2'), symmetry=1, barrier=(6.03662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03926,'amu*angstrom^2'), symmetry=1, barrier=(46.8867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262789,'amu*angstrom^2'), symmetry=1, barrier=(6.04204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03387,'amu*angstrom^2'), symmetry=1, barrier=(46.7627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889086,0.0826447,-7.83955e-05,-2.54422e-08,7.22865e-11,-87350.9,29.2258], Tmin=(100,'K'), Tmax=(478.249,'K')), NASAPolynomial(coeffs=[9.0264,0.039541,-2.1475e-05,4.32109e-09,-3.0767e-13,-88414.7,-7.01355], Tmin=(478.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-727.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=CC(=O)F(335)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC=C(F)OOCF(3236)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-608.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,326,540,652,719,1357,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.226993,'amu*angstrom^2'), symmetry=1, barrier=(5.21902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226976,'amu*angstrom^2'), symmetry=1, barrier=(5.21861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24262,'amu*angstrom^2'), symmetry=1, barrier=(51.5623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3151.58,'J/mol'), sigma=(5.06004,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.27 K, Pc=55.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06891,0.072491,-0.000116297,1.08471e-07,-4.06571e-11,-73105,25.8181], Tmin=(100,'K'), Tmax=(777.321,'K')), NASAPolynomial(coeffs=[6.04495,0.034313,-1.83642e-05,3.6723e-09,-2.60394e-13,-73498.8,5.50817], Tmin=(777.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-608.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=C[C]F-2(1596)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OOC(F)=CF(212)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {8,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-420.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.525346,'amu*angstrom^2'), symmetry=1, barrier=(12.0787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13913,'amu*angstrom^2'), symmetry=1, barrier=(26.1909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73442,0.020817,0.000110574,-3.93802e-07,3.62983e-10,-50623.8,10.4607], Tmin=(10,'K'), Tmax=(403.043,'K')), NASAPolynomial(coeffs=[6.33783,0.0252554,-1.86223e-05,6.27928e-09,-7.87518e-13,-51079.6,-2.78415], Tmin=(403.043,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-420.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OOC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C(F)OOC=COF(6948)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {3,S} {8,S}
7  C u0 p0 c0 {4,S} {8,D} {11,S}
8  C u0 p0 c0 {6,S} {7,D} {12,S}
9  C u0 p0 c0 {1,S} {5,S} {10,D}
10 C u0 p0 c0 {2,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-339.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,194,682,905,1196,1383,3221,239.458,243.085],'cm^-1')),
        HinderedRotor(inertia=(0.208805,'amu*angstrom^2'), symmetry=1, barrier=(8.68646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13543,'amu*angstrom^2'), symmetry=1, barrier=(46.5414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210402,'amu*angstrom^2'), symmetry=1, barrier=(8.67701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55909,'amu*angstrom^2'), symmetry=1, barrier=(22.9289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0709549,0.0959291,-0.000132292,9.4005e-08,-2.66518e-11,-40741.6,33.1662], Tmin=(100,'K'), Tmax=(861.087,'K')), NASAPolynomial(coeffs=[14.4666,0.0283981,-1.46544e-05,2.92842e-09,-2.0955e-13,-43245.2,-34.7949], Tmin=(861.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FC=COOOC(F)=CF(6949)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {4,S} {5,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u0 p0 c0 {2,S} {7,D} {12,S}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-417.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,326,540,652,719,1357,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590762,0.0822872,-0.000103559,7.12412e-08,-2.02814e-11,-50151.5,34.0473], Tmin=(100,'K'), Tmax=(844.288,'K')), NASAPolynomial(coeffs=[10.7629,0.0340946,-1.79391e-05,3.63424e-09,-2.62661e-13,-51869.2,-13.306], Tmin=(844.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=CC(F)OC(F)C(=O)F(6950)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {10,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p0 c0 {3,S} {5,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1117.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479607,0.0849165,-0.000106774,7.24963e-08,-2.0296e-11,-134325,29.2373], Tmin=(100,'K'), Tmax=(858.913,'K')), NASAPolynomial(coeffs=[11.3809,0.0341489,-1.81143e-05,3.68175e-09,-2.66614e-13,-136198,-21.6975], Tmin=(858.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1117.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COCsFO)"""),
)

species(
    label = 'F[C]1OOC(F)[CH]OC1F(6951)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
9  C u1 p0 c0 {4,S} {8,S} {13,S}
10 C u1 p0 c0 {3,S} {6,S} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-548.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([401,573,564,712,601,775,1078,1160,1259,1391,1337,1437,3048,3250,2950,1000,395,473,707,1436,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.17511,0.0474193,9.09613e-05,-1.97376e-07,9.18339e-11,-65797,28.9714], Tmin=(100,'K'), Tmax=(912.199,'K')), NASAPolynomial(coeffs=[42.8848,-0.0246534,1.75026e-05,-3.38918e-09,2.17835e-13,-78510.1,-201.434], Tmin=(912.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFHO) + ring(Cycloheptane) + radical(CCsJOCs) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C[C](F)OOC(F)=CF(6952)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-397.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180,928.394],'cm^-1')),
        HinderedRotor(inertia=(0.228464,'amu*angstrom^2'), symmetry=1, barrier=(5.25285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227672,'amu*angstrom^2'), symmetry=1, barrier=(5.23464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228848,'amu*angstrom^2'), symmetry=1, barrier=(5.26166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227664,'amu*angstrom^2'), symmetry=1, barrier=(5.23443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.278238,0.114565,-0.000226109,2.27853e-07,-8.56048e-11,-47628.9,35.4142], Tmin=(100,'K'), Tmax=(845.326,'K')), NASAPolynomial(coeffs=[3.52284,0.0495932,-2.74454e-05,5.42316e-09,-3.7653e-13,-46592.8,27.6444], Tmin=(845.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-397.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OOC(F)[CH]F(6351)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {10,S}
9  C u1 p0 c0 {3,S} {7,S} {12,S}
10 C u0 p0 c0 {6,D} {8,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-579.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,280,501,1494,1531,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.81178,'amu*angstrom^2'), symmetry=1, barrier=(41.6564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8152,'amu*angstrom^2'), symmetry=1, barrier=(41.7349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80504,'amu*angstrom^2'), symmetry=1, barrier=(41.5013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63671,'amu*angstrom^2'), symmetry=1, barrier=(37.6312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81886,'amu*angstrom^2'), symmetry=1, barrier=(41.8191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0569599,0.097903,-0.000143393,1.16053e-07,-3.84935e-11,-69611.8,31.621], Tmin=(100,'K'), Tmax=(732.585,'K')), NASAPolynomial(coeffs=[10.9457,0.0378278,-2.03875e-05,4.11742e-09,-2.952e-13,-71224,-18.0371], Tmin=(732.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OO[C](F)CF(6215)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-584.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.93823,'amu*angstrom^2'), symmetry=1, barrier=(21.5717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95458,'amu*angstrom^2'), symmetry=1, barrier=(44.9397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95458,'amu*angstrom^2'), symmetry=1, barrier=(44.9397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95459,'amu*angstrom^2'), symmetry=1, barrier=(44.9399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95472,'amu*angstrom^2'), symmetry=1, barrier=(44.9428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3609.11,'J/mol'), sigma=(5.7886,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.73 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823099,0.0849361,-7.65299e-05,-4.17222e-08,8.93841e-11,-70142.6,29.1265], Tmin=(100,'K'), Tmax=(472.756,'K')), NASAPolynomial(coeffs=[9.40092,0.0399798,-2.15262e-05,4.29658e-09,-3.03843e-13,-71262.3,-9.09473], Tmin=(472.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(CsCsF1sO2s) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OOC(F)[CH]F(6915)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {10,S} {12,S}
9  C u1 p0 c0 {3,S} {7,S} {13,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-560.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,355,410,600,1181,1341,1420,3056,334,575,1197,1424,3202,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.929087,'amu*angstrom^2'), symmetry=1, barrier=(21.3615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931461,'amu*angstrom^2'), symmetry=1, barrier=(21.4161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.923538,'amu*angstrom^2'), symmetry=1, barrier=(21.234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73244,'amu*angstrom^2'), symmetry=1, barrier=(39.8322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73748,'amu*angstrom^2'), symmetry=1, barrier=(39.9482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.811056,0.115093,-0.000190343,1.62122e-07,-5.44705e-11,-67262.4,33.7731], Tmin=(100,'K'), Tmax=(784.738,'K')), NASAPolynomial(coeffs=[14.2026,0.0314882,-1.70102e-05,3.37871e-09,-2.37607e-13,-69400.9,-33.6315], Tmin=(784.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C]C(F)OO[C](F)CF(6352)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {4,S} {10,S} {13,S}
9  C u1 p0 c0 {3,S} {5,S} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-564.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,551,1088,1226,1380,1420,1481,3057,3119,355,410,600,1181,1341,1420,3056,395,473,707,1436,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.71679,'amu*angstrom^2'), symmetry=1, barrier=(39.4723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71529,'amu*angstrom^2'), symmetry=1, barrier=(39.438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830144,'amu*angstrom^2'), symmetry=1, barrier=(19.0866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.829818,'amu*angstrom^2'), symmetry=1, barrier=(19.0791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830797,'amu*angstrom^2'), symmetry=1, barrier=(19.1017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.743092,0.11585,-0.000200844,1.78906e-07,-6.19432e-11,-67758.6,33.9776], Tmin=(100,'K'), Tmax=(819.95,'K')), NASAPolynomial(coeffs=[12.2022,0.0344802,-1.86616e-05,3.68408e-09,-2.57043e-13,-69269.1,-22.1719], Tmin=(819.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(CsCsF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]CC(F)OOC(F)=[C]F(6953)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-322.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,2750,2850,1437.5,1250,1305,750,350,293,496,537,1218,167,640,1190,180,180,1812.88,1812.98],'cm^-1')),
        HinderedRotor(inertia=(0.162968,'amu*angstrom^2'), symmetry=1, barrier=(3.74695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16314,'amu*angstrom^2'), symmetry=1, barrier=(3.7509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16342,'amu*angstrom^2'), symmetry=1, barrier=(49.7412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16429,'amu*angstrom^2'), symmetry=1, barrier=(49.7614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.137096,0.108251,-0.000203333,2.02555e-07,-7.66115e-11,-38611.4,35.1407], Tmin=(100,'K'), Tmax=(827.235,'K')), NASAPolynomial(coeffs=[4.37835,0.0487202,-2.70329e-05,5.38823e-09,-3.77663e-13,-38068.6,22.0087], Tmin=(827.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-322.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCOJ) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O[CH]C(F)OOC(F)=[C]F(6954)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,293,496,537,1218,167,640,1190,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.824147,0.122797,-0.000232488,2.18452e-07,-7.74957e-11,-43982.2,36.5228], Tmin=(100,'K'), Tmax=(851.654,'K')), NASAPolynomial(coeffs=[9.73782,0.0378949,-2.07863e-05,4.07178e-09,-2.80539e-13,-44501.2,-5.22175], Tmin=(851.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-F1sO2sH)(O2s-H)(H)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'OC=C(F)OOC(F)=CF(6359)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u0 p0 c0 {6,S} {7,D} {11,S}
10 C u0 p0 c0 {3,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-656.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,262,390,483,597,572,732,631,807,1275,1439,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.320967,'amu*angstrom^2'), symmetry=1, barrier=(7.37965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62604,'amu*angstrom^2'), symmetry=1, barrier=(14.3939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.320392,'amu*angstrom^2'), symmetry=1, barrier=(7.36643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0842,'amu*angstrom^2'), symmetry=1, barrier=(47.9199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.225913,0.101484,-0.000161319,1.37894e-07,-4.73644e-11,-78770.8,33.2539], Tmin=(100,'K'), Tmax=(748.668,'K')), NASAPolynomial(coeffs=[11.7116,0.0338819,-1.82174e-05,3.64774e-09,-2.58942e-13,-80451.2,-20.167], Tmin=(748.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
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
    label = 'O=C[CH]OOC(F)=CF(6955)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {3,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u0 p0 c0 {5,D} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {7,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-359.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,326,540,652,719,1357,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,840.385],'cm^-1')),
        HinderedRotor(inertia=(0.202055,'amu*angstrom^2'), symmetry=1, barrier=(4.64564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202192,'amu*angstrom^2'), symmetry=1, barrier=(4.64878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25046,'amu*angstrom^2'), symmetry=1, barrier=(51.7425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.24657,'amu*angstrom^2'), symmetry=1, barrier=(51.6532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440431,0.0874638,-0.000141957,1.29898e-07,-4.75535e-11,-43173,30.2663], Tmin=(100,'K'), Tmax=(781.668,'K')), NASAPolynomial(coeffs=[7.63032,0.0372927,-2.0007e-05,3.99287e-09,-2.82484e-13,-43888.3,-0.0352919], Tmin=(781.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-359.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(OCJC=O)"""),
)

species(
    label = 'O=CC(F)OO[C]=CF(6956)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-301.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,247.287],'cm^-1')),
        HinderedRotor(inertia=(1.78202,'amu*angstrom^2'), symmetry=1, barrier=(40.9721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521897,'amu*angstrom^2'), symmetry=1, barrier=(11.9994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78204,'amu*angstrom^2'), symmetry=1, barrier=(40.9727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7821,'amu*angstrom^2'), symmetry=1, barrier=(40.9739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588261,0.0827448,-0.000115112,8.91164e-08,-2.85308e-11,-36090.3,28.8212], Tmin=(100,'K'), Tmax=(755.233,'K')), NASAPolynomial(coeffs=[9.78467,0.03404,-1.8383e-05,3.73599e-09,-2.69499e-13,-37479.4,-12.965], Tmin=(755.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-301.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-F1sH)(O2s-O2s))"""),
)

species(
    label = '[CH]=C(F)OOC(F)C=O(6957)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-297.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,293,496,537,1218,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.749887,'amu*angstrom^2'), symmetry=1, barrier=(17.2414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156415,'amu*angstrom^2'), symmetry=1, barrier=(3.59628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02319,'amu*angstrom^2'), symmetry=1, barrier=(46.5171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03013,'amu*angstrom^2'), symmetry=1, barrier=(46.6766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79649,0.0812414,-0.000103012,5.26112e-08,2.3365e-12,-35708,28.7853], Tmin=(100,'K'), Tmax=(543.095,'K')), NASAPolynomial(coeffs=[9.29223,0.0347302,-1.89109e-05,3.82916e-09,-2.74552e-13,-36867.7,-9.19608], Tmin=(543.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-297.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'O=C(F)[CH]F(215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-442.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(1.5365,'amu*angstrom^2'), symmetry=1, barrier=(35.3272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91969,0.0052474,7.61145e-05,-1.81228e-07,1.26791e-10,-53179.7,10.2212], Tmin=(10,'K'), Tmax=(489.215,'K')), NASAPolynomial(coeffs=[4.0714,0.0192239,-1.3397e-05,4.3331e-09,-5.27008e-13,-53376.6,7.73666], Tmin=(489.215,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-442.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)C=O(506)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-324.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.586986,'amu*angstrom^2'), symmetry=1, barrier=(13.496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82948,0.0149886,6.15997e-05,-2.22156e-07,2.15915e-10,-39081.1,10.2785], Tmin=(10,'K'), Tmax=(359.038,'K')), NASAPolynomial(coeffs=[4.10807,0.0226824,-1.56543e-05,5.05177e-09,-6.15171e-13,-39170.7,8.25062], Tmin=(359.038,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-324.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]C(F)CDO""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)=CF(214)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-249.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.358734,'amu*angstrom^2'), symmetry=1, barrier=(8.24799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3441.78,'J/mol'), sigma=(5.35945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.60 K, Pc=50.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75962,0.0203814,3.21936e-05,-1.12324e-07,8.28e-11,-30007.8,12.0802], Tmin=(10,'K'), Tmax=(521.492,'K')), NASAPolynomial(coeffs=[5.63612,0.021414,-1.5147e-05,4.91851e-09,-5.9756e-13,-30413.3,2.2378], Tmin=(521.492,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-249.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C[CH]F(880)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u0 p0 c0 {2,D} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-191.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.21522,'amu*angstrom^2'), symmetry=1, barrier=(35.2035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96722,0.00209998,6.04576e-05,-1.2409e-07,8.01311e-11,-23018.4,8.77531], Tmin=(10,'K'), Tmax=(478.869,'K')), NASAPolynomial(coeffs=[2.63637,0.0192853,-1.23829e-05,3.78104e-09,-4.41625e-13,-22960.5,13.4894], Tmin=(478.869,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-191.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHFCF[Z](72)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-58.4765,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([18.6764,93.2568,111.933],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([203.406,408.59,750.559,801.84,1025.92,1150.34,1361.98,1782.01,3238.53],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9644,0.00217073,4.76141e-05,-9.63411e-08,5.92621e-11,-7031.25,9.1333], Tmin=(10,'K'), Tmax=(529.917,'K')), NASAPolynomial(coeffs=[3.26037,0.0150322,-1.01557e-05,3.21347e-09,-3.84694e-13,-7062.61,11.0829], Tmin=(529.917,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-58.4765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""F[C]DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)C=O(2473)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-328.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.18317,'amu*angstrom^2'), symmetry=1, barrier=(27.2033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603055,'amu*angstrom^2'), symmetry=1, barrier=(13.8654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15439,0.0440801,-5.80197e-05,4.31461e-08,-1.32601e-11,-39422,18.3048], Tmin=(100,'K'), Tmax=(786.755,'K')), NASAPolynomial(coeffs=[7.18333,0.0185099,-9.26418e-06,1.82889e-09,-1.29953e-13,-40213.3,-4.7504], Tmin=(786.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[CH]OOC(F)=CF(3245)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u1 p0 c0 {3,S} {5,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-409.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,194,682,905,1196,1383,3221,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.226767,'amu*angstrom^2'), symmetry=1, barrier=(5.21382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22726,'amu*angstrom^2'), symmetry=1, barrier=(5.22516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46124,'amu*angstrom^2'), symmetry=1, barrier=(56.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91761,0.0785622,-0.000142421,1.3701e-07,-5.06695e-11,-49160.5,26.5988], Tmin=(100,'K'), Tmax=(821.716,'K')), NASAPolynomial(coeffs=[6.02591,0.0319352,-1.75834e-05,3.5008e-09,-2.45504e-13,-49265.4,7.42757], Tmin=(821.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-409.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sHO2s)"""),
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
    label = 'O=C[C](F)OOC(F)=CF(6349)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {1,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {5,S} {10,D}
9  C u0 p0 c0 {6,D} {7,S} {11,S}
10 C u0 p0 c0 {3,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-587.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,280,501,1494,1531,326,540,652,719,1357,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,1714.56],'cm^-1')),
        HinderedRotor(inertia=(2.26774,'amu*angstrom^2'), symmetry=1, barrier=(52.1398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370232,'amu*angstrom^2'), symmetry=1, barrier=(8.51236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26706,'amu*angstrom^2'), symmetry=1, barrier=(52.1242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26571,'amu*angstrom^2'), symmetry=1, barrier=(52.093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303638,0.0927505,-0.0001585,1.50582e-07,-5.61706e-11,-70581.4,31.3389], Tmin=(100,'K'), Tmax=(800.747,'K')), NASAPolynomial(coeffs=[6.43111,0.0411159,-2.23887e-05,4.47031e-09,-3.15398e-13,-70888.6,7.3481], Tmin=(800.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OOC(F)=CF(6958)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {12,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-568.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,326,540,652,719,1357,194,682,905,1196,1383,3221,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.318285,'amu*angstrom^2'), symmetry=1, barrier=(7.318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319027,'amu*angstrom^2'), symmetry=1, barrier=(7.33505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695043,'amu*angstrom^2'), symmetry=1, barrier=(15.9804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.35515,'amu*angstrom^2'), symmetry=1, barrier=(54.1496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.291489,0.10781,-0.000196666,1.83056e-07,-6.5209e-11,-68238.6,32.9364], Tmin=(100,'K'), Tmax=(835.335,'K')), NASAPolynomial(coeffs=[9.31136,0.0354688,-1.94333e-05,3.8353e-09,-2.66664e-13,-68923.3,-6.15952], Tmin=(835.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-568.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)OOC(F)=[C]F(6959)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-457.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,293,496,537,1218,167,640,1190,180,1110],'cm^-1')),
        HinderedRotor(inertia=(0.243117,'amu*angstrom^2'), symmetry=1, barrier=(5.58975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240557,'amu*angstrom^2'), symmetry=1, barrier=(5.53087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09488,'amu*angstrom^2'), symmetry=1, barrier=(48.1655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08944,'amu*angstrom^2'), symmetry=1, barrier=(48.0402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0277753,0.0989513,-0.00017188,1.61189e-07,-5.93226e-11,-54897.5,33.1484], Tmin=(100,'K'), Tmax=(795.472,'K')), NASAPolynomial(coeffs=[8.07791,0.0388867,-2.16873e-05,4.365e-09,-3.09134e-13,-55558.6,0.0480106], Tmin=(795.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-457.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O=CC(F)O[O+]=[C-]C(F)F(6960)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p1 c+1 {4,S} {10,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p1 c-1 {5,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-604.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.193985,0.107169,-0.000192876,1.85788e-07,-6.87152e-11,-72566.1,43.605], Tmin=(100,'K'), Tmax=(826.709,'K')), NASAPolynomial(coeffs=[6.22346,0.0452478,-2.45117e-05,4.84719e-09,-3.38588e-13,-72572.2,20.2462], Tmin=(826.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(F)(F)OOC(F)C=O(6961)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p1 c0 {8,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-527.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368821,0.105646,-0.000165607,1.39596e-07,-4.7679e-11,-63264.5,31.2059], Tmin=(100,'K'), Tmax=(714.641,'K')), NASAPolynomial(coeffs=[11.9827,0.0365034,-2.04626e-05,4.17818e-09,-3.00549e-13,-65029.6,-24.2317], Tmin=(714.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=CC(F)OOC(F)[C]F(6962)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {10,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p1 c0 {3,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-565.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.43585,'amu*angstrom^2'), symmetry=1, barrier=(33.0129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4361,'amu*angstrom^2'), symmetry=1, barrier=(33.0188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43605,'amu*angstrom^2'), symmetry=1, barrier=(33.0175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43601,'amu*angstrom^2'), symmetry=1, barrier=(33.0167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43582,'amu*angstrom^2'), symmetry=1, barrier=(33.0123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227569,0.091266,-0.00011732,8.00766e-08,-2.24221e-11,-67943.5,30.9861], Tmin=(100,'K'), Tmax=(860.19,'K')), NASAPolynomial(coeffs=[12.3474,0.0349065,-1.90389e-05,3.90605e-09,-2.84159e-13,-70028.6,-25.6597], Tmin=(860.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-565.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + group(Cds-OdCsH) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C=COOC(F)=CF(6938)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-316.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,194,682,905,1196,1383,3221,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.330755,'amu*angstrom^2'), symmetry=1, barrier=(7.60471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331792,'amu*angstrom^2'), symmetry=1, barrier=(7.62855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74675,'amu*angstrom^2'), symmetry=1, barrier=(40.1613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504755,0.0841757,-0.000136517,1.18333e-07,-4.07786e-11,-37917,29.4175], Tmin=(100,'K'), Tmax=(784.691,'K')), NASAPolynomial(coeffs=[10.1235,0.0277979,-1.47039e-05,2.91141e-09,-2.04817e-13,-39200.4,-13.214], Tmin=(784.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d)"""),
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
    label = 'C#COOC(F)C=O(6963)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {4,D} {5,S} {10,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-143.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.56322,'amu*angstrom^2'), symmetry=1, barrier=(35.9415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5649,'amu*angstrom^2'), symmetry=1, barrier=(35.9801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56521,'amu*angstrom^2'), symmetry=1, barrier=(35.9872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56327,'amu*angstrom^2'), symmetry=1, barrier=(35.9427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.993907,0.071731,-8.76039e-05,5.60936e-08,-1.46633e-11,-17154.8,24.0981], Tmin=(100,'K'), Tmax=(920.446,'K')), NASAPolynomial(coeffs=[11.3691,0.0266429,-1.4126e-05,2.87438e-09,-2.08505e-13,-19064.8,-25.0963], Tmin=(920.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'O=CC(F)OOC#CF(6964)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {2,S} {8,T}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-243.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,2175,525,239,401,1367,180],'cm^-1')),
        HinderedRotor(inertia=(1.44727,'amu*angstrom^2'), symmetry=1, barrier=(33.2757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45614,'amu*angstrom^2'), symmetry=1, barrier=(33.4794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447995,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.569143,0.0832942,-0.000126903,1.06152e-07,-3.6358e-11,-29141.1,27.934], Tmin=(100,'K'), Tmax=(709.959,'K')), NASAPolynomial(coeffs=[9.74685,0.0315959,-1.76964e-05,3.62403e-09,-2.61586e-13,-30444.5,-13.2011], Tmin=(709.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Ct-CtOs) + group(CtCF)"""),
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
    label = 'O=CC(F)O[O+]=[C-]F(2476)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p1 c+1 {3,S} {8,D}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
7  C u0 p0 c0 {5,D} {6,S} {10,S}
8  C u0 p1 c-1 {2,S} {4,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-618.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42871,0.0625753,-8.53649e-05,6.71604e-08,-2.21947e-11,-74259.8,23.6263], Tmin=(100,'K'), Tmax=(728.711,'K')), NASAPolynomial(coeffs=[7.62393,0.0285692,-1.53665e-05,3.12248e-09,-2.25351e-13,-75162.7,-4.30145], Tmin=(728.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-618.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(CJ2_singlet-FO)"""),
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
    E0 = (-383.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-96.9598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-18.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (132.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (10.6598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-322.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-226.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-52.5887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-218.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-222.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-213.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-217.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (25.0768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-19.7905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-158.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (35.1752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (94.0932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (97.3627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-383.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-118.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-64.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-54.8651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-43.7527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-32.2647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (76.5014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-119.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-72.9941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-157.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-70.8905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (169.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (69.4767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-145.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CC(F)OOC(F)=CF(6221)'],
    products = ['O=CC(=O)F(335)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(21.6381,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'FC=C(F)OOCF(3236)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.63854e-16,'m^3/(mol*s)'), n=6.80628, Ea=(308.187,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C[C]F-2(1596)', 'OOC(F)=CF(212)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(44.2643,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FC=C(F)OOC=COF(6948)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(149.789,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['FC=COOOC(F)=CF(6949)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(106.375,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=CC(F)OOC(F)=CF(6221)'],
    products = ['O=CC(F)OC(F)C(=O)F(6950)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(82.5351,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]1OOC(F)[CH]OC1F(6951)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C[C](F)OOC(F)=CF(6952)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C[C](F)OOC(F)[CH]F(6351)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OO[C](F)CF(6215)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]C(F)OOC(F)[CH]F(6915)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C]C(F)OO[C](F)CF(6352)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC(F)OOC(F)=[C]F(6953)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O[CH]C(F)OOC(F)=[C]F(6954)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OC=C(F)OOC(F)=CF(6359)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(75.1006,'s^-1'), n=3.11041, Ea=(175.051,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=C[CH]OOC(F)=CF(6955)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'O=CC(F)OO[C]=CF(6956)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', '[CH]=C(F)OOC(F)C=O(6957)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.11549e+09,'m^3/(mol*s)'), n=-0.68237, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.29431919638206655, var=1.0853977775937997, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C(F)[CH]F(215)', '[O]C(F)C=O(506)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(60.888,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC(F)=CF(214)', 'O=C[CH]F(880)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.53331,'m^3/(mol*s)'), n=1.7252, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_Sp-4R!H-2C_N-4R!H-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_Sp-4R!H-2C_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CHFCF[Z](72)', '[O]OC(F)C=O(2473)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCO(15)', 'F[CH]OOC(F)=CF(3245)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'O=C[C](F)OOC(F)=CF(6349)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(10.0672,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', 'O=[C]C(F)OOC(F)=CF(6958)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.85197e+07,'m^3/(mol*s)'), n=0.213787, Ea=(2.23984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00012681105802580086, var=0.012046587691161825, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_Ext-8R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_Ext-8R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(5)', 'O=CC(F)OOC(F)=[C]F(6959)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_N-Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_N-Sp-4C-2C_Ext-4C-R
Ea raised from -14.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=CC(F)O[O+]=[C-]C(F)F(6960)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.82066e+10,'s^-1'), n=0.827, Ea=(162.407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(F)(F)OOC(F)C=O(6961)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.82066e+10,'s^-1'), n=0.827, Ea=(132.014,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=CC(F)OOC(F)[C]F(6962)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.33334e+12,'s^-1'), n=-1.40567e-07, Ea=(86.0109,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'O=C=COOC(F)=CF(6938)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(204.235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F2(78)', 'C#COOC(F)C=O(6963)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction31',
    reactants = ['HF(38)', 'O=CC(F)OOC#CF(6964)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(271.613,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CHF(40)', 'O=CC(F)O[O+]=[C-]F(2476)'],
    products = ['O=CC(F)OOC(F)=CF(6221)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.1e+24,'cm^3/(mol*s)'), n=-3.8, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for CF
Exact match found for rate rule [CF]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #1843',
    isomers = [
        'O=CC(F)OOC(F)=CF(6221)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1843',
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

