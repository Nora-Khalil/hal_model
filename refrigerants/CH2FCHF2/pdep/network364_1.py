species(
    label = 'O=C=C(F)OOC(F)C=O(1919)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {5,D} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-590.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,197,221,431,657,2120,512.5,787.5,180,1753.44],'cm^-1')),
        HinderedRotor(inertia=(1.65169,'amu*angstrom^2'), symmetry=1, barrier=(37.9755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.017404,'amu*angstrom^2'), symmetry=1, barrier=(37.9665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65202,'amu*angstrom^2'), symmetry=1, barrier=(37.9833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458384,'amu*angstrom^2'), symmetry=1, barrier=(10.5392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903327,0.0811034,-8.91056e-05,6.76415e-09,4.38341e-11,-70926.5,29.8584], Tmin=(100,'K'), Tmax=(496.696,'K')), NASAPolynomial(coeffs=[9.3023,0.0357803,-1.96243e-05,3.96022e-09,-2.82391e-13,-72036.1,-7.55541], Tmin=(496.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-590.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
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
    label = 'O=C=C(F)OOCF(1903)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {5,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-472.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,548,1085,1183,1302,1466,1520,3060,3119,197,221,431,657,2120,512.5,787.5,180,1700.89],'cm^-1')),
        HinderedRotor(inertia=(0.362037,'amu*angstrom^2'), symmetry=1, barrier=(8.32394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0200088,'amu*angstrom^2'), symmetry=1, barrier=(41.1028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78853,'amu*angstrom^2'), symmetry=1, barrier=(41.1218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3027.97,'J/mol'), sigma=(5.38123,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=472.96 K, Pc=44.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23448,0.0682986,-0.000111342,1.03315e-07,-3.83624e-11,-56686.9,25.952], Tmin=(100,'K'), Tmax=(779.331,'K')), NASAPolynomial(coeffs=[6.38293,0.0304322,-1.64378e-05,3.29239e-09,-2.33459e-13,-57141.9,4.62646], Tmin=(779.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsFHHO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
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
    label = 'O=C=C(F)OO(2244)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-257.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,197,221,431,657,2120,512.5,787.5,361.244],'cm^-1')),
        HinderedRotor(inertia=(0.00128875,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00128865,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11824,0.0469859,-8.14597e-05,7.59519e-08,-2.75708e-11,-30928.4,18.9221], Tmin=(100,'K'), Tmax=(815.358,'K')), NASAPolynomial(coeffs=[5.85129,0.0183802,-9.90015e-06,1.96099e-09,-1.37408e-13,-31195.1,3.77219], Tmin=(815.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C=C(F)OOC=COF(2465)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {5,S} {7,D} {12,S}
9  C u0 p0 c0 {1,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-203.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,2995,3025,975,1000,1300,1375,400,500,1630,1680,197,221,431,657,2120,512.5,787.5,411.223,450.248,2340.71],'cm^-1')),
        HinderedRotor(inertia=(1.62048,'amu*angstrom^2'), symmetry=1, barrier=(37.258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258981,'amu*angstrom^2'), symmetry=1, barrier=(37.2577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647826,'amu*angstrom^2'), symmetry=1, barrier=(14.8948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103523,'amu*angstrom^2'), symmetry=1, barrier=(14.8947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0701427,0.0920674,-0.000128716,9.10115e-08,-2.54779e-11,-24322.5,33.3854], Tmin=(100,'K'), Tmax=(874.405,'K')), NASAPolynomial(coeffs=[14.8583,0.0244187,-1.26679e-05,2.53374e-09,-1.81354e-13,-26908.7,-35.9741], Tmin=(874.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C=C(F)OOOC=CF(2466)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {3,S} {4,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {1,S} {7,D} {12,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-281.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,197,221,431,657,2120,512.5,787.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.740969,0.0783195,-9.96177e-05,6.77723e-08,-1.88995e-11,-33732.8,34.2337], Tmin=(100,'K'), Tmax=(864.37,'K')), NASAPolynomial(coeffs=[11.1604,0.0301024,-1.59438e-05,3.23724e-09,-2.34256e-13,-35534,-14.5154], Tmin=(864.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=CC(F)OC(=O)C(=O)F(2380)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1060.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970076,0.0740418,-8.47003e-05,5.21553e-08,-1.3471e-11,-127438,28.2163], Tmin=(100,'K'), Tmax=(917.634,'K')), NASAPolynomial(coeffs=[10.3298,0.0332418,-1.80064e-05,3.70129e-09,-2.70073e-13,-129156,-16.1345], Tmin=(917.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1060.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-OdCsH) + group(COCFO)"""),
)

species(
    label = 'O=C1O[CH]C(F)OO[C]1F(2467)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {8,S} {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {3,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {6,D} {10,S}
10 C u1 p0 c0 {2,S} {5,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-579.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,2950,1000,280,501,1494,1531,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798188,0.022302,0.000148131,-2.47603e-07,1.07345e-10,-69576,26.9514], Tmin=(100,'K'), Tmax=(917.316,'K')), NASAPolynomial(coeffs=[39.6772,-0.0221392,1.62497e-05,-3.09728e-09,1.93597e-13,-81971.9,-185.949], Tmin=(917.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cycloheptane) + radical(CCsJOC(O)) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C[C](F)OOC(F)=C=O(2468)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {3,S} {7,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-260.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,197,221,431,657,2120,512.5,787.5,180,1851.17,1851.36,1852.25],'cm^-1')),
        HinderedRotor(inertia=(0.301119,'amu*angstrom^2'), symmetry=1, barrier=(6.92332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299917,'amu*angstrom^2'), symmetry=1, barrier=(6.89569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10584,'amu*angstrom^2'), symmetry=1, barrier=(48.4174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10924,'amu*angstrom^2'), symmetry=1, barrier=(48.4955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.108422,0.110322,-0.000220984,2.22505e-07,-8.32545e-11,-31211,35.533], Tmin=(100,'K'), Tmax=(848.349,'K')), NASAPolynomial(coeffs=[3.83564,0.0457572,-2.55456e-05,5.0497e-09,-3.50141e-13,-30226,26.9031], Tmin=(848.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OO[C](F)C=O(2256)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
8  C u1 p0 c0 {2,S} {4,S} {9,S}
9  C u0 p0 c0 {5,D} {8,S} {12,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-515.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,280,501,1494,1531,2782.5,750,1395,475,1775,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.91507,'amu*angstrom^2'), symmetry=1, barrier=(44.0312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92736,'amu*angstrom^2'), symmetry=1, barrier=(44.3138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93108,'amu*angstrom^2'), symmetry=1, barrier=(44.3993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.557682,'amu*angstrom^2'), symmetry=1, barrier=(12.8222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91897,'amu*angstrom^2'), symmetry=1, barrier=(44.1209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.129346,0.101151,-0.000171496,1.54202e-07,-5.4364e-11,-61849.5,32.0431], Tmin=(100,'K'), Tmax=(809.546,'K')), NASAPolynomial(coeffs=[10.0019,0.0349021,-1.87447e-05,3.70567e-09,-2.59487e-13,-62959.3,-11.4172], Tmin=(809.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=C[C](F)OO[C](F)C=O(1914)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u1 p0 c0 {2,S} {4,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-534.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,253,307,474,528,1484,1504,1502,1560,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,277.548],'cm^-1')),
        HinderedRotor(inertia=(0.941161,'amu*angstrom^2'), symmetry=1, barrier=(51.3828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940658,'amu*angstrom^2'), symmetry=1, barrier=(51.3846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94094,'amu*angstrom^2'), symmetry=1, barrier=(51.3849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939978,'amu*angstrom^2'), symmetry=1, barrier=(51.3831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941187,'amu*angstrom^2'), symmetry=1, barrier=(51.3827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3968.58,'J/mol'), sigma=(6.11625,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=619.88 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493724,0.0857462,-0.000132079,1.20146e-07,-4.47372e-11,-64193.4,30.3466], Tmin=(100,'K'), Tmax=(755.742,'K')), NASAPolynomial(coeffs=[6.98817,0.0407898,-2.18447e-05,4.37587e-09,-3.11206e-13,-64872.8,2.83287], Tmin=(755.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-534.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OOC(F)[C]=O(2469)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,S} {12,S}
9  C u1 p0 c0 {5,D} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-496.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,308,402,315,505,557,643,1127,1235,1327,1355,1420,1420,3041,3071,1850,1860,440,470,900,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.709264,'amu*angstrom^2'), symmetry=1, barrier=(16.3074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708951,'amu*angstrom^2'), symmetry=1, barrier=(16.3002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64665,'amu*angstrom^2'), symmetry=1, barrier=(37.8598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64462,'amu*angstrom^2'), symmetry=1, barrier=(37.8131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709212,'amu*angstrom^2'), symmetry=1, barrier=(16.3062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729691,0.116289,-0.000210035,1.8736e-07,-6.38184e-11,-59506.5,33.6584], Tmin=(100,'K'), Tmax=(849.62,'K')), NASAPolynomial(coeffs=[12.8601,0.0292934,-1.58117e-05,3.07602e-09,-2.11202e-13,-60985,-24.8011], Tmin=(849.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-496.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(COj(Cs-F1sO2sH)(O2d)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=C[C](F)OOC(F)=[C]O(2265)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {10,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u0 p0 c0 {6,D} {7,S} {11,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-363.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,280,501,1494,1531,293,496,537,1218,2782.5,750,1395,475,1775,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.135642,'amu*angstrom^2'), symmetry=1, barrier=(3.11867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247436,'amu*angstrom^2'), symmetry=1, barrier=(3.07208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16987,'amu*angstrom^2'), symmetry=1, barrier=(49.8897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17053,'amu*angstrom^2'), symmetry=1, barrier=(49.9047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17024,'amu*angstrom^2'), symmetry=1, barrier=(49.898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0474438,0.0978295,-0.000167751,1.54357e-07,-5.56359e-11,-43556.5,35.084], Tmin=(100,'K'), Tmax=(806.831,'K')), NASAPolynomial(coeffs=[8.66021,0.03673,-1.99511e-05,3.96697e-09,-2.78795e-13,-44347.4,-0.907435], Tmin=(806.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CsCOF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]C(F)OOC(F)=[C]O(2470)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {10,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u1 p0 c0 {6,D} {7,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-343.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,355,410,600,1181,1341,1420,3056,293,496,537,1218,1855,455,950,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.320997,'amu*angstrom^2'), symmetry=1, barrier=(7.38035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3211,'amu*angstrom^2'), symmetry=1, barrier=(7.38273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321002,'amu*angstrom^2'), symmetry=1, barrier=(7.38048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7041,'amu*angstrom^2'), symmetry=1, barrier=(39.1807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70434,'amu*angstrom^2'), symmetry=1, barrier=(39.1861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.549852,0.112924,-0.000206096,1.8718e-07,-6.48977e-11,-41213.7,36.6888], Tmin=(100,'K'), Tmax=(844.195,'K')), NASAPolynomial(coeffs=[11.5239,0.0311118,-1.70127e-05,3.33603e-09,-2.30403e-13,-42375.5,-14.3227], Tmin=(844.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(COj(Cs-F1sO2sH)(O2d)) + radical(C=CJO)"""),
)

species(
    label = 'O=C=C(F)OOC(F)=CO(2260)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {11,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-519.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,326,540,652,719,1357,3010,987.5,1337.5,450,1655,197,221,431,657,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.243203,'amu*angstrom^2'), symmetry=1, barrier=(5.59172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243316,'amu*angstrom^2'), symmetry=1, barrier=(5.59432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235996,'amu*angstrom^2'), symmetry=1, barrier=(33.1391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215315,0.0932282,-0.000137186,9.75478e-08,-2.32033e-11,-62364.5,32.444], Tmin=(100,'K'), Tmax=(605.951,'K')), NASAPolynomial(coeffs=[11.8916,0.0303103,-1.64871e-05,3.31742e-09,-2.3633e-13,-64039.5,-20.1829], Tmin=(605.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-519.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C=C(F)OO[CH]C=O(2471)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {2,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {5,D} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-223.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,197,221,431,657,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,180,742.644],'cm^-1')),
        HinderedRotor(inertia=(0.00462131,'amu*angstrom^2'), symmetry=1, barrier=(1.81076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0785284,'amu*angstrom^2'), symmetry=1, barrier=(1.80552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132291,'amu*angstrom^2'), symmetry=1, barrier=(51.8188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25371,'amu*angstrom^2'), symmetry=1, barrier=(51.8172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605758,0.0832748,-0.000137017,1.24767e-07,-4.52733e-11,-26754.9,30.4011], Tmin=(100,'K'), Tmax=(783.922,'K')), NASAPolynomial(coeffs=[7.96827,0.033412,-1.80806e-05,3.61296e-09,-2.5555e-13,-27531.4,-0.916793], Tmin=(783.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(OCJC=O)"""),
)

species(
    label = 'O=C=[C]OOC(F)C=O(2472)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u0 p0 c0 {4,D} {6,S} {11,S}
8  C u1 p0 c0 {3,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-162.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,1685,370,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.500181,'amu*angstrom^2'), symmetry=1, barrier=(11.5001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32603,'amu*angstrom^2'), symmetry=1, barrier=(7.49607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79492,'amu*angstrom^2'), symmetry=1, barrier=(41.2688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79126,'amu*angstrom^2'), symmetry=1, barrier=(41.1846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524408,0.0841295,-0.000137284,1.20677e-07,-4.2206e-11,-19436.5,31.9453], Tmin=(100,'K'), Tmax=(784.365,'K')), NASAPolynomial(coeffs=[9.58492,0.0292137,-1.56068e-05,3.10109e-09,-2.18578e-13,-20589.9,-7.85784], Tmin=(784.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]C(=O)F(1543)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {3,D} {4,S}
"""),
    E0 = (-325.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([233,496,705,1150,2014,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.06462,'amu*angstrom^2'), symmetry=1, barrier=(24.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90371,0.00648742,6.66758e-05,-1.78696e-07,1.33913e-10,-39126.2,10.0186], Tmin=(10,'K'), Tmax=(477.098,'K')), NASAPolynomial(coeffs=[4.98548,0.0142778,-1.08248e-05,3.66779e-09,-4.58866e-13,-39421.3,3.58924], Tmin=(477.098,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-325.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OD[C]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OC(F)=C=O(2248)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u1 p2 c0 {2,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (-105.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,197,221,431,657,2120,512.5,787.5,2787.33],'cm^-1')),
        HinderedRotor(inertia=(0.0346032,'amu*angstrom^2'), symmetry=1, barrier=(3.05444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37136,0.043506,-8.64301e-05,8.52867e-08,-3.11222e-11,-12657.6,18.7731], Tmin=(100,'K'), Tmax=(866.808,'K')), NASAPolynomial(coeffs=[4.28985,0.0166729,-8.88159e-06,1.71369e-09,-1.1673e-13,-12314.7,13.688], Tmin=(866.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
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
    label = 'O=C=[C]F(289)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (33.6712,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(58.9933,'amu')),
        NonlinearRotor(inertia=([3.58812,117.154,120.742],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([280.764,365.894,563.98,866.433,1429.66,2068.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94745,0.00353157,4.23087e-05,-1.10022e-07,8.19717e-11,4052.13,8.19715], Tmin=(10,'K'), Tmax=(472.411,'K')), NASAPolynomial(coeffs=[4.41104,0.00925722,-6.51505e-06,2.12227e-09,-2.60232e-13,3900.64,5.16846], Tmin=(472.411,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(33.6712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCD[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C=C(F)OO[CH]F(1938)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {8,D}
7 C u1 p0 c0 {2,S} {4,S} {9,S}
8 C u0 p0 c0 {5,D} {6,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-273.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,197,221,431,657,580,1155,1237,1373,3147,2120,512.5,787.5,180,1031.33],'cm^-1')),
        HinderedRotor(inertia=(0.269071,'amu*angstrom^2'), symmetry=1, barrier=(6.18646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00817956,'amu*angstrom^2'), symmetry=1, barrier=(6.17956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26103,'amu*angstrom^2'), symmetry=1, barrier=(51.9854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08558,0.0743408,-0.000137365,1.31735e-07,-4.83364e-11,-32742.5,26.7242], Tmin=(100,'K'), Tmax=(826.293,'K')), NASAPolynomial(coeffs=[6.35109,0.0280775,-1.56707e-05,3.12423e-09,-2.18852e-13,-32903.5,6.61711], Tmin=(826.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsFHHO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(CsF1sHO2s)"""),
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
    label = 'O=C=C(F)OO[C](F)C=O(2255)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-451.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,280,501,1494,1531,197,221,431,657,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.121297,'amu*angstrom^2'), symmetry=1, barrier=(2.78886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0438828,'amu*angstrom^2'), symmetry=1, barrier=(52.1522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26971,'amu*angstrom^2'), symmetry=1, barrier=(52.1852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0438685,'amu*angstrom^2'), symmetry=1, barrier=(52.1849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.470098,0.0885471,-0.000153507,1.45381e-07,-5.38613e-11,-54163.3,31.4697], Tmin=(100,'K'), Tmax=(803.658,'K')), NASAPolynomial(coeffs=[6.76445,0.0372435,-2.04674e-05,4.09163e-09,-2.88567e-13,-54529.9,6.49219], Tmin=(803.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OOC(F)=C=O(2474)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u1 p0 c0 {5,D} {7,S}
10 C u0 p0 c0 {6,D} {8,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-431.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,197,221,431,657,1855,455,950,2120,512.5,787.5,180,1765.1],'cm^-1')),
        HinderedRotor(inertia=(0.273646,'amu*angstrom^2'), symmetry=1, barrier=(6.29166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273454,'amu*angstrom^2'), symmetry=1, barrier=(6.28725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94169,'amu*angstrom^2'), symmetry=1, barrier=(44.6433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94218,'amu*angstrom^2'), symmetry=1, barrier=(44.6544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.123526,0.10359,-0.000191625,1.77821e-07,-6.29064e-11,-51820.6,33.0618], Tmin=(100,'K'), Tmax=(839.867,'K')), NASAPolynomial(coeffs=[9.63179,0.0316193,-1.75255e-05,3.45991e-09,-2.40111e-13,-52559.6,-6.9434], Tmin=(839.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)O[O+]=[C-]C(=O)F(2475)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p1 c+1 {3,S} {10,D}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {5,D} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {6,D} {10,S}
10 C u0 p1 c-1 {4,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-537.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0897775,0.104502,-0.000191594,1.83485e-07,-6.71091e-11,-64454.6,43.2514], Tmin=(100,'K'), Tmax=(833.327,'K')), NASAPolynomial(coeffs=[6.90976,0.0406451,-2.2183e-05,4.38129e-09,-3.05172e-13,-64570.6,17.0628], Tmin=(833.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COCFO) + group(CsJ2_singlet-CsH)"""),
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
    label = 'O=C=COOC(F)=C=O(2463)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {4,D} {6,D}
9  C u0 p0 c0 {5,D} {7,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-179.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,197,221,431,657,2110,2130,495,530,650,925,218.63,218.812],'cm^-1')),
        HinderedRotor(inertia=(0.131199,'amu*angstrom^2'), symmetry=1, barrier=(4.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01047,'amu*angstrom^2'), symmetry=1, barrier=(34.3609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00707723,'amu*angstrom^2'), symmetry=1, barrier=(34.3624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668943,0.0800021,-0.000131642,1.13306e-07,-3.85542e-11,-21498.8,29.5562], Tmin=(100,'K'), Tmax=(788.969,'K')), NASAPolynomial(coeffs=[10.4628,0.0239148,-1.2776e-05,2.53114e-09,-1.77851e-13,-22844,-14.1029], Tmin=(788.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + missing(Cdd-CdO2d)"""),
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
    E0 = (-286.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1.46437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (101.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (230.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (109.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-243.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-295.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (45.8355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-192.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-211.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-187.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-39.9124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-34.8487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-60.4246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (133.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (194.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-292.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-12.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-10.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (43.5591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (54.6714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (66.1595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-60.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (27.5337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-305.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=C(F)OOC(F)C=O(1919)'],
    products = ['O=CC(=O)F(335)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(20.3866,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'O=C=C(F)OOCF(1903)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.63854e-16,'m^3/(mol*s)'), n=6.80628, Ea=(308.187,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C[C]F-2(1596)', 'O=C=C(F)OO(2244)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(38.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=C(F)OOC=COF(2465)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(149.789,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C=C(F)OOOC=CF(2466)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(106.375,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C=C(F)OOC(F)C=O(1919)'],
    products = ['O=CC(F)OC(=O)C(=O)F(2380)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(63.3761,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C1O[CH]C(F)OO[C]1F(2467)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R7JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C[C](F)OOC(F)=C=O(2468)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]C(F)OO[C](F)C=O(2256)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OO[C](F)C=O(1914)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]C(F)OOC(F)[C]=O(2469)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C[C](F)OOC(F)=[C]O(2265)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]C(F)OOC(F)=[C]O(2470)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C=C(F)OOC(F)=CO(2260)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(75.1006,'s^-1'), n=3.11041, Ea=(175.051,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'O=C=C(F)OO[CH]C=O(2471)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=C=[C]OOC(F)C=O(2472)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]C(=O)F(1543)', '[O]C(F)C=O(506)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(73.1646,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC(F)=C=O(2248)', 'O=C[CH]F(880)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.53331,'m^3/(mol*s)'), n=1.7252, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_Sp-4R!H-2C_N-4R!H-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_Sp-4R!H-2C_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C=[C]F(289)', '[O]OC(F)C=O(2473)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.004135,'m^3/(mol*s)'), n=2.525, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_N-Sp-4R!H-2C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H-4R!H_N-Sp-4R!H-2C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HCO(15)', 'O=C=C(F)OO[CH]F(1938)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', 'O=C=C(F)OO[C](F)C=O(2255)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(10.0672,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'O=[C]C(F)OOC(F)=C=O(2474)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.85197e+07,'m^3/(mol*s)'), n=0.213787, Ea=(2.23984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00012681105802580086, var=0.012046587691161825, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_Ext-8R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_Ext-8R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=CC(F)O[O+]=[C-]C(=O)F(2475)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(192.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=C=COOC(F)=C=O(2463)'],
    products = ['O=C=C(F)OOC(F)C=O(1919)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(204.235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C=C(F)OOC(F)C=O(1919)'],
    products = ['CO(13)', 'O=CC(F)O[O+]=[C-]F(2476)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.15699e+14,'s^-1'), n=0.0573689, Ea=(0.639264,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='RFC=C=O',), comment="""Estimated from node RFC=C=O"""),
)

network(
    label = 'PDepNetwork #364',
    isomers = [
        'O=C=C(F)OOC(F)C=O(1919)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #364',
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

