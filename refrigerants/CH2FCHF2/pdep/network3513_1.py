species(
    label = 'O=C(F)OC(O)C(=O)F(13079)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-1165.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,482,664,788,1296,1923,264.58,264.759,265.014,1070.77],'cm^-1')),
        HinderedRotor(inertia=(0.00240082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00615519,'amu*angstrom^2'), symmetry=1, barrier=(9.28753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514388,'amu*angstrom^2'), symmetry=1, barrier=(25.574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514667,'amu*angstrom^2'), symmetry=1, barrier=(25.5735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29994,0.0631584,-6.81823e-05,3.7544e-08,-8.42671e-12,-140078,28.9111], Tmin=(100,'K'), Tmax=(1062.05,'K')), NASAPolynomial(coeffs=[11.787,0.023661,-1.23973e-05,2.5266e-09,-1.83799e-13,-142305,-22.314], Tmin=(1062.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1165.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + group(COFOO)"""),
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
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'O=C(F)OC(O)F(5365)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-1020.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,482,664,788,1296,1923,222.971,223.169,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.265891,'amu*angstrom^2'), symmetry=1, barrier=(9.39111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817569,'amu*angstrom^2'), symmetry=1, barrier=(28.8891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608996,'amu*angstrom^2'), symmetry=1, barrier=(21.5108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3648.6,'J/mol'), sigma=(5.56031,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.90 K, Pc=48.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71583,0.0500006,-5.30385e-05,2.7775e-08,-5.77561e-12,-122651,21.0141], Tmin=(100,'K'), Tmax=(1161.05,'K')), NASAPolynomial(coeffs=[12.0233,0.0144896,-7.16047e-06,1.43206e-09,-1.03373e-13,-125045,-30.2527], Tmin=(1161.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1020.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO)"""),
)

species(
    label = 'O=C(F)C(O)OF(12539)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {2,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-657.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,220.55,220.573],'cm^-1')),
        HinderedRotor(inertia=(0.00346516,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111806,'amu*angstrom^2'), symmetry=1, barrier=(3.8595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618985,'amu*angstrom^2'), symmetry=1, barrier=(21.3696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70374,0.0545929,-7.11829e-05,4.86337e-08,-1.34331e-11,-78975.3,23.2397], Tmin=(100,'K'), Tmax=(877.761,'K')), NASAPolynomial(coeffs=[9.5936,0.0186385,-9.74065e-06,1.96791e-09,-1.42004e-13,-80360.4,-13.7956], Tmin=(877.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-657.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsOsOsH) + group(COCsFO)"""),
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
    label = 'O=C=COC(=O)F(13334)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {7,D} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {3,D}
7 C u0 p0 c0 {4,D} {5,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-544.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,482,664,788,1296,1923,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04794,'amu*angstrom^2'), symmetry=1, barrier=(24.0941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04663,'amu*angstrom^2'), symmetry=1, barrier=(24.0642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32069,0.0544445,-6.55648e-05,3.66851e-08,-7.81407e-12,-65342,20.0955], Tmin=(100,'K'), Tmax=(1164.35,'K')), NASAPolynomial(coeffs=[15.6233,0.0053095,-2.26565e-06,4.42221e-10,-3.22907e-14,-68672.7,-51.0825], Tmin=(1164.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-544.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(COFOO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C(F)OOC(F)=CO(13335)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {10,S}
9  C u0 p0 c0 {2,S} {4,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-861.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,326,540,652,719,1357,3010,987.5,1337.5,450,1655,482,664,788,1296,1923,180],'cm^-1')),
        HinderedRotor(inertia=(1.442,'amu*angstrom^2'), symmetry=1, barrier=(33.1545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43874,'amu*angstrom^2'), symmetry=1, barrier=(33.0794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44431,'amu*angstrom^2'), symmetry=1, barrier=(33.2074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43998,'amu*angstrom^2'), symmetry=1, barrier=(33.1079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164618,0.0723302,-6.43916e-05,1.86564e-08,8.61713e-13,-103507,30.3864], Tmin=(100,'K'), Tmax=(1023.86,'K')), NASAPolynomial(coeffs=[22.3788,0.00823347,-3.72793e-06,8.00448e-10,-6.32455e-14,-109246,-83.1159], Tmin=(1023.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-861.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OC=C(F)OO(13336)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {7,D}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-759.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,482,664,788,1296,1923,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.145496,'amu*angstrom^2'), symmetry=1, barrier=(3.34523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34868,'amu*angstrom^2'), symmetry=1, barrier=(31.0088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538416,'amu*angstrom^2'), symmetry=1, barrier=(12.3793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34884,'amu*angstrom^2'), symmetry=1, barrier=(31.0125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502116,0.0845738,-0.000130416,1.06663e-07,-3.52276e-11,-91229.3,28.2989], Tmin=(100,'K'), Tmax=(738.869,'K')), NASAPolynomial(coeffs=[10.9349,0.0280982,-1.57715e-05,3.22911e-09,-2.32883e-13,-92771.1,-18.8768], Tmin=(738.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-759.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COFOO)"""),
)

species(
    label = 'OC1O[C](F)OO[C]1F(13337)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
8  C u1 p0 c0 {1,S} {4,S} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-542.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,1000,395,473,707,1436,482,586,761,1411,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.627074,0.0752139,-8.15951e-05,4.23454e-08,-7.92495e-12,-65079.2,30.8792], Tmin=(100,'K'), Tmax=(1614.03,'K')), NASAPolynomial(coeffs=[18.4817,0.0060669,2.91755e-06,-9.26598e-10,7.31244e-14,-68409.3,-61.6653], Tmin=(1614.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-542.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFHOO) + ring(124trioxane) + radical(CsCsF1sO2s) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]C(F)[C](O)OC(=O)F(13338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {8,S} {9,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-784.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,391,562,707,872,1109,1210,1289,3137,360,370,350,482,664,788,1296,1923,285.231,285.246,285.26,285.29],'cm^-1')),
        HinderedRotor(inertia=(0.520656,'amu*angstrom^2'), symmetry=1, barrier=(30.0568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00571204,'amu*angstrom^2'), symmetry=1, barrier=(30.0571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000804058,'amu*angstrom^2'), symmetry=1, barrier=(9.12932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520536,'amu*angstrom^2'), symmetry=1, barrier=(30.0565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867135,0.0741138,-9.78593e-05,6.64417e-08,-1.80944e-11,-94205.3,31.0998], Tmin=(100,'K'), Tmax=(892.626,'K')), NASAPolynomial(coeffs=[12.2399,0.0231495,-1.22152e-05,2.4761e-09,-1.79001e-13,-96235.6,-22.475], Tmin=(892.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-784.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(COFOO) + radical(O2sj(Cs-CsF1sH)) + radical(Cs_P)"""),
)

species(
    label = '[O]C(F)O[C](O)C(=O)F(13339)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u0 p0 c0 {2,S} {6,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-738.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,451,553,637,1069,1180,1265,1301,3056,360,370,350,611,648,830,1210,1753,222.731,224.147,224.894,226.351],'cm^-1')),
        HinderedRotor(inertia=(0.265735,'amu*angstrom^2'), symmetry=1, barrier=(9.41673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00335176,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00333952,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00896137,'amu*angstrom^2'), symmetry=1, barrier=(31.4131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.886545,0.0743009,-0.000102859,7.49745e-08,-2.2039e-11,-88662.7,31.9942], Tmin=(100,'K'), Tmax=(827.772,'K')), NASAPolynomial(coeffs=[11.1461,0.0247249,-1.30237e-05,2.62494e-09,-1.88643e-13,-90361.3,-15.5632], Tmin=(827.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFHOO) + group(COCsFO) + radical(O2sj(Cs-F1sO2sH)) + radical(Cs_P)"""),
)

species(
    label = '[O]C(F)C([O])OC(=O)F(13340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-763.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,482,664,788,1296,1923,180,180,180,1558.3,1559.85],'cm^-1')),
        HinderedRotor(inertia=(0.150251,'amu*angstrom^2'), symmetry=1, barrier=(3.45457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76339,'amu*angstrom^2'), symmetry=1, barrier=(40.5438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150812,'amu*angstrom^2'), symmetry=1, barrier=(3.46747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41985,0.0664201,-7.45649e-05,2.18485e-08,1.78672e-11,-91769.7,29.3686], Tmin=(100,'K'), Tmax=(516.135,'K')), NASAPolynomial(coeffs=[7.804,0.0318035,-1.7148e-05,3.46766e-09,-2.48705e-13,-92626.7,0.873706], Tmin=(516.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(COFOO) + radical(CCOJ) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O=C(F)[C](O)O[C](O)F(13341)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {3,S} {4,S} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {5,S}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-788.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,482,586,761,1411,611,648,830,1210,1753,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683157,0.0779629,-0.000106398,7.31018e-08,-1.99333e-11,-94752.3,33.0598], Tmin=(100,'K'), Tmax=(895.549,'K')), NASAPolynomial(coeffs=[13.3987,0.021169,-1.12722e-05,2.28889e-09,-1.65501e-13,-97029.8,-26.8828], Tmin=(895.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-788.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFHOO) + group(COCsFO) + radical(Cs_P) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = '[O]C(OC(=O)F)[C](O)F(13342)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
8  C u1 p0 c0 {1,S} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-795.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,395,473,707,1436,482,664,788,1296,1923,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541516,0.0865992,-0.00015339,1.44014e-07,-5.25046e-11,-95589.8,32.3192], Tmin=(100,'K'), Tmax=(811.414,'K')), NASAPolynomial(coeffs=[7.60667,0.0325633,-1.79913e-05,3.59635e-09,-2.53143e-13,-96104,3.60656], Tmin=(811.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-795.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(COFOO) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)OC([O])C(=O)F(13343)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {3,S} {5,S} {11,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-699.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,451,553,637,1069,1180,1265,1301,3056,486,617,768,1157,1926,222.715,222.741,222.767,222.817,2217.63],'cm^-1')),
        HinderedRotor(inertia=(0.00339379,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00805595,'amu*angstrom^2'), symmetry=1, barrier=(28.1133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339931,'amu*angstrom^2'), symmetry=1, barrier=(11.9757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00884,0.0716507,-9.99857e-05,7.61914e-08,-2.36708e-11,-84038.4,31.6446], Tmin=(100,'K'), Tmax=(781.948,'K')), NASAPolynomial(coeffs=[9.78491,0.026758,-1.38699e-05,2.77251e-09,-1.98097e-13,-85410.9,-8.53629], Tmin=(781.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFHOO) + group(COCsFO) + radical(C=OCOJ) + radical(O2sj(Cs-F1sO2sH))"""),
)

species(
    label = '[O]C(O[C](O)F)C(=O)F(13344)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {9,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-750.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,482,586,761,1411,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.838681,0.0748719,-0.000101736,7.15957e-08,-2.01979e-11,-90129.4,32.5941], Tmin=(100,'K'), Tmax=(863.333,'K')), NASAPolynomial(coeffs=[11.9605,0.0233413,-1.22023e-05,2.45688e-09,-1.7669e-13,-92049.7,-19.4278], Tmin=(863.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-750.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFHOO) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = 'O=C(F)OC(O)=C(O)F(13345)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-1041.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,350,440,435,1725,326,540,652,719,1357,482,664,788,1296,1923,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.76589,'amu*angstrom^2'), symmetry=1, barrier=(17.6093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765938,'amu*angstrom^2'), symmetry=1, barrier=(17.6104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765696,'amu*angstrom^2'), symmetry=1, barrier=(17.6049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76572,'amu*angstrom^2'), symmetry=1, barrier=(17.6054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0760927,0.0946972,-0.000143328,1.04376e-07,-2.94873e-11,-125078,26.1202], Tmin=(100,'K'), Tmax=(873.413,'K')), NASAPolynomial(coeffs=[17.0079,0.0164567,-8.95684e-06,1.81174e-09,-1.29943e-13,-128063,-53.9879], Tmin=(873.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1041.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COFOO)"""),
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
    label = 'O=[C]C(O)OC(=O)F(13346)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,D}
8  C u1 p0 c0 {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-756.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,482,664,788,1296,1923,1855,455,950,329.28,329.28,1301.16],'cm^-1')),
        HinderedRotor(inertia=(0.00155476,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0590371,'amu*angstrom^2'), symmetry=1, barrier=(4.54238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336138,'amu*angstrom^2'), symmetry=1, barrier=(25.8643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336203,'amu*angstrom^2'), symmetry=1, barrier=(25.8643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70788,0.0534983,-5.5292e-05,2.9079e-08,-6.25676e-12,-90940.8,27.1975], Tmin=(100,'K'), Tmax=(1101.59,'K')), NASAPolynomial(coeffs=[10.6836,0.0209061,-1.0912e-05,2.22068e-09,-1.61376e-13,-92918.3,-16.9737], Tmin=(1101.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(COFOO) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]OC(O)C(=O)F(13347)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7  C u0 p0 c0 {1,S} {4,D} {6,S}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-746.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,1855,455,950,180,301.169,560.139],'cm^-1')),
        HinderedRotor(inertia=(0.0440734,'amu*angstrom^2'), symmetry=1, barrier=(9.96842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00582153,'amu*angstrom^2'), symmetry=1, barrier=(37.2048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163462,'amu*angstrom^2'), symmetry=1, barrier=(37.2297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61869,'amu*angstrom^2'), symmetry=1, barrier=(37.2169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59084,0.0575877,-6.67702e-05,4.13732e-08,-1.05906e-11,-89688.6,25.8118], Tmin=(100,'K'), Tmax=(933.81,'K')), NASAPolynomial(coeffs=[9.50607,0.0236817,-1.23048e-05,2.48825e-09,-1.80011e-13,-91166.8,-11.8323], Tmin=(933.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(=O)F(341)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[CH]O(1552)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-466.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,634.702],'cm^-1')),
        HinderedRotor(inertia=(0.178809,'amu*angstrom^2'), symmetry=1, barrier=(51.1156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178808,'amu*angstrom^2'), symmetry=1, barrier=(51.1155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9486,0.00327696,8.62766e-05,-1.78856e-07,1.15095e-10,-56066.9,9.87278], Tmin=(10,'K'), Tmax=(486.336,'K')), NASAPolynomial(coeffs=[2.12805,0.0278262,-1.89743e-05,5.9068e-09,-6.93439e-13,-56003.1,16.1793], Tmin=(486.336,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-466.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)[CH]O""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CFO(51)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)C(=O)F(12721)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-548.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,180,3056.44],'cm^-1')),
        HinderedRotor(inertia=(0.287431,'amu*angstrom^2'), symmetry=1, barrier=(6.6086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257157,'amu*angstrom^2'), symmetry=1, barrier=(12.1649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26506,0.0412616,-5.57075e-05,4.24648e-08,-1.329e-11,-65858.1,21.8266], Tmin=(100,'K'), Tmax=(775.281,'K')), NASAPolynomial(coeffs=[7.02534,0.0167011,-8.18786e-06,1.60206e-09,-1.13152e-13,-66596.2,0.0727112], Tmin=(775.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(F)[CH]OC(=O)F(2174)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {3,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-756.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,482,664,788,1296,1923,180,180,180,993.435],'cm^-1')),
        HinderedRotor(inertia=(0.437388,'amu*angstrom^2'), symmetry=1, barrier=(19.3387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841242,'amu*angstrom^2'), symmetry=1, barrier=(19.3418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0276088,'amu*angstrom^2'), symmetry=1, barrier=(19.3403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27907,0.0603605,-7.18957e-05,4.10935e-08,-9.19399e-12,-90933.4,24.5936], Tmin=(100,'K'), Tmax=(1090.49,'K')), NASAPolynomial(coeffs=[14.1284,0.0132274,-7.06182e-06,1.45698e-09,-1.07003e-13,-93735.8,-38.5102], Tmin=(1090.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COFOO) + radical(CCsJOC(O))"""),
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
    label = '[O]C(OC(=O)F)C(=O)F(13348)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-921.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,482,664,788,1296,1923,360.234,360.234,360.253,360.279,360.294],'cm^-1')),
        HinderedRotor(inertia=(0.00129887,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00129878,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00129891,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47443,0.0609041,-7.38533e-05,4.75698e-08,-1.26207e-11,-110771,29.0731], Tmin=(100,'K'), Tmax=(903.546,'K')), NASAPolynomial(coeffs=[9.75707,0.0242364,-1.29796e-05,2.65467e-09,-1.93108e-13,-112268,-10.0457], Tmin=(903.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-921.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + group(COFOO) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C(F)O[CH]O(5834)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-595.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,482,664,788,1296,1923,265.047,558.524,2496.93],'cm^-1')),
        HinderedRotor(inertia=(0.833188,'amu*angstrom^2'), symmetry=1, barrier=(19.1566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833213,'amu*angstrom^2'), symmetry=1, barrier=(19.1572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1153,'amu*angstrom^2'), symmetry=1, barrier=(25.517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81287,0.0464932,-5.47881e-05,3.13845e-08,-6.97397e-12,-71562.3,19.9356], Tmin=(100,'K'), Tmax=(1106.43,'K')), NASAPolynomial(coeffs=[12.1838,0.00899983,-3.95791e-06,7.57392e-10,-5.37145e-14,-73857.2,-31.1472], Tmin=(1106.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(COFOO) + radical(OCJO)"""),
)

species(
    label = 'O=C(F)O[C](O)C(=O)F(13349)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-960.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,611,648,830,1210,1753,482,664,788,1296,1923,270.9,271.436,271.87,2094.75],'cm^-1')),
        HinderedRotor(inertia=(0.213297,'amu*angstrom^2'), symmetry=1, barrier=(11.1881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0229948,'amu*angstrom^2'), symmetry=1, barrier=(27.816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0022676,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530524,'amu*angstrom^2'), symmetry=1, barrier=(27.8137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31257,0.0640913,-7.89508e-05,4.97746e-08,-1.27111e-11,-115394,29.56], Tmin=(100,'K'), Tmax=(943.77,'K')), NASAPolynomial(coeffs=[11.2535,0.0219584,-1.19856e-05,2.47105e-09,-1.8059e-13,-117270,-17.8236], Tmin=(943.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-960.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + group(COFOO) + radical(Cs_P)"""),
)

species(
    label = 'O=C(O)F(343)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=C(O)OC(=O)F(13320)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {8,D}
7 C u0 p0 c0 {1,S} {2,S} {4,D}
8 C u0 p0 c0 {5,D} {6,D}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-733.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,482,664,788,1296,1923,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3057,'amu*angstrom^2'), symmetry=1, barrier=(30.0205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30077,'amu*angstrom^2'), symmetry=1, barrier=(29.9072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29405,'amu*angstrom^2'), symmetry=1, barrier=(29.7528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71035,0.0858968,-0.000104397,5.28997e-08,-9.4263e-12,-87999.8,29.9055], Tmin=(100,'K'), Tmax=(1697.2,'K')), NASAPolynomial(coeffs=[28.1163,-0.0134995,9.16953e-06,-1.81224e-09,1.20977e-13,-93933,-117.421], Tmin=(1697.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + group(COFOO) + missing(Cdd-CdO2d)"""),
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
    E0 = (-516.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-231.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-149.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (80.0039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-218.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-145.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-6.84313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-226.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-138.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-164.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-244.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-234.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-138.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-189.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-342.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-148.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-137.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-311.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-202.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-192.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-174.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-250.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-212.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-444.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-284.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)OC(O)C(=O)F(13079)'],
    products = ['HF(38)', 'CO2(14)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.29418e+10,'s^-1'), n=0.745914, Ea=(113.426,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'O=C(F)OC(O)F(5365)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(371.724,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO(13)', 'O=C(F)C(O)OF(12539)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(91.1518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OF(482)', 'O=C=COC(=O)F(13334)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.37874e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd_HNd;H_OH]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)OOC(F)=CO(13335)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(107.353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)OC=C(F)OO(13336)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(78.1838,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC1O[C](F)OO[C]1F(13337)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)[C](O)OC(=O)F(13338)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)O[C](O)C(=O)F(13339)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)C([O])OC(=O)F(13340)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C(F)[C](O)O[C](O)F(13341)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(OC(=O)F)[C](O)F(13342)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)OC([O])C(=O)F(13343)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(O[C](O)F)C(=O)F(13344)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C(F)OC(O)=C(O)F(13345)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(162.644,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=[C]C(O)OC(=O)F(13346)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'O=[C]OC(O)C(=O)F(13347)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(=O)F(341)', 'O=C(F)[CH]O(1552)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CFO(51)', '[O]C(O)C(=O)F(12721)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OH(7)', 'O=C(F)[CH]OC(=O)F(2174)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]C(OC(=O)F)C(=O)F(13348)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CFO(51)', 'O=C(F)O[CH]O(5834)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'O=C(F)O[C](O)C(=O)F(13349)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(F)OC(O)C(=O)F(13079)'],
    products = ['O=C(O)F(343)', 'O=CC(=O)F(335)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.82035e+06,'s^-1'), n=1.67734, Ea=(185.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'O=C=C(O)OC(=O)F(13320)'],
    products = ['O=C(F)OC(O)C(=O)F(13079)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(194.502,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #3513',
    isomers = [
        'O=C(F)OC(O)C(=O)F(13079)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3513',
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

