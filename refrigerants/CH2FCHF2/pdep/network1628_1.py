species(
    label = '[O]OC(F)=C(F)OOC(F)F(5346)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-766.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,262,390,483,597,572,732,631,807,1275,1439,180],'cm^-1')),
        HinderedRotor(inertia=(0.269966,'amu*angstrom^2'), symmetry=1, barrier=(6.20705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269808,'amu*angstrom^2'), symmetry=1, barrier=(6.20341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269923,'amu*angstrom^2'), symmetry=1, barrier=(6.20607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26982,'amu*angstrom^2'), symmetry=1, barrier=(6.2037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.344756,0.115749,-0.000232504,2.31514e-07,-8.57461e-11,-92086.3,37.0597], Tmin=(100,'K'), Tmax=(849.658,'K')), NASAPolynomial(coeffs=[5.1131,0.044092,-2.48574e-05,4.92034e-09,-3.41004e-13,-91354.7,21.381], Tmin=(849.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-766.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)C(=O)F(4704)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-582.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.519617,'amu*angstrom^2'), symmetry=1, barrier=(11.947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54035,'amu*angstrom^2'), symmetry=1, barrier=(35.4156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3788.12,'J/mol'), sigma=(5.7366,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.70 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66217,0.0558669,-8.74757e-05,7.36244e-08,-2.48646e-11,-69953.7,20.9883], Tmin=(100,'K'), Tmax=(757.994,'K')), NASAPolynomial(coeffs=[8.40288,0.0184184,-9.65371e-06,1.91163e-09,-1.34843e-13,-70921.6,-9.30828], Tmin=(757.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=C(F)OO(5900)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {3,S} {9,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {3,S} {8,D}
8 C u0 p0 c0 {2,S} {4,S} {7,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-319.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439],'cm^-1')),
        HinderedRotor(inertia=(0.225341,'amu*angstrom^2'), symmetry=1, barrier=(5.18104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226784,'amu*angstrom^2'), symmetry=1, barrier=(5.2142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225602,'amu*angstrom^2'), symmetry=1, barrier=(5.18703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01178,0.0821821,-0.000172851,1.75775e-07,-6.54256e-11,-38387.2,27.38], Tmin=(100,'K'), Tmax=(862.066,'K')), NASAPolynomial(coeffs=[3.3526,0.0322772,-1.80803e-05,3.54822e-09,-2.43742e-13,-37340,24.8489], Tmin=(862.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=C(F)OOF(5901)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {9,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,S} {9,S}
6 O u0 p2 c0 {3,S} {4,S}
7 O u1 p2 c0 {5,S}
8 C u0 p0 c0 {2,S} {4,S} {9,D}
9 C u0 p0 c0 {1,S} {5,S} {8,D}
"""),
    E0 = (-224.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,277,555,632,262,390,483,597,572,732,631,807,1275,1439,180,2414.5],'cm^-1')),
        HinderedRotor(inertia=(0.197565,'amu*angstrom^2'), symmetry=1, barrier=(4.54242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199104,'amu*angstrom^2'), symmetry=1, barrier=(4.57778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201256,'amu*angstrom^2'), symmetry=1, barrier=(4.62727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648507,0.0918965,-0.000195802,1.98327e-07,-7.35225e-11,-26913.8,30.0649], Tmin=(100,'K'), Tmax=(861.086,'K')), NASAPolynomial(coeffs=[4.02106,0.0333715,-1.91934e-05,3.79155e-09,-2.61048e-13,-25905.7,23.5249], Tmin=(861.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(OC(F)F)C(=O)F(5902)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {11,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {5,S} {12,S}
11 C u0 p0 c0 {4,S} {7,D} {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1216.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.691372,0.116417,-0.0002127,1.93486e-07,-6.72831e-11,-146189,34.0712], Tmin=(100,'K'), Tmax=(840.044,'K')), NASAPolynomial(coeffs=[11.7779,0.0321446,-1.77637e-05,3.50037e-09,-2.42566e-13,-147405,-18.6838], Tmin=(840.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1216.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFHO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = 'O=C(F)[C](F)OOC(F)F(5903)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
9  C u1 p0 c0 {3,S} {6,S} {10,S}
10 C u0 p0 c0 {4,S} {7,D} {9,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-1041.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,519,593,830,1115,1166,1389,1413,3114,280,501,1494,1531,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.43221,'amu*angstrom^2'), symmetry=1, barrier=(55.9212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42979,'amu*angstrom^2'), symmetry=1, barrier=(55.8657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75609,'amu*angstrom^2'), symmetry=1, barrier=(17.384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43839,'amu*angstrom^2'), symmetry=1, barrier=(56.0634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16074,0.0755925,-7.47332e-05,-2.25781e-08,6.81657e-11,-125226,27.2447], Tmin=(100,'K'), Tmax=(478.495,'K')), NASAPolynomial(coeffs=[9.04283,0.0336711,-1.84566e-05,3.7031e-09,-2.62451e-13,-126255,-7.83916], Tmin=(478.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1041.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFHO) + group(COCsFO) + radical(CsCOF1sO2s)"""),
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
    label = 'FC1=C(F)OOO1(5904)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u0 p0 c0 {1,S} {4,S} {6,D}
"""),
    E0 = (-242.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11875,0.00894401,3.87226e-05,-5.40041e-08,1.99784e-11,-29142.3,19.1533], Tmin=(100,'K'), Tmax=(1010.02,'K')), NASAPolynomial(coeffs=[8.91177,0.0124164,-5.66293e-06,1.18561e-09,-9.10774e-14,-31659.8,-15.5223], Tmin=(1010.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CdCFO) + group(CdCFO) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + ring(Cyclopentane)"""),
)

species(
    label = 'F[C]1OOC1(F)OOC(F)F(5905)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {5,S} {10,S}
8  O u0 p2 c0 {6,S} {11,S}
9  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {8,S} {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-828.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282386,0.103354,-0.000169844,1.48128e-07,-5.16231e-11,-99490.1,34.0574], Tmin=(100,'K'), Tmax=(756.744,'K')), NASAPolynomial(coeffs=[11.7649,0.0333388,-1.85006e-05,3.73494e-09,-2.65846e-13,-101132,-19.507], Tmin=(756.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-828.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + group(CsFFHO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](OOC(F)F)C1(F)OO1(5906)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {7,S} {11,S}
9  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {8,S} {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-839.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.324927,0.104114,-0.000169526,1.4637e-07,-5.05842e-11,-100776,33.5364], Tmin=(100,'K'), Tmax=(750.464,'K')), NASAPolynomial(coeffs=[12.0902,0.0332164,-1.83729e-05,3.70481e-09,-2.63642e-13,-102506,-21.9088], Tmin=(750.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-839.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFHO) + ring(O2s-O2s-Cs(C-F)) + radical(CsCsF1sO2s)"""),
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
    label = '[O]OC(F)=C(F)OO[CH]F(5907)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {5,S} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-335.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.258236,'amu*angstrom^2'), symmetry=1, barrier=(5.93735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257419,'amu*angstrom^2'), symmetry=1, barrier=(5.91856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258156,'amu*angstrom^2'), symmetry=1, barrier=(5.93552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257813,'amu*angstrom^2'), symmetry=1, barrier=(5.92764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0188235,0.109509,-0.00022864,2.31376e-07,-8.60961e-11,-40201.3,35.175], Tmin=(100,'K'), Tmax=(856.992,'K')), NASAPolynomial(coeffs=[3.85121,0.041977,-2.38525e-05,4.71188e-09,-3.25223e-13,-39048.1,27.7002], Tmin=(856.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]O[C]=C(F)OOC(F)F(5908)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-339.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370,180,1709.2],'cm^-1')),
        HinderedRotor(inertia=(0.282711,'amu*angstrom^2'), symmetry=1, barrier=(6.50007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282933,'amu*angstrom^2'), symmetry=1, barrier=(6.50519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282511,'amu*angstrom^2'), symmetry=1, barrier=(6.49548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281884,'amu*angstrom^2'), symmetry=1, barrier=(6.48108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270087,0.102107,-0.000211395,2.15721e-07,-8.0992e-11,-40773.3,37.2627], Tmin=(100,'K'), Tmax=(854.549,'K')), NASAPolynomial(coeffs=[2.8954,0.0425521,-2.38894e-05,4.71413e-09,-3.25785e-13,-39496.2,35.1076], Tmin=(854.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]OC(F)=[C]OOC(F)F(5909)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-339.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370,180,1709.2],'cm^-1')),
        HinderedRotor(inertia=(0.282711,'amu*angstrom^2'), symmetry=1, barrier=(6.50007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282933,'amu*angstrom^2'), symmetry=1, barrier=(6.50519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282511,'amu*angstrom^2'), symmetry=1, barrier=(6.49548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281884,'amu*angstrom^2'), symmetry=1, barrier=(6.48108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270087,0.102107,-0.000211395,2.15721e-07,-8.0992e-11,-40773.3,37.2627], Tmin=(100,'K'), Tmax=(854.549,'K')), NASAPolynomial(coeffs=[2.8954,0.0425521,-2.38894e-05,4.71413e-09,-3.25785e-13,-39496.2,35.1076], Tmin=(854.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]O[C](F)C(=O)F(4794)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
"""),
    E0 = (-443.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,280,501,1494,1531,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(4.11501,'amu*angstrom^2'), symmetry=1, barrier=(94.6122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.22048,'amu*angstrom^2'), symmetry=1, barrier=(97.0371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93926,0.0519595,-9.19904e-05,8.61301e-08,-3.09897e-11,-53220.8,20.2084], Tmin=(100,'K'), Tmax=(837.273,'K')), NASAPolynomial(coeffs=[5.80387,0.0200068,-1.05786e-05,2.06414e-09,-1.42885e-13,-53395.1,5.07392], Tmin=(837.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
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
    label = '[O]OC(F)=C(F)O[O](5910)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u1 p2 c0 {3,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {2,S} {3,S} {8,D}
8 C u0 p0 c0 {1,S} {4,S} {7,D}
"""),
    E0 = (-167.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,262,390,483,597,572,732,631,807,1275,1439],'cm^-1')),
        HinderedRotor(inertia=(0.10821,'amu*angstrom^2'), symmetry=1, barrier=(2.48797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108861,'amu*angstrom^2'), symmetry=1, barrier=(2.50293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27952,0.0785304,-0.000177246,1.84454e-07,-6.87705e-11,-20116.9,26.4856], Tmin=(100,'K'), Tmax=(877.806,'K')), NASAPolynomial(coeffs=[1.70016,0.0307295,-1.71559e-05,3.32351e-09,-2.24961e-13,-18423,34.5804], Tmin=(877.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)F(2774)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-427.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114],'cm^-1')),
        HinderedRotor(inertia=(0.586675,'amu*angstrom^2'), symmetry=1, barrier=(13.4888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3217.18,'J/mol'), sigma=(5.19785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.52 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92725,0.00470123,7.24333e-05,-1.66658e-07,1.13367e-10,-51470.2,11.559], Tmin=(10,'K'), Tmax=(497.316,'K')), NASAPolynomial(coeffs=[3.77619,0.0196241,-1.39226e-05,4.52949e-09,-5.50671e-13,-51624.7,10.478], Tmin=(497.316,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-427.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)=[C]F(4600)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u1 p0 c0 {2,S} {5,D}
"""),
    E0 = (0.671019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,293,496,537,1218,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.147304,'amu*angstrom^2'), symmetry=1, barrier=(3.38682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94084,0.0545835,-0.000112998,1.10274e-07,-3.96452e-11,145.962,18.1125], Tmin=(100,'K'), Tmax=(865.493,'K')), NASAPolynomial(coeffs=[5.76978,0.015686,-8.83972e-06,1.73986e-09,-1.19296e-13,277.257,4.78062], Tmin=(865.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.671019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,-1038.59,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,-1040.82,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C]=C(F)OOC(F)F(3479)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-571.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.29546,'amu*angstrom^2'), symmetry=1, barrier=(6.7932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0038354,'amu*angstrom^2'), symmetry=1, barrier=(6.82622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30504,'amu*angstrom^2'), symmetry=1, barrier=(52.9973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631788,0.0865961,-0.000162241,1.57473e-07,-5.82917e-11,-68630.4,29.3763], Tmin=(100,'K'), Tmax=(825.957,'K')), NASAPolynomial(coeffs=[6.21224,0.0336312,-1.89447e-05,3.78877e-09,-2.658e-13,-68667.4,8.8771], Tmin=(825.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-571.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    label = '[O]OC(F)=C(F)OO[C](F)F(5911)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {1,S} {7,S} {9,D}
11 C u1 p0 c0 {3,S} {4,S} {6,S}
"""),
    E0 = (-548.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.241434,'amu*angstrom^2'), symmetry=1, barrier=(5.55105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242456,'amu*angstrom^2'), symmetry=1, barrier=(5.57454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242283,'amu*angstrom^2'), symmetry=1, barrier=(5.57056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24156,'amu*angstrom^2'), symmetry=1, barrier=(5.55393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.281864,0.113326,-0.000228726,2.2517e-07,-8.24215e-11,-65866.5,37.3324], Tmin=(100,'K'), Tmax=(853.107,'K')), NASAPolynomial(coeffs=[6.29383,0.0392179,-2.23305e-05,4.41801e-09,-3.0548e-13,-65413.7,15.883], Tmin=(853.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
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
    label = '[O]OC#COOC(F)F(5912)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {5,S} {8,T}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-130.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(1.45551,'amu*angstrom^2'), symmetry=1, barrier=(33.465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46971,'amu*angstrom^2'), symmetry=1, barrier=(33.7915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05141,'amu*angstrom^2'), symmetry=1, barrier=(24.1741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43879,'amu*angstrom^2'), symmetry=1, barrier=(33.0807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387015,0.0879023,-0.000151357,1.34243e-07,-4.65445e-11,-15557.8,30.9166], Tmin=(100,'K'), Tmax=(806.533,'K')), NASAPolynomial(coeffs=[10.4902,0.026008,-1.43219e-05,2.85074e-09,-2.00109e-13,-16804.1,-13.2762], Tmin=(806.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(CsFFHO) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)=C(F)OO[C](F)F(5913)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {7,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {1,S} {7,S} {9,D}
11 C u1 p0 c0 {3,S} {4,S} {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-700.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,262,390,483,597,572,732,631,807,1275,1439,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.184858,'amu*angstrom^2'), symmetry=1, barrier=(4.25026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.78291,'amu*angstrom^2'), symmetry=1, barrier=(63.9845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184613,'amu*angstrom^2'), symmetry=1, barrier=(4.24462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185684,'amu*angstrom^2'), symmetry=1, barrier=(4.26923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77592,'amu*angstrom^2'), symmetry=1, barrier=(63.8239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.539668,0.116867,-0.000223996,2.16189e-07,-7.90414e-11,-84137.2,37.4978], Tmin=(100,'K'), Tmax=(834.498,'K')), NASAPolynomial(coeffs=[7.86532,0.0409069,-2.3338e-05,4.6626e-09,-3.25927e-13,-84297.9,5.91141], Tmin=(834.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-700.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'F[CH]OOC(F)=C(F)OOF(5914)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {1,S} {7,S} {9,D}
11 C u1 p0 c0 {3,S} {6,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-391.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,262,390,483,597,572,732,631,807,1275,1439,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(0.545186,'amu*angstrom^2'), symmetry=1, barrier=(12.5349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94934,'amu*angstrom^2'), symmetry=1, barrier=(44.8192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.541578,'amu*angstrom^2'), symmetry=1, barrier=(12.4519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.644116,0.122815,-0.000247042,2.45175e-07,-9.09006e-11,-46998.4,38.0403], Tmin=(100,'K'), Tmax=(841.841,'K')), NASAPolynomial(coeffs=[6.11186,0.0447237,-2.59515e-05,5.1946e-09,-3.62539e-13,-46506.2,16.2889], Tmin=(841.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-391.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sHO2s)"""),
)

species(
    label = 'FOOC(F)=[C]OOC(F)F(5915)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {4,S} {6,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-396.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.344111,'amu*angstrom^2'), symmetry=1, barrier=(7.91178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343983,'amu*angstrom^2'), symmetry=1, barrier=(7.90884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12929,'amu*angstrom^2'), symmetry=1, barrier=(48.9567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12596,'amu*angstrom^2'), symmetry=1, barrier=(48.8799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.354608,0.115408,-0.000229786,2.29526e-07,-8.58125e-11,-47570.4,40.1257], Tmin=(100,'K'), Tmax=(838.702,'K')), NASAPolynomial(coeffs=[5.14833,0.0453123,-2.59963e-05,5.19876e-09,-3.63259e-13,-46951.2,23.7396], Tmin=(838.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'FOO[C]=C(F)OOC(F)F(5916)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-396.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.344111,'amu*angstrom^2'), symmetry=1, barrier=(7.91178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343983,'amu*angstrom^2'), symmetry=1, barrier=(7.90884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12929,'amu*angstrom^2'), symmetry=1, barrier=(48.9567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12596,'amu*angstrom^2'), symmetry=1, barrier=(48.8799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.354608,0.115408,-0.000229786,2.29526e-07,-8.58125e-11,-47570.4,40.1257], Tmin=(100,'K'), Tmax=(838.702,'K')), NASAPolynomial(coeffs=[5.14833,0.0453123,-2.59963e-05,5.19876e-09,-3.63259e-13,-46951.2,23.7396], Tmin=(838.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-420.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-119.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (260.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-363.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-432.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-334.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-305.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-385.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (72.2712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (67.5939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (67.5939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-432.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-89.9714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-92.6162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-245.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-2.30268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (195.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-303.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (3.60107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (0.306171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (22.7386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    products = ['CF2O(49)', '[O]OC(F)C(=O)F(4704)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(11.7801,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', '[O]OC(F)=C(F)OO(5900)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(69.7321,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[O]OC(F)=C(F)OOF(5901)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(12.0472,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    products = ['[O]OC(F)(OC(F)F)C(=O)F(5902)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(68.3004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[C](F)OOC(F)F(5903)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(32.1045,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 32.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    products = ['[O]C(F)F(364)', 'FC1=C(F)OOO1(5904)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(97.9041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 93.7 to 97.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    products = ['F[C]1OOC1(F)OOC(F)F(5905)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    products = ['F[C](OOC(F)F)C1(F)OO1(5906)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O]OC(F)=C(F)OO[CH]F(5907)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]O[C]=C(F)OOC(F)F(5908)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]OC(F)=[C]OOC(F)F(5909)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)F(364)', '[O]O[C](F)C(=O)F(4794)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(102.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Ea raised from 102.1 to 102.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHF2(82)', '[O]OC(F)=C(F)O[O](5910)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC(F)F(2774)', '[O]OC(F)=[C]F(4600)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', 'F[C]=C(F)OOC(F)F(3479)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(F)=C(F)OO[C](F)F(5911)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(78)', '[O]OC#COOC(F)F(5912)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OOC(F)=C(F)OO[C](F)F(5913)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2000,'s^-1'), n=1.9, Ea=(62.3416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H_OOCs4;C_rad_out_single;XH_out] for rate rule [R7H_OOCs4;C_rad_out_noH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]OOC(F)=C(F)OOF(5914)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(60.9255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FOOC(F)=[C]OOC(F)F(5915)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.3079,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FOO[C]=C(F)OOC(F)F(5916)'],
    products = ['[O]OC(F)=C(F)OOC(F)F(5346)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(84.7403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1628',
    isomers = [
        '[O]OC(F)=C(F)OOC(F)F(5346)',
    ],
    reactants = [
        ('CF2O(49)', '[O]OC(F)C(=O)F(4704)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1628',
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

