species(
    label = '[O]OC(F)COOC=C(F)F(9173)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {11,D} {15,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-599.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.258434,'amu*angstrom^2'), symmetry=1, barrier=(5.94191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258321,'amu*angstrom^2'), symmetry=1, barrier=(5.93931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98355,'amu*angstrom^2'), symmetry=1, barrier=(45.6058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19569,'amu*angstrom^2'), symmetry=1, barrier=(27.4912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.835686,0.116689,-0.000185728,1.61849e-07,-5.6603e-11,-71902.4,37.6354], Tmin=(100,'K'), Tmax=(771.244,'K')), NASAPolynomial(coeffs=[11.6184,0.0424229,-2.24725e-05,4.4665e-09,-3.15485e-13,-73535.7,-17.3483], Tmin=(771.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = 'O=CC(F)F(990)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.595284,'amu*angstrom^2'), symmetry=1, barrier=(13.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93683,0.00748989,0.000222302,-1.39306e-06,2.75418e-09,-67109.3,9.32225], Tmin=(10,'K'), Tmax=(169.939,'K')), NASAPolynomial(coeffs=[3.97443,0.0203448,-1.24424e-05,3.60772e-09,-4.02215e-13,-67130.4,8.62375], Tmin=(169.939,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-557.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)C=O(2464)',
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
        HinderedRotor(inertia=(1.18318,'amu*angstrom^2'), symmetry=1, barrier=(27.2035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603053,'amu*angstrom^2'), symmetry=1, barrier=(13.8654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3728.91,'J/mol'), sigma=(5.84942,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.45 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15437,0.0440803,-5.80203e-05,4.31468e-08,-1.32604e-11,-39422,18.3048], Tmin=(100,'K'), Tmax=(786.661,'K')), NASAPolynomial(coeffs=[7.18331,0.0185099,-9.26421e-06,1.8289e-09,-1.29954e-13,-40213.3,-4.75028], Tmin=(786.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OC(F)OOC=C(F)F(11532)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-603.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,509,613,660,1171,1360,1414,3084,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.575252,'amu*angstrom^2'), symmetry=1, barrier=(13.2262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91683,'amu*angstrom^2'), symmetry=1, barrier=(44.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575678,'amu*angstrom^2'), symmetry=1, barrier=(13.236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91774,'amu*angstrom^2'), symmetry=1, barrier=(44.0926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.227546,0.102342,-0.000171372,1.50611e-07,-5.22296e-11,-72489.6,33.1771], Tmin=(100,'K'), Tmax=(792.234,'K')), NASAPolynomial(coeffs=[11.5004,0.0320827,-1.7431e-05,3.47157e-09,-2.44472e-13,-74001.2,-18.4842], Tmin=(792.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-603.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHOO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = 'HO2(11)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1064.4,1465.7,3224.93],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[C]COOC=C(F)F-2(5228)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,D}
9  C u0 p1 c0 {3,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-259.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.247795,'amu*angstrom^2'), symmetry=1, barrier=(5.69729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24739,'amu*angstrom^2'), symmetry=1, barrier=(5.68798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12363,'amu*angstrom^2'), symmetry=1, barrier=(48.8265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11606,'amu*angstrom^2'), symmetry=1, barrier=(48.6523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445252,0.0855724,-0.00013183,1.13645e-07,-3.97752e-11,-31091.6,30.8932], Tmin=(100,'K'), Tmax=(753.698,'K')), NASAPolynomial(coeffs=[9.35003,0.0332036,-1.74371e-05,3.46655e-09,-2.45428e-13,-32288.8,-8.5864], Tmin=(753.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CJ2_singlet-FCs)"""),
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
    label = '[O]OCOOC=C(F)F(6588)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-368.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(1.02793,'amu*angstrom^2'), symmetry=1, barrier=(23.6341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0299,'amu*angstrom^2'), symmetry=1, barrier=(23.6793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02724,'amu*angstrom^2'), symmetry=1, barrier=(23.6182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02933,'amu*angstrom^2'), symmetry=1, barrier=(23.6664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.14972,0.0831557,-9.3994e-05,5.02031e-08,-1.02737e-11,-44112,32.3834], Tmin=(100,'K'), Tmax=(1207.56,'K')), NASAPolynomial(coeffs=[21.5623,0.0112354,-4.65673e-06,8.82138e-10,-6.28789e-14,-49355.7,-76.4597], Tmin=(1207.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(CF)OOC=C(F)F(11533)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {6,S} {9,S} {12,S}
9  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {5,S} {11,D} {15,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-602.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.83435,0.116718,-0.000186488,1.6278e-07,-5.70381e-11,-72327,38.8546], Tmin=(100,'K'), Tmax=(767.451,'K')), NASAPolynomial(coeffs=[11.6896,0.0421224,-2.24736e-05,4.48122e-09,-3.17132e-13,-73974.9,-16.4634], Tmin=(767.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-602.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)COC(F)(F)C=O(11534)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {11,D}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
11 C u0 p0 c0 {6,D} {10,S} {15,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-953.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.871896,0.117102,-0.000187507,1.61578e-07,-5.54888e-11,-114509,35.8069], Tmin=(100,'K'), Tmax=(783.866,'K')), NASAPolynomial(coeffs=[12.4665,0.0394658,-2.06251e-05,4.06893e-09,-2.85847e-13,-116306,-23.4191], Tmin=(783.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-953.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = '[O]C(F)COOC=C(F)F(11535)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-596.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,2402.1],'cm^-1')),
        HinderedRotor(inertia=(0.483537,'amu*angstrom^2'), symmetry=1, barrier=(11.1175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83285,'amu*angstrom^2'), symmetry=1, barrier=(42.1409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.482805,'amu*angstrom^2'), symmetry=1, barrier=(11.1006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82795,'amu*angstrom^2'), symmetry=1, barrier=(42.0281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.111208,0.0985684,-0.000147163,1.23557e-07,-4.24053e-11,-71557.8,34.4222], Tmin=(100,'K'), Tmax=(732.72,'K')), NASAPolynomial(coeffs=[10.5083,0.0385176,-1.99756e-05,3.96518e-09,-2.81023e-13,-73058.2,-13.1278], Tmin=(732.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-596.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O=C[C](F)F(787)',
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
    label = 'FC1COOO1(1271)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-275.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27851,0.0247619,2.97546e-05,-6.83338e-08,3.25693e-11,-33006.5,13.5315], Tmin=(100,'K'), Tmax=(880.079,'K')), NASAPolynomial(coeffs=[14.6769,0.00466582,2.21299e-06,-6.61872e-10,4.91965e-14,-36592.9,-52.6769], Tmin=(880.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(Cs-CsOsHH) + ring(123trioxolane)"""),
)

species(
    label = 'FC=COOC=C(F)F(2505)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-459.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,194,682,905,1196,1383,3221,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.0836088,'amu*angstrom^2'), symmetry=1, barrier=(21.7354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22505,'amu*angstrom^2'), symmetry=1, barrier=(5.17434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58333,'amu*angstrom^2'), symmetry=1, barrier=(59.3959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918625,0.0747774,-9.67782e-05,7.10186e-08,-2.18264e-11,-55130.2,29.1101], Tmin=(100,'K'), Tmax=(781.67,'K')), NASAPolynomial(coeffs=[8.97024,0.0335749,-1.77109e-05,3.5833e-09,-2.58452e-13,-56388.9,-7.75092], Tmin=(781.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'FC1COO[CH]C(F)(F)OO1(11536)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
11 C u1 p0 c0 {7,S} {10,S} {15,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-704.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.101731,0.0776461,-5.12585e-05,1.28497e-08,-6.43161e-13,-84629.6,27.4386], Tmin=(100,'K'), Tmax=(1400.73,'K')), NASAPolynomial(coeffs=[22.2103,0.0286535,-1.45601e-05,2.88719e-09,-2.04346e-13,-92324.6,-92.879], Tmin=(1400.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFFO) + group(Cs-CsOsHH) + ring(Cyclooctane) + radical(CCsJOOC)"""),
)

species(
    label = 'F[C](F)C1OOCC(F)OO1(11537)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {7,S} {11,S} {15,S}
11 C u1 p0 c0 {2,S} {3,S} {10,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-689.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.478605,0.0562741,7.6659e-05,-1.76319e-07,8.12443e-11,-82756,33.5334], Tmin=(100,'K'), Tmax=(927.419,'K')), NASAPolynomial(coeffs=[40.4236,-0.011315,9.96621e-06,-1.85375e-09,1.08141e-13,-95022.7,-185.946], Tmin=(927.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-689.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCsFFH) + ring(Cycloheptane) + radical(Csj(Cs-O2sO2sH)(F1s)(F1s))"""),
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
    label = '[O]O[CH]COOC=C(F)F(11538)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {7,S} {13,S}
9  C u0 p0 c0 {4,S} {10,D} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-184.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,2411.99],'cm^-1')),
        HinderedRotor(inertia=(0.33572,'amu*angstrom^2'), symmetry=1, barrier=(7.71887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335123,'amu*angstrom^2'), symmetry=1, barrier=(7.70515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.077,'amu*angstrom^2'), symmetry=1, barrier=(47.7543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08689,'amu*angstrom^2'), symmetry=1, barrier=(47.9817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.510771,0.112423,-0.000194619,1.7921e-07,-6.38297e-11,-22089.1,36.4466], Tmin=(100,'K'), Tmax=(832.757,'K')), NASAPolynomial(coeffs=[8.84946,0.0420122,-2.19492e-05,4.27885e-09,-2.96497e-13,-22765.6,-1.69956], Tmin=(832.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[O]OC(F)COOC=[C]F(11539)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
9  C u0 p0 c0 {4,S} {10,D} {14,S}
10 C u1 p0 c0 {2,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-147.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.416954,'amu*angstrom^2'), symmetry=1, barrier=(9.58658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98404,'amu*angstrom^2'), symmetry=1, barrier=(45.6169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416126,'amu*angstrom^2'), symmetry=1, barrier=(9.56755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98666,'amu*angstrom^2'), symmetry=1, barrier=(45.6773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416656,'amu*angstrom^2'), symmetry=1, barrier=(9.57975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.669383,0.114579,-0.000194913,1.75414e-07,-6.16584e-11,-17623.6,36.3227], Tmin=(100,'K'), Tmax=(817.376,'K')), NASAPolynomial(coeffs=[10.7007,0.0393322,-2.08481e-05,4.10063e-09,-2.86134e-13,-18827.5,-12.2326], Tmin=(817.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = '[O]CC(F)O[O](2435)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-192.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.151867,'amu*angstrom^2'), symmetry=1, barrier=(3.49173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151695,'amu*angstrom^2'), symmetry=1, barrier=(3.48777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88497,0.0547742,-9.52073e-05,9.34154e-08,-3.51331e-11,-23131.5,20.6623], Tmin=(100,'K'), Tmax=(835.28,'K')), NASAPolynomial(coeffs=[3.66998,0.0280018,-1.44021e-05,2.80114e-09,-1.94133e-13,-22793.9,16.1775], Tmin=(835.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC=C(F)F(785)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-270.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.857306,'amu*angstrom^2'), symmetry=1, barrier=(19.7111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3441.78,'J/mol'), sigma=(5.35945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.60 K, Pc=50.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84,0.010514,0.000109777,-2.90322e-07,2.13215e-10,-32466.3,11.2834], Tmin=(10,'K'), Tmax=(492.518,'K')), NASAPolynomial(coeffs=[6.26831,0.0212777,-1.585e-05,5.39864e-09,-6.82851e-13,-33075.2,-2.46567], Tmin=(492.518,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-270.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)O[O](216)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-48.0821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.486055,'amu*angstrom^2'), symmetry=1, barrier=(11.1754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485983,'amu*angstrom^2'), symmetry=1, barrier=(11.1737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3366.65,'J/mol'), sigma=(5.62898,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.86 K, Pc=42.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33352,0.040021,-5.12824e-05,3.75377e-08,-1.13901e-11,-5725.96,17.2582], Tmin=(100,'K'), Tmax=(795.131,'K')), NASAPolynomial(coeffs=[6.83717,0.0173642,-8.53954e-06,1.69956e-09,-1.21736e-13,-6442.14,-3.43674], Tmin=(795.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'CF2CH(73)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-92.0165,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([43.4896,46.129,89.6186],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([420.994,524.562,552.014,640.304,746.121,978.508,1261.79,1779.04,3347.18],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95511,0.00251604,4.75146e-05,-8.97685e-08,4.98372e-11,-11065.2,8.58309], Tmin=(10,'K'), Tmax=(611.191,'K')), NASAPolynomial(coeffs=[3.6485,0.0155109,-1.13452e-05,3.84899e-09,-4.87703e-13,-11233,8.2324], Tmin=(611.191,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-92.0165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[CH]DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OCC(F)O[O](9752)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u1 p2 c0 {2,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-195.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,180],'cm^-1')),
        HinderedRotor(inertia=(0.428614,'amu*angstrom^2'), symmetry=1, barrier=(9.85467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.814622,'amu*angstrom^2'), symmetry=1, barrier=(18.7298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427176,'amu*angstrom^2'), symmetry=1, barrier=(9.82161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11101,0.0708361,-0.000118584,1.05571e-07,-3.64788e-11,-23366.7,24.1017], Tmin=(100,'K'), Tmax=(841.837,'K')), NASAPolynomial(coeffs=[7.9624,0.0247539,-1.237e-05,2.36921e-09,-1.62541e-13,-24040.9,-4.92563], Tmin=(841.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + radical(ROOJ) + radical(ROOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]COOC=C(F)F(5259)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {6,S} {13,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-432.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,1561.66],'cm^-1')),
        HinderedRotor(inertia=(0.459853,'amu*angstrom^2'), symmetry=1, barrier=(10.5729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460675,'amu*angstrom^2'), symmetry=1, barrier=(10.5918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75683,'amu*angstrom^2'), symmetry=1, barrier=(40.3929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75281,'amu*angstrom^2'), symmetry=1, barrier=(40.3006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3300.81,'J/mol'), sigma=(5.3802,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=515.58 K, Pc=48.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.063971,0.095483,-0.000152789,1.35076e-07,-4.77761e-11,-51837,31.9185], Tmin=(100,'K'), Tmax=(783.142,'K')), NASAPolynomial(coeffs=[9.50595,0.0364947,-1.91911e-05,3.80077e-09,-2.67754e-13,-52985.8,-9.21846], Tmin=(783.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-432.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = '[O]O[CH]F(287)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-10.1775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(0.627598,'amu*angstrom^2'), symmetry=1, barrier=(14.4297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87617,0.0275723,-4.75711e-05,4.25871e-08,-1.46862e-11,-1186.35,12.8044], Tmin=(100,'K'), Tmax=(836.078,'K')), NASAPolynomial(coeffs=[5.8498,0.00835167,-4.12766e-06,8.02216e-10,-5.55856e-14,-1509.04,0.0345696], Tmin=(836.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.1775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[CH2]OOC=C(F)F(1854)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,D} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {5,D}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-206.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.221354,'amu*angstrom^2'), symmetry=1, barrier=(5.08937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127167,'amu*angstrom^2'), symmetry=1, barrier=(51.6655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25137,'amu*angstrom^2'), symmetry=1, barrier=(51.7634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69193,0.0559924,-6.72883e-05,4.66586e-08,-1.37484e-11,-24700.5,22.0546], Tmin=(100,'K'), Tmax=(808.78,'K')), NASAPolynomial(coeffs=[7.48746,0.0273299,-1.41303e-05,2.84201e-09,-2.04627e-13,-25638,-4.67565], Tmin=(808.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsJOOC)"""),
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
    label = '[O]OC(F)[CH]OOC=C(F)F(11540)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
9  C u1 p0 c0 {4,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-412.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.236855,'amu*angstrom^2'), symmetry=1, barrier=(5.44577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10945,'amu*angstrom^2'), symmetry=1, barrier=(48.5003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.105,'amu*angstrom^2'), symmetry=1, barrier=(48.3982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237239,'amu*angstrom^2'), symmetry=1, barrier=(5.45459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10993,'amu*angstrom^2'), symmetry=1, barrier=(48.5116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.860157,0.118939,-0.000201362,1.81027e-07,-6.3934e-11,-49479.6,38.0161], Tmin=(100,'K'), Tmax=(804.328,'K')), NASAPolynomial(coeffs=[11.167,0.0409716,-2.21023e-05,4.38486e-09,-3.07856e-13,-50827.1,-13.7383], Tmin=(804.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-412.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CCsJOOC)"""),
)

species(
    label = '[O]O[C](F)COOC=C(F)F(6594)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
9  C u1 p0 c0 {1,S} {6,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.309139,'amu*angstrom^2'), symmetry=1, barrier=(7.10772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309368,'amu*angstrom^2'), symmetry=1, barrier=(7.11298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18003,'amu*angstrom^2'), symmetry=1, barrier=(50.1232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309181,'amu*angstrom^2'), symmetry=1, barrier=(7.10868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17696,'amu*angstrom^2'), symmetry=1, barrier=(50.0525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.999278,0.125586,-0.000227432,2.11554e-07,-7.53593e-11,-48502.2,39.0225], Tmin=(100,'K'), Tmax=(836.803,'K')), NASAPolynomial(coeffs=[9.92279,0.0424552,-2.29891e-05,4.51974e-09,-3.13689e-13,-49247.5,-5.25535], Tmin=(836.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)COO[C]=C(F)F(11541)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {11,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {11,D}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-359.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,562,600,623,1070,1265,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70395,0.11942,-0.000216209,2.05444e-07,-7.48506e-11,-43078.1,40.2872], Tmin=(100,'K'), Tmax=(830.384,'K')), NASAPolynomial(coeffs=[7.83206,0.0459195,-2.49448e-05,4.92765e-09,-3.43565e-13,-43379.3,7.41479], Tmin=(830.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-359.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = '[O]OC=COOC=C(F)F(11542)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u0 p0 c0 {3,S} {10,D} {11,S}
9  C u0 p0 c0 {5,S} {7,D} {13,S}
10 C u0 p0 c0 {1,S} {2,S} {8,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-226.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.610394,'amu*angstrom^2'), symmetry=1, barrier=(14.0341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609864,'amu*angstrom^2'), symmetry=1, barrier=(14.022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05482,'amu*angstrom^2'), symmetry=1, barrier=(24.2523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609392,'amu*angstrom^2'), symmetry=1, barrier=(14.0111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0491018,0.0913249,-0.000120591,8.12826e-08,-2.17281e-11,-27108.5,35.4241], Tmin=(100,'K'), Tmax=(914.352,'K')), NASAPolynomial(coeffs=[15.1144,0.0254188,-1.24722e-05,2.45146e-09,-1.7429e-13,-29863.5,-35.9085], Tmin=(914.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)COOC#CF(11543)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
9  C u0 p0 c0 {5,S} {10,T}
10 C u0 p0 c0 {2,S} {9,T}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-110.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,2175,525,239,401,1367,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.216904,'amu*angstrom^2'), symmetry=1, barrier=(4.98705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65347,'amu*angstrom^2'), symmetry=1, barrier=(38.0165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65418,'amu*angstrom^2'), symmetry=1, barrier=(38.0329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664379,'amu*angstrom^2'), symmetry=1, barrier=(15.2754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.504557,0.110461,-0.000189155,1.71088e-07,-6.07614e-11,-13084.5,33.8364], Tmin=(100,'K'), Tmax=(799.75,'K')), NASAPolynomial(coeffs=[10.6894,0.0375381,-2.06179e-05,4.11892e-09,-2.90285e-13,-14333.4,-14.2809], Tmin=(799.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(CtCF) + radical(ROOJ)"""),
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
    label = 'O=C(F)COOC=C(F)F(11544)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-757.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.472029,0.0853558,-0.000109295,7.66413e-08,-2.22501e-11,-91016.1,31.709], Tmin=(100,'K'), Tmax=(828.563,'K')), NASAPolynomial(coeffs=[10.8279,0.0353607,-1.87845e-05,3.81444e-09,-2.75873e-13,-92732.2,-16.3042], Tmin=(828.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-757.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'OOC(F)[CH]OOC=C(F)F(11545)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {4,S} {15,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {5,S} {8,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-564.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04198,0.121515,-0.00019296,1.66954e-07,-5.83052e-11,-67753.5,37.9137], Tmin=(100,'K'), Tmax=(746.155,'K')), NASAPolynomial(coeffs=[12.4529,0.0431748,-2.34186e-05,4.70461e-09,-3.34675e-13,-69600.4,-22.1212], Tmin=(746.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOOC)"""),
)

species(
    label = 'OOC(F)COO[C]=C(F)F(11546)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {11,S}
7  O u0 p2 c0 {5,S} {15,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {11,D}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-511.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,562,600,623,1070,1265,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.943093,0.122736,-0.000210697,1.95521e-07,-7.11422e-11,-61349.6,40.3862], Tmin=(100,'K'), Tmax=(803.254,'K')), NASAPolynomial(coeffs=[9.30316,0.0477867,-2.60581e-05,5.19781e-09,-3.66168e-13,-62223.8,-1.99675], Tmin=(803.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C=CJO)"""),
)

species(
    label = 'FOO[CH]COOC=C(F)F(11547)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {3,S} {6,S}
8  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 C u0 p0 c0 {5,S} {11,D} {15,S}
11 C u0 p0 c0 {1,S} {2,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-242.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(1.73054,'amu*angstrom^2'), symmetry=1, barrier=(39.7886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73092,'amu*angstrom^2'), symmetry=1, barrier=(39.7972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73054,'amu*angstrom^2'), symmetry=1, barrier=(39.7884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73051,'amu*angstrom^2'), symmetry=1, barrier=(39.7879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575267,'amu*angstrom^2'), symmetry=1, barrier=(13.2265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09211,0.124785,-0.000209906,1.89368e-07,-6.72509e-11,-29024.3,39.2288], Tmin=(100,'K'), Tmax=(802.801,'K')), NASAPolynomial(coeffs=[11.0308,0.0448612,-2.40971e-05,4.77827e-09,-3.35595e-13,-30341.7,-12.6764], Tmin=(802.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOO)"""),
)

species(
    label = 'F[C]=COOCC(F)OOF(11548)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {3,S} {5,S}
8  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {6,S} {11,D} {15,S}
11 C u1 p0 c0 {2,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-204.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,3010,987.5,1337.5,450,1655,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(1.28107,'amu*angstrom^2'), symmetry=1, barrier=(29.4544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28157,'amu*angstrom^2'), symmetry=1, barrier=(29.4657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28128,'amu*angstrom^2'), symmetry=1, barrier=(29.4591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2818,'amu*angstrom^2'), symmetry=1, barrier=(29.4711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28149,'amu*angstrom^2'), symmetry=1, barrier=(29.464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25446,0.127388,-0.0002115,1.86845e-07,-6.55056e-11,-24422.4,39.0455], Tmin=(100,'K'), Tmax=(781.479,'K')), NASAPolynomial(coeffs=[12.7809,0.0424008,-2.31391e-05,4.62984e-09,-3.27377e-13,-26214.6,-22.638], Tmin=(781.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    E0 = (-312.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (133.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (36.7082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (145.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-51.1711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-257.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-104.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-295.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-204.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-301.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-254.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (137.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (174.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-316.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-69.0653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-38.0801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-191.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (32.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (48.0743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (56.1878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (101.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-38.7574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (127.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-192.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-233.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-153.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (94.5679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (100.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['O=CC(F)F(990)', '[O]OC(F)C=O(2464)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(37.8783,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[O]OC(F)OOC=C(F)F(11532)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.26474e-05,'m^3/(mol*s)'), n=2.95311, Ea=(69.4891,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['HO2(11)', 'F[C]COOC=C(F)F-2(5228)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(44.7096,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(40)', '[O]OCOOC=C(F)F(6588)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.12553e-07,'m^3/(mol*s)'), n=3.34134, Ea=(125.462,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['[O]OC(CF)OOC=C(F)F(11533)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['[O]OC(F)COC(F)(F)C=O(11534)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(92.6483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', '[O]C(F)COOC=C(F)F(11535)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['O=C[C](F)F(787)', 'FC1COOO1(1271)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['HO2(11)', 'FC=COOC=C(F)F(2505)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.00406e+10,'s^-1'), n=0.563333, Ea=(145.788,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['FC1COO[CH]C(F)(F)OO1(11536)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.61579e+10,'s^-1'), n=0.429849, Ea=(48.6953,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['F[C](F)C1OOCC(F)OO1(11537)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;doublebond_intra;radadd_intra_O] + [R8;doublebond_intra;radadd_intra] for rate rule [R8;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]O[CH]COOC=C(F)F(11538)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC(F)COOC=[C]F(11539)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C[C](F)F(787)', '[O]CC(F)O[O](2435)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(19.0511,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC=C(F)F(785)', '[CH2]C(F)O[O](216)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2CH(73)', '[O]OCC(F)O[O](9752)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', 'F[CH]COOC=C(F)F(5259)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]O[CH]F(287)', '[CH2]OOC=C(F)F(1854)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]OC(F)[CH]OOC=C(F)F(11540)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]O[C](F)COOC=C(F)F(6594)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]OC(F)COO[C]=C(F)F(11541)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]OC=COOC=C(F)F(11542)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(219.876,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]OC(F)COOC#CF(11543)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(269.181,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['OH(7)', 'O=C(F)COOC=C(F)F(11544)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OC(F)COOC=C(F)F(9173)'],
    products = ['OOC(F)[CH]OOC=C(F)F(11545)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.94e+12,'s^-1'), n=-0.24, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 291 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_OOH/H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OOC(F)COO[C]=C(F)F(11546)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.81284e+06,'s^-1'), n=1.31318, Ea=(108.582,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R7H_OOCs4;Y_rad_out;XH_out] for rate rule [R7H_OOCs4;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.1622776601683795
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['FOO[CH]COOC=C(F)F(11547)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(88.2866,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C]=COOCC(F)OOF(11548)'],
    products = ['[O]OC(F)COOC=C(F)F(9173)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(56.2044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #3061',
    isomers = [
        '[O]OC(F)COOC=C(F)F(9173)',
    ],
    reactants = [
        ('O=CC(F)F(990)', '[O]OC(F)C=O(2464)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3061',
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

