species(
    label = '[O]OC(F)C(F)OOC=C(F)F(12368)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
11 C u0 p0 c0 {6,S} {12,D} {15,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-810.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.966844,'amu*angstrom^2'), symmetry=1, barrier=(22.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84221,'amu*angstrom^2'), symmetry=1, barrier=(42.356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83501,'amu*angstrom^2'), symmetry=1, barrier=(42.1904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284641,'amu*angstrom^2'), symmetry=1, barrier=(6.54446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83613,'amu*angstrom^2'), symmetry=1, barrier=(42.2163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.891078,0.117905,-0.000170658,1.31984e-07,-4.14226e-11,-97281.7,38.3997], Tmin=(100,'K'), Tmax=(775.111,'K')), NASAPolynomial(coeffs=[13.941,0.041363,-2.25311e-05,4.57969e-09,-3.29893e-13,-99581,-29.3777], Tmin=(775.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-810.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(=O)F(5430)',
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
        HinderedRotor(inertia=(0.519618,'amu*angstrom^2'), symmetry=1, barrier=(11.947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54034,'amu*angstrom^2'), symmetry=1, barrier=(35.4154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3788.12,'J/mol'), sigma=(5.7366,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.70 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66212,0.0558676,-8.74785e-05,7.3629e-08,-2.48671e-11,-69953.7,20.9885], Tmin=(100,'K'), Tmax=(757.894,'K')), NASAPolynomial(coeffs=[8.40295,0.0184183,-9.65364e-06,1.91161e-09,-1.34842e-13,-70921.6,-9.30864], Tmin=(757.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = 'F[C]C(F)OOC=C(F)F-2(4623)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 C u0 p1 c0 {4,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-485.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.41161,'amu*angstrom^2'), symmetry=1, barrier=(32.4556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40911,'amu*angstrom^2'), symmetry=1, barrier=(32.3983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39839,'amu*angstrom^2'), symmetry=1, barrier=(32.1518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40578,'amu*angstrom^2'), symmetry=1, barrier=(32.3217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240842,0.0909294,-0.000132507,1.04149e-07,-3.33613e-11,-58231.3,31.8149], Tmin=(100,'K'), Tmax=(758.798,'K')), NASAPolynomial(coeffs=[11.209,0.033107,-1.81957e-05,3.70994e-09,-2.67631e-13,-59895.7,-18.0719], Tmin=(758.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]OC(OOC=C(F)F)C(F)F(12433)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {9,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {7,S} {10,S} {13,S}
10 C u0 p0 c0 {1,S} {2,S} {9,S} {14,S}
11 C u0 p0 c0 {6,S} {12,D} {15,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-824.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29667,0.128348,-0.000213566,1.88393e-07,-6.58568e-11,-98952.8,41.2993], Tmin=(100,'K'), Tmax=(783.635,'K')), NASAPolynomial(coeffs=[13.0213,0.0419987,-2.28894e-05,4.57621e-09,-3.23347e-13,-100790,-21.687], Tmin=(783.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-824.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(F)OC(F)(F)C=O(12434)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {12,S}
12 C u0 p0 c0 {7,D} {11,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1164.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.953272,0.118679,-0.00017399,1.34239e-07,-4.16672e-11,-139887,36.6613], Tmin=(100,'K'), Tmax=(786.052,'K')), NASAPolynomial(coeffs=[14.8407,0.038309,-2.06238e-05,4.16728e-09,-2.98979e-13,-142370,-35.7334], Tmin=(786.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1164.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = '[O]C(F)C(F)OOC=C(F)F(12435)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-807.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,391,562,707,872,1109,1210,1289,3137,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.108765,'amu*angstrom^2'), symmetry=1, barrier=(2.50073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17669,'amu*angstrom^2'), symmetry=1, barrier=(50.0464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874156,'amu*angstrom^2'), symmetry=1, barrier=(20.0986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17683,'amu*angstrom^2'), symmetry=1, barrier=(50.0497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.127658,0.0993424,-0.000130613,9.17948e-08,-2.63761e-11,-96938.8,35.0459], Tmin=(100,'K'), Tmax=(841.013,'K')), NASAPolynomial(coeffs=[12.9213,0.0372795,-1.99205e-05,4.0495e-09,-2.92906e-13,-99133.6,-25.6485], Tmin=(841.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-807.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(O2sj(Cs-CsF1sH))"""),
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
    label = 'FC1OOOC1F(4568)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-475.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9147,0.0309237,2.38484e-05,-6.7797e-08,3.32792e-11,-57125.8,14.9279], Tmin=(100,'K'), Tmax=(891.009,'K')), NASAPolynomial(coeffs=[17.316,0.00257078,2.91405e-06,-7.56524e-10,5.35491e-14,-61489.4,-66.6829], Tmin=(891.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(123trioxolane)"""),
)

species(
    label = 'FC=C(F)OOC=C(F)F(4618)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u0 p0 c0 {6,S} {10,D} {11,S}
9  C u0 p0 c0 {2,S} {7,D} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-646.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.0999628,'amu*angstrom^2'), symmetry=1, barrier=(2.29834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100081,'amu*angstrom^2'), symmetry=1, barrier=(2.30105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.90997,'amu*angstrom^2'), symmetry=1, barrier=(66.9059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.210419,0.0941671,-0.000160864,1.49957e-07,-5.4977e-11,-77609.3,32.3446], Tmin=(100,'K'), Tmax=(797.895,'K')), NASAPolynomial(coeffs=[7.78726,0.037901,-2.07179e-05,4.14413e-09,-2.92649e-13,-78236.4,1.14812], Tmin=(797.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-646.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'FC1OO[CH]C(F)(F)OOC1F(12436)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {6,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {11,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-905.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.475431,0.0839817,-5.79812e-05,1.46256e-08,-4.81526e-13,-108749,29.5591], Tmin=(100,'K'), Tmax=(1319.27,'K')), NASAPolynomial(coeffs=[23.8903,0.0279555,-1.45752e-05,2.94708e-09,-2.11936e-13,-116731,-100.63], Tmin=(1319.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-905.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(CsCFHO) + group(CsCFFO) + group(Cs-CsOsHH) + ring(Cyclooctane) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C1OOC(F)C(F)OO1(12437)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {7,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
11 C u0 p0 c0 {6,S} {8,S} {12,S} {15,S}
12 C u1 p0 c0 {3,S} {4,S} {11,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-890.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.838939,0.0623803,7.10301e-05,-1.76293e-07,8.22567e-11,-106875,35.6113], Tmin=(100,'K'), Tmax=(929.48,'K')), NASAPolynomial(coeffs=[43.0977,-0.0134683,1.07003e-05,-1.9561e-09,1.13126e-13,-119934,-199.456], Tmin=(929.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFHO) + group(CsCsFFH) + ring(Cycloheptane) + radical(Csj(Cs-O2sO2sH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[O]O[CH]C(F)OOC=C(F)F(12438)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-410.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.254099,'amu*angstrom^2'), symmetry=1, barrier=(5.84223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09375,'amu*angstrom^2'), symmetry=1, barrier=(48.1394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254151,'amu*angstrom^2'), symmetry=1, barrier=(5.84344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09533,'amu*angstrom^2'), symmetry=1, barrier=(48.1758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0938,'amu*angstrom^2'), symmetry=1, barrier=(48.1406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.898463,0.120242,-0.000205505,1.85677e-07,-6.56831e-11,-49221.1,38.0076], Tmin=(100,'K'), Tmax=(808.861,'K')), NASAPolynomial(coeffs=[11.0733,0.0412441,-2.22985e-05,4.42161e-09,-3.101e-13,-50510.2,-13.2069], Tmin=(808.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CCsJOOH)"""),
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
    label = '[O]OC(F)C(F)OOC=[C]F(12439)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-358.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.73656,'amu*angstrom^2'), symmetry=1, barrier=(39.9269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553869,'amu*angstrom^2'), symmetry=1, barrier=(12.7345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553772,'amu*angstrom^2'), symmetry=1, barrier=(12.7323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73434,'amu*angstrom^2'), symmetry=1, barrier=(39.8759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73657,'amu*angstrom^2'), symmetry=1, barrier=(39.9273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.941983,0.118718,-0.000192012,1.64719e-07,-5.65088e-11,-42993.7,37.8445], Tmin=(100,'K'), Tmax=(754.613,'K')), NASAPolynomial(coeffs=[13.3842,0.0376067,-2.05004e-05,4.11384e-09,-2.91995e-13,-45008.6,-26.2618], Tmin=(754.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = '[O]OC(F)C([O])F(2090)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {7,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-403.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.434452,'amu*angstrom^2'), symmetry=1, barrier=(9.98892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938913,'amu*angstrom^2'), symmetry=1, barrier=(21.5875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61322,0.0562002,-7.44598e-05,5.23728e-08,-1.4851e-11,-48394.3,22.2333], Tmin=(100,'K'), Tmax=(857.871,'K')), NASAPolynomial(coeffs=[9.58516,0.0190298,-9.4675e-06,1.86665e-09,-1.32691e-13,-49762.1,-15.0047], Tmin=(857.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)[CH]F(502)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-229.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.420111,'amu*angstrom^2'), symmetry=1, barrier=(9.65918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418454,'amu*angstrom^2'), symmetry=1, barrier=(9.62109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.65,'J/mol'), sigma=(5.51737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.71 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66563,0.0577517,-0.000101955,9.1467e-08,-3.1227e-11,-27473.1,20.2392], Tmin=(100,'K'), Tmax=(863.103,'K')), NASAPolynomial(coeffs=[7.65514,0.0169238,-8.28529e-06,1.57108e-09,-1.06541e-13,-28020.1,-4.95473], Tmin=(863.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = '[O]OC(F)C(F)O[O](12440)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {3,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-406.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,180],'cm^-1')),
        HinderedRotor(inertia=(0.732742,'amu*angstrom^2'), symmetry=1, barrier=(16.8472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.731088,'amu*angstrom^2'), symmetry=1, barrier=(16.8091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12952,'amu*angstrom^2'), symmetry=1, barrier=(48.9618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86642,0.0745831,-0.000113957,9.19683e-08,-2.96946e-11,-48738,24.8335], Tmin=(100,'K'), Tmax=(758.633,'K')), NASAPolynomial(coeffs=[10.64,0.0230422,-1.20321e-05,2.38504e-09,-1.68639e-13,-50220.7,-19.617], Tmin=(758.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(ROOJ) + radical(ROOJ)"""),
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
    label = 'F[CH]C(F)OOC=C(F)F(4579)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {7,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-638.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67729,'amu*angstrom^2'), symmetry=1, barrier=(38.5643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67858,'amu*angstrom^2'), symmetry=1, barrier=(38.5938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4336,'amu*angstrom^2'), symmetry=1, barrier=(9.96932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67989,'amu*angstrom^2'), symmetry=1, barrier=(38.624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.85,'J/mol'), sigma=(5.29792,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.55 K, Pc=50.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0936353,0.098574,-0.00014272,1.10721e-07,-3.48592e-11,-76642.2,32.4285], Tmin=(100,'K'), Tmax=(772.882,'K')), NASAPolynomial(coeffs=[12.244,0.0347172,-1.87792e-05,3.80642e-09,-2.73737e-13,-78549.2,-23.914], Tmin=(772.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-638.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'F[CH]OOC=C(F)F(1677)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,D} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {6,D}
8  C u1 p0 c0 {3,S} {5,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-414.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.25081,'amu*angstrom^2'), symmetry=1, barrier=(5.76662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66128,'amu*angstrom^2'), symmetry=1, barrier=(38.1961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998036,0.07318,-0.000119119,1.06928e-07,-3.84225e-11,-49799.7,26.6702], Tmin=(100,'K'), Tmax=(774.496,'K')), NASAPolynomial(coeffs=[8.01183,0.0282627,-1.52892e-05,3.06084e-09,-2.16921e-13,-50625.4,-3.69169], Tmin=(774.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsF1sHO2s)"""),
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
    label = '[O]OC(F)[C](F)OOC=C(F)F(12441)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {9,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {5,S} {9,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-615.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,395,473,707,1436,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.427478,'amu*angstrom^2'), symmetry=1, barrier=(9.82857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427387,'amu*angstrom^2'), symmetry=1, barrier=(9.82647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55538,'amu*angstrom^2'), symmetry=1, barrier=(35.7613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88245,'amu*angstrom^2'), symmetry=1, barrier=(43.2812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89063,'amu*angstrom^2'), symmetry=1, barrier=(43.4693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3215,0.130364,-0.000227018,2.04421e-07,-7.18559e-11,-73870.2,40.7188], Tmin=(100,'K'), Tmax=(806.319,'K')), NASAPolynomial(coeffs=[12.7691,0.0404348,-2.24637e-05,4.48961e-09,-3.15868e-13,-75491.5,-20.1897], Tmin=(806.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]O[C](F)C(F)OOC=C(F)F(12442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {7,S} {9,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-615.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,395,473,707,1436,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.427478,'amu*angstrom^2'), symmetry=1, barrier=(9.82857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427387,'amu*angstrom^2'), symmetry=1, barrier=(9.82647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55538,'amu*angstrom^2'), symmetry=1, barrier=(35.7613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88245,'amu*angstrom^2'), symmetry=1, barrier=(43.2812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89063,'amu*angstrom^2'), symmetry=1, barrier=(43.4693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3215,0.130364,-0.000227018,2.04421e-07,-7.18559e-11,-73870.2,40.7188], Tmin=(100,'K'), Tmax=(806.319,'K')), NASAPolynomial(coeffs=[12.7691,0.0404348,-2.24637e-05,4.48961e-09,-3.15868e-13,-75491.5,-20.1897], Tmin=(806.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)C(F)OO[C]=C(F)F(12443)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {12,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-570.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,562,600,623,1070,1265,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02632,0.124205,-0.000215854,1.98462e-07,-7.14623e-11,-68446.1,41.9838], Tmin=(100,'K'), Tmax=(800.409,'K')), NASAPolynomial(coeffs=[10.6623,0.0439276,-2.44362e-05,4.90157e-09,-3.46086e-13,-69616.9,-7.42988], Tmin=(800.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-570.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = '[O]OC=C(F)OOC=C(F)F(12444)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {5,S} {10,D}
9  C u0 p0 c0 {4,S} {11,D} {12,S}
10 C u0 p0 c0 {6,S} {8,D} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {9,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-384.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,326,540,652,719,1357,2995,3025,975,1000,1300,1375,400,500,1630,1680,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.231982,'amu*angstrom^2'), symmetry=1, barrier=(5.33371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231436,'amu*angstrom^2'), symmetry=1, barrier=(5.32117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233068,'amu*angstrom^2'), symmetry=1, barrier=(5.3587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.320698,0.109711,-0.000199075,1.8971e-07,-6.94016e-11,-46158.2,38.7471], Tmin=(100,'K'), Tmax=(826.237,'K')), NASAPolynomial(coeffs=[7.50733,0.0423239,-2.3199e-05,4.60359e-09,-3.21946e-13,-46445.2,8.56694], Tmin=(826.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-384.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=COOC=C(F)F(6598)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {5,S} {10,D} {13,S}
9  C u0 p0 c0 {4,S} {11,D} {12,S}
10 C u0 p0 c0 {1,S} {6,S} {8,D}
11 C u0 p0 c0 {2,S} {3,S} {9,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-384.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.232311,'amu*angstrom^2'), symmetry=1, barrier=(5.34129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232068,'amu*angstrom^2'), symmetry=1, barrier=(5.3357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232159,'amu*angstrom^2'), symmetry=1, barrier=(5.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.320698,0.109711,-0.000199075,1.8971e-07,-6.94016e-11,-46158.2,38.7471], Tmin=(100,'K'), Tmax=(826.237,'K')), NASAPolynomial(coeffs=[7.50733,0.0423239,-2.3199e-05,4.60359e-09,-3.21946e-13,-46445.2,8.56694], Tmin=(826.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-384.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(F)OOC#CF(12445)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {6,S} {11,T}
11 C u0 p0 c0 {3,S} {10,T}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-321.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,2175,525,239,401,1367,180,304.295],'cm^-1')),
        HinderedRotor(inertia=(1.59821,'amu*angstrom^2'), symmetry=1, barrier=(36.746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60053,'amu*angstrom^2'), symmetry=1, barrier=(36.7994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5994,'amu*angstrom^2'), symmetry=1, barrier=(36.7734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651414,'amu*angstrom^2'), symmetry=1, barrier=(14.9773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59697,'amu*angstrom^2'), symmetry=1, barrier=(36.7174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714414,0.113754,-0.000182737,1.54913e-07,-5.28014e-11,-38457.3,35.1396], Tmin=(100,'K'), Tmax=(717.869,'K')), NASAPolynomial(coeffs=[13.2216,0.0360923,-2.04413e-05,4.1743e-09,-2.99755e-13,-40457.9,-27.4724], Tmin=(717.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(CtCF) + radical(ROOJ)"""),
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
    label = 'O=C(F)C(F)OOC=C(F)F(12446)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {9,D}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-986.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10102,0.0992486,-0.000144343,1.13792e-07,-3.66273e-11,-118497,33.8822], Tmin=(100,'K'), Tmax=(754.759,'K')), NASAPolynomial(coeffs=[11.6871,0.0367726,-2.01743e-05,4.11181e-09,-2.96583e-13,-120277,-19.6718], Tmin=(754.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-986.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'OOC(F)[C](F)OOC=C(F)F(12447)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {6,S} {11,S}
8  O u0 p2 c0 {5,S} {15,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {6,S} {9,S}
11 C u0 p0 c0 {7,S} {12,D} {14,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-767.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,500,795,815,487,638,688,1119,1325,1387,3149,395,473,707,1436,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51303,0.133069,-0.000219138,1.91138e-07,-6.66171e-11,-92143.7,40.6503], Tmin=(100,'K'), Tmax=(754.833,'K')), NASAPolynomial(coeffs=[14.0799,0.0425925,-2.37523e-05,4.80253e-09,-3.42105e-13,-94274.2,-28.7101], Tmin=(754.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OOC(F)C(F)OO[C]=C(F)F(12448)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {12,S}
8  O u0 p2 c0 {6,S} {15,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-722.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,562,600,623,1070,1265,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368046,0.113415,-0.000137473,3.87179e-08,3.83514e-11,-86755.9,39.0506], Tmin=(100,'K'), Tmax=(510.043,'K')), NASAPolynomial(coeffs=[11.9372,0.0461756,-2.57903e-05,5.23245e-09,-3.73974e-13,-88391.8,-15.7619], Tmin=(510.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-722.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C=CJO)"""),
)

species(
    label = 'FOOC(F)[CH]OOC=C(F)F(12449)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {6,S} {11,S}
8  O u0 p2 c0 {4,S} {5,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {6,S} {9,S} {14,S}
11 C u0 p0 c0 {7,S} {12,D} {15,S}
12 C u0 p0 c0 {2,S} {3,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-469.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,500,795,815,277,555,632,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(1.55358,'amu*angstrom^2'), symmetry=1, barrier=(35.7199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55496,'amu*angstrom^2'), symmetry=1, barrier=(35.7515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55448,'amu*angstrom^2'), symmetry=1, barrier=(35.7406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55492,'amu*angstrom^2'), symmetry=1, barrier=(35.7506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5553,'amu*angstrom^2'), symmetry=1, barrier=(35.7595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43441,0.131611,-0.000217434,1.91761e-07,-6.74909e-11,-56278.8,40.7007], Tmin=(100,'K'), Tmax=(765.731,'K')), NASAPolynomial(coeffs=[13.2024,0.0441213,-2.44422e-05,4.92604e-09,-3.50115e-13,-58197,-23.8949], Tmin=(765.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOOC)"""),
)

species(
    label = 'FOO[CH]C(F)OOC=C(F)F(12450)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u0 p0 c0 {6,S} {12,D} {15,S}
12 C u0 p0 c0 {2,S} {3,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-468.453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(1.55124,'amu*angstrom^2'), symmetry=1, barrier=(35.666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55101,'amu*angstrom^2'), symmetry=1, barrier=(35.6608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54832,'amu*angstrom^2'), symmetry=1, barrier=(35.5989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55336,'amu*angstrom^2'), symmetry=1, barrier=(35.7148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55154,'amu*angstrom^2'), symmetry=1, barrier=(35.6729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45712,0.13232,-0.000219739,1.9444e-07,-6.85365e-11,-56157.1,40.7096], Tmin=(100,'K'), Tmax=(769.702,'K')), NASAPolynomial(coeffs=[13.1561,0.0442707,-2.45531e-05,4.94699e-09,-3.51401e-13,-58048.1,-23.636], Tmin=(769.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOO)"""),
)

species(
    label = 'F[C]=COOC(F)C(F)OOF(12451)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {4,S} {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {7,S} {12,D} {15,S}
12 C u1 p0 c0 {3,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-415.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,3010,987.5,1337.5,450,1655,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(1.49258,'amu*angstrom^2'), symmetry=1, barrier=(34.3174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49325,'amu*angstrom^2'), symmetry=1, barrier=(34.3328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49039,'amu*angstrom^2'), symmetry=1, barrier=(34.267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49051,'amu*angstrom^2'), symmetry=1, barrier=(34.2697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4918,'amu*angstrom^2'), symmetry=1, barrier=(34.2994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3571,0.129226,-0.000198964,1.60909e-07,-5.23661e-11,-49799.7,39.9754], Tmin=(100,'K'), Tmax=(750.473,'K')), NASAPolynomial(coeffs=[15.1466,0.0412634,-2.31513e-05,4.73179e-09,-3.40836e-13,-52276.9,-34.9085], Tmin=(750.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    E0 = (-415.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (34.0984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-68.5928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-144.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-351.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-197.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-389.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-273.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-395.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-348.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (28.5306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (26.3888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (80.2735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-409.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-132.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-131.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-280.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-58.8252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-37.6213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-37.6213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (7.56881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (130.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-102.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-102.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (33.2355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-286.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-317.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-247.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-33.8529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-9.7862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (6.88232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['[O]OC(F)C(=O)F(5430)', 'O=CC(F)F(990)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(28.6779,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[O]OC(F)OOC=C(F)F(11532)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.25106e-07,'m^3/(mol*s)'), n=3.34134, Ea=(132.994,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['HO2(11)', 'F[C]C(F)OOC=C(F)F-2(4623)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(47.8902,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['[O]OC(OOC=C(F)F)C(F)F(12433)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['[O]OC(F)C(F)OC(F)(F)C=O(12434)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(92.6483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', '[O]C(F)C(F)OOC=C(F)F(12435)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['O=C[C](F)F(787)', 'FC1OOOC1F(4568)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['HO2(11)', 'FC=C(F)OOC=C(F)F(4618)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(170.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['FC1OO[CH]C(F)(F)OOC1F(12436)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.61579e+10,'s^-1'), n=0.429849, Ea=(48.6953,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['F[C](F)C1OOC(F)C(F)OO1(12437)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;doublebond_intra;radadd_intra_O] + [R8;doublebond_intra;radadd_intra] for rate rule [R8;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]O[CH]C(F)OOC=C(F)F(12438)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]OC(F)[CH]OOC=C(F)F(11540)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC(F)C(F)OOC=[C]F(12439)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C[C](F)F(787)', '[O]OC(F)C([O])F(2090)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(18.7449,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]OC=C(F)F(785)', '[O]OC(F)[CH]F(502)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2CH(73)', '[O]OC(F)C(F)O[O](12440)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62995e+14,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', 'F[CH]C(F)OOC=C(F)F(4579)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]O[CH]F(287)', 'F[CH]OOC=C(F)F(1677)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]OC(F)[C](F)OOC=C(F)F(12441)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]O[C](F)C(F)OOC=C(F)F(12442)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]OC(F)C(F)OO[C]=C(F)F(12443)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F2(78)', '[O]OC=COOC=C(F)F(11542)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]OC=C(F)OOC=C(F)F(12444)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(196.849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', '[O]OC(F)=COOC=C(F)F(6598)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(196.849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', '[O]OC(F)C(F)OOC#CF(12445)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(269.181,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    products = ['OH(7)', 'O=C(F)C(F)OOC=C(F)F(12446)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OOC(F)[C](F)OOC=C(F)F(12447)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;O_H_out] for rate rule [R4H_SSS;C_rad_out_noH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OOC(F)C(F)OO[C]=C(F)F(12448)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.81284e+06,'s^-1'), n=1.31318, Ea=(108.582,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R7H_OOCs4;Y_rad_out;XH_out] for rate rule [R7H_OOCs4;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.1622776601683795
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['FOOC(F)[CH]OOC=C(F)F(12449)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(69.3539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction30',
    reactants = ['FOO[CH]C(F)OOC=C(F)F(12450)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.4141,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]=COOC(F)C(F)OOF(12451)'],
    products = ['[O]OC(F)C(F)OOC=C(F)F(12368)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(56.2044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #4021',
    isomers = [
        '[O]OC(F)C(F)OOC=C(F)F(12368)',
    ],
    reactants = [
        ('[O]OC(F)C(=O)F(5430)', 'O=CC(F)F(990)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4021',
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

