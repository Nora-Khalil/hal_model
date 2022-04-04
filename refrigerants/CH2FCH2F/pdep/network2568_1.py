species(
    label = 'OC(OC(F)C(F)F)C(F)F(6763)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {16,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
11 C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {7,S}
"""),
    E0 = (-1500.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,176,294,457,589,547,707,1086,1160,1127,1157,1302,1442,1365,1447,3042,3152,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.6379,0.110258,-0.000143295,9.77235e-08,-2.69615e-11,-180261,34.9077], Tmin=(100,'K'), Tmax=(878.513,'K')), NASAPolynomial(coeffs=[15.2483,0.037924,-1.97858e-05,3.99483e-09,-2.88186e-13,-183052,-39.6757], Tmin=(878.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1500.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F))"""),
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
    label = 'F[C]C(F)F-2(1828)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p1 c0 {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-295.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([257,409,1143,1361,2944,617,898,1187,180,1242.56,1243.85],'cm^-1')),
        HinderedRotor(inertia=(0.0989364,'amu*angstrom^2'), symmetry=1, barrier=(2.27474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89361,0.016532,-3.68777e-06,-7.04441e-09,3.72589e-12,-35574,10.6199], Tmin=(10,'K'), Tmax=(1003.12,'K')), NASAPolynomial(coeffs=[6.83938,0.0116737,-6.72299e-06,1.81811e-09,-1.88858e-13,-36511.5,-5.32834], Tmin=(1003.12,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-295.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(O)C(F)F(6564)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {9,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-853.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,180],'cm^-1')),
        HinderedRotor(inertia=(1.01563,'amu*angstrom^2'), symmetry=1, barrier=(23.3513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01598,'amu*angstrom^2'), symmetry=1, barrier=(23.3593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01588,'amu*angstrom^2'), symmetry=1, barrier=(23.3571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0488,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.80855,0.0119068,0.000150141,-3.53338e-07,2.36034e-10,-102606,10.8049], Tmin=(10,'K'), Tmax=(531.266,'K')), NASAPolynomial(coeffs=[6.21374,0.0338789,-2.50634e-05,8.52892e-09,-1.07724e-12,-103427,-4.60097], Tmin=(531.266,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-853.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(O)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(OC(F)F)C(F)F(6011)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1293.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,519,593,830,1115,1166,1389,1413,3114,245.904,245.912,245.922,245.929],'cm^-1')),
        HinderedRotor(inertia=(0.161786,'amu*angstrom^2'), symmetry=1, barrier=(6.94289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00278739,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161763,'amu*angstrom^2'), symmetry=1, barrier=(6.94252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386783,'amu*angstrom^2'), symmetry=1, barrier=(16.6022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.056,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3552.3,'J/mol'), sigma=(5.67295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=554.86 K, Pc=44.15 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320095,0.0879944,-0.000128918,1.03092e-07,-3.3429e-11,-155407,29.8161], Tmin=(100,'K'), Tmax=(751.948,'K')), NASAPolynomial(coeffs=[10.9328,0.0315383,-1.62958e-05,3.24007e-09,-2.30198e-13,-157003,-18.3582], Tmin=(751.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1293.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CsFFHO)"""),
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
    label = 'OCOC(F)C(F)F(6741)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {8,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1043.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,2750,2850,1437.5,1250,1305,750,350,306.344,306.398,306.4,306.455],'cm^-1')),
        HinderedRotor(inertia=(0.30938,'amu*angstrom^2'), symmetry=1, barrier=(20.6077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309201,'amu*angstrom^2'), symmetry=1, barrier=(20.6073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30896,'amu*angstrom^2'), symmetry=1, barrier=(20.6081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30947,'amu*angstrom^2'), symmetry=1, barrier=(20.6084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3532.06,'J/mol'), sigma=(5.75458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.70 K, Pc=42.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140204,0.0699582,-4.37105e-05,-9.64689e-09,1.2892e-11,-125405,27.2818], Tmin=(100,'K'), Tmax=(954.428,'K')), NASAPolynomial(coeffs=[23.6542,0.00728002,-1.57639e-06,2.98434e-10,-2.70369e-14,-131527,-93.6214], Tmin=(954.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1043.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-OsOsHH)"""),
)

species(
    label = 'OC([C]F)OC(F)C(F)F(8740)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {8,S} {14,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {13,S}
10 C u0 p1 c0 {4,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-912.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,617,898,1187,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0960197,0.0972215,-0.000130381,9.15181e-08,-2.58747e-11,-109560,33.4805], Tmin=(100,'K'), Tmax=(859.885,'K')), NASAPolynomial(coeffs=[13.9178,0.0320316,-1.66615e-05,3.35048e-09,-2.40859e-13,-111970,-32.0124], Tmin=(859.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-912.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + group(Cs-CsOsOsH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)OC(F)C(F)F(6014)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1279.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,509,613,660,1171,1360,1414,3084,283.186,283.26,283.33,283.416],'cm^-1')),
        HinderedRotor(inertia=(0.207118,'amu*angstrom^2'), symmetry=1, barrier=(11.7955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.392571,'amu*angstrom^2'), symmetry=1, barrier=(22.3636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206958,'amu*angstrom^2'), symmetry=1, barrier=(11.7935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.584788,'amu*angstrom^2'), symmetry=1, barrier=(33.2978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.056,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3552.3,'J/mol'), sigma=(5.67295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=554.86 K, Pc=44.15 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340597,0.0854905,-0.000106295,6.8069e-08,-1.75069e-11,-153794,27.1], Tmin=(100,'K'), Tmax=(942.801,'K')), NASAPolynomial(coeffs=[14.0114,0.0274884,-1.40119e-05,2.81287e-09,-2.02719e-13,-156372,-38.0484], Tmin=(942.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1279.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFHOO)"""),
)

species(
    label = 'OC(OCF)C(F)F(6737)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {13,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1060.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,548,1085,1183,1302,1466,1520,3060,3119,382.596,382.597,382.597,382.597],'cm^-1')),
        HinderedRotor(inertia=(0.0745367,'amu*angstrom^2'), symmetry=1, barrier=(7.74246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0745367,'amu*angstrom^2'), symmetry=1, barrier=(7.74246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0745368,'amu*angstrom^2'), symmetry=1, barrier=(7.74247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195952,'amu*angstrom^2'), symmetry=1, barrier=(20.3545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3532.06,'J/mol'), sigma=(5.75458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=551.70 K, Pc=42.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901024,0.0743144,-9.3328e-05,6.56493e-08,-1.9163e-11,-127471,26.7874], Tmin=(100,'K'), Tmax=(824.436,'K')), NASAPolynomial(coeffs=[9.64737,0.0318785,-1.61184e-05,3.21446e-09,-2.30268e-13,-128913,-13.7201], Tmin=(824.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1060.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CsFHHO)"""),
)

species(
    label = 'OC(OC(F)[C]F)C(F)F(8741)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {10,S} {13,S}
10 C u0 p1 c0 {4,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-932.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,2750,2850,1437.5,1250,1305,750,350,617,898,1187,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218765,0.100956,-0.000144469,1.09536e-07,-3.34971e-11,-111949,33.4509], Tmin=(100,'K'), Tmax=(796.781,'K')), NASAPolynomial(coeffs=[13.1797,0.0336947,-1.78484e-05,3.59491e-09,-2.57796e-13,-114085,-28.1456], Tmin=(796.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-932.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CsCFHO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)C(F)OC(F)C(F)F(8742)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {10,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {9,S} {15,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {7,S}
"""),
    E0 = (-1486.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,205,261,317,423,493,563,514,600,686,1105,1152,1199,1293,1365,1437,1381,1422,1463,3031,3097,3163,235,523,627,1123,1142,1372,1406,3097,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.48565,0.103398,-0.000115827,6.57306e-08,-1.50342e-11,-178579,32.8787], Tmin=(100,'K'), Tmax=(1051.1,'K')), NASAPolynomial(coeffs=[17.5429,0.0347884,-1.79153e-05,3.62849e-09,-2.6325e-13,-182369,-54.9969], Tmin=(1051.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1486.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = 'OC=C(F)F(1567)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-502.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.552091,'amu*angstrom^2'), symmetry=1, barrier=(12.6937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3292.1,'J/mol'), sigma=(5.16209,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.22 K, Pc=54.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85967,0.00956102,9.95844e-05,-2.79091e-07,2.18422e-10,-60486.4,9.26794], Tmin=(10,'K'), Tmax=(463.499,'K')), NASAPolynomial(coeffs=[6.00669,0.0174593,-1.1501e-05,3.70008e-09,-4.58213e-13,-60969.3,-2.50148], Tmin=(463.499,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-502.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)C(F)F(1614)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {9,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-880.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,180],'cm^-1')),
        HinderedRotor(inertia=(1.07214,'amu*angstrom^2'), symmetry=1, barrier=(24.6507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07211,'amu*angstrom^2'), symmetry=1, barrier=(24.6499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84462,0.00955472,0.000130159,-2.95984e-07,1.92886e-10,-105845,11.1136], Tmin=(10,'K'), Tmax=(538.658,'K')), NASAPolynomial(coeffs=[5.31308,0.0316785,-2.34231e-05,7.92737e-09,-9.9477e-13,-106482,0.489879], Tmin=(538.658,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-880.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC(F)C(F)F(3101)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-739.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.797804,'amu*angstrom^2'), symmetry=1, barrier=(18.3431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67499,'amu*angstrom^2'), symmetry=1, barrier=(38.5113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62455,0.0310837,4.25757e-05,-1.64434e-07,1.26508e-10,-88891.9,11.8283], Tmin=(10,'K'), Tmax=(506.268,'K')), NASAPolynomial(coeffs=[6.43956,0.0309358,-2.24458e-05,7.38648e-09,-9.0465e-13,-89460,-2.63266], Tmin=(506.268,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-739.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""FOC(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC=CF(455)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-317.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(1.37732,'amu*angstrom^2'), symmetry=1, barrier=(31.6673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95033,0.00268027,6.73468e-05,-1.15374e-07,5.90231e-11,-38204.1,8.04413], Tmin=(10,'K'), Tmax=(637.71,'K')), NASAPolynomial(coeffs=[1.89968,0.0278733,-2.09144e-05,7.21532e-09,-9.21467e-13,-38193.3,15.049], Tmin=(637.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-317.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787929,2.90875e-06,-1.47516e-09,2.12833e-13,-30281.6,-0.311362], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81028e-07,1.07243e-10,-5.79385e-15,-29958.6,5.9104], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC(F)=COC(F)C(F)F(6332)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1159.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3084,'amu*angstrom^2'), symmetry=1, barrier=(30.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30931,'amu*angstrom^2'), symmetry=1, barrier=(30.1035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30888,'amu*angstrom^2'), symmetry=1, barrier=(30.0938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.310191,0.0890815,-9.74557e-05,5.08771e-08,-1.03013e-11,-139243,29.0591], Tmin=(100,'K'), Tmax=(1209.55,'K')), NASAPolynomial(coeffs=[21.4913,0.0169833,-8.04402e-06,1.59596e-09,-1.15352e-13,-144517,-80.2681], Tmin=(1209.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1159.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'OF(140)',
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
    label = 'FC=COC(F)C(F)F(8743)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {4,S} {8,D} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-966.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29264,'amu*angstrom^2'), symmetry=1, barrier=(29.7204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29164,'amu*angstrom^2'), symmetry=1, barrier=(29.6973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.288,'amu*angstrom^2'), symmetry=1, barrier=(29.6136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287054,0.0825013,-8.42936e-05,4.07962e-08,-7.59226e-12,-116098,27.9049], Tmin=(100,'K'), Tmax=(1322.22,'K')), NASAPolynomial(coeffs=[22.6955,0.0129739,-5.41745e-06,1.02648e-09,-7.27318e-14,-122176,-89.3917], Tmin=(1322.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-966.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'CHFCHF[Z](59)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(OF)C(F)F(4136)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-719.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.958732,'amu*angstrom^2'), symmetry=1, barrier=(22.0431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958327,'amu*angstrom^2'), symmetry=1, barrier=(22.0338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960571,'amu*angstrom^2'), symmetry=1, barrier=(22.0854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69892,0.0210936,0.000171449,-5.05684e-07,4.0418e-10,-86478.5,11.6887], Tmin=(10,'K'), Tmax=(457.754,'K')), NASAPolynomial(coeffs=[7.73788,0.0334981,-2.54998e-05,8.78503e-09,-1.11573e-12,-87348,-10.0993], Tmin=(457.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-719.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(OF)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(O[CH]C(F)F)C(F)F(8744)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {7,S} {15,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
10 C u1 p0 c0 {5,S} {9,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1108.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.849673,0.115893,-0.000185632,1.5738e-07,-5.29544e-11,-133128,35.6172], Tmin=(100,'K'), Tmax=(787.184,'K')), NASAPolynomial(coeffs=[13.441,0.0357239,-1.84766e-05,3.62862e-09,-2.54144e-13,-135143,-28.421], Tmin=(787.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1108.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H))"""),
)

species(
    label = 'OC(OC(F)[CH]F)C(F)F(8745)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {15,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
10 C u1 p0 c0 {4,S} {8,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1085.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,487,638,688,1119,1325,1387,3149,235,523,627,1123,1142,1372,1406,3097,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.562525,0.108723,-0.000155184,1.16887e-07,-3.53991e-11,-130360,34.097], Tmin=(100,'K'), Tmax=(805.48,'K')), NASAPolynomial(coeffs=[14.2248,0.0352869,-1.8421e-05,3.68875e-09,-2.63679e-13,-132742,-34.0437], Tmin=(805.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1085.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'OC([CH]F)OC(F)C(F)F(8746)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {15,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {13,S}
10 C u1 p0 c0 {4,S} {7,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (-1080.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.919897,0.118265,-0.000192341,1.65903e-07,-5.6774e-11,-129756,34.9236], Tmin=(100,'K'), Tmax=(785.88,'K')), NASAPolynomial(coeffs=[13.0495,0.0375156,-1.98037e-05,3.91871e-09,-2.75513e-13,-131654,-27.2096], Tmin=(785.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1080.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
)

species(
    label = '[O]C(F)C(F)F(1662)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-651.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.650368,'amu*angstrom^2'), symmetry=1, barrier=(14.9532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79706,0.0159808,9.87453e-05,-3.10799e-07,2.69371e-10,-78322,11.5156], Tmin=(10,'K'), Tmax=(408.816,'K')), NASAPolynomial(coeffs=[4.76899,0.0273062,-1.92554e-05,6.29186e-09,-7.71788e-13,-78575.6,5.56657], Tmin=(408.816,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-651.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[O]C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[CH]C(F)F(1906)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {8,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-481.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,309.104],'cm^-1')),
        HinderedRotor(inertia=(0.27531,'amu*angstrom^2'), symmetry=1, barrier=(18.6693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471239,'amu*angstrom^2'), symmetry=1, barrier=(31.9681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89233,0.00668522,0.000104616,-2.31512e-07,1.50155e-10,-57892.5,11.6044], Tmin=(10,'K'), Tmax=(527.066,'K')), NASAPolynomial(coeffs=[3.95879,0.0284142,-2.04983e-05,6.77504e-09,-8.32891e-13,-58208.3,8.39674], Tmin=(527.066,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-481.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""O[CH]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF2-CHF(67)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-481.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.191999,'amu*angstrom^2'), symmetry=1, barrier=(4.41444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0324,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03194,0.0198962,-1.06494e-05,1.96763e-09,3.06708e-14,-57935.2,10.3725], Tmin=(10,'K'), Tmax=(1633.48,'K')), NASAPolynomial(coeffs=[15.0609,-0.000860555,2.67159e-06,-1.12649e-09,1.457e-13,-62372.2,-50.8003], Tmin=(1633.48,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-481.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[CH]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)C(F)F(4954)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-620.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,180,1741.63],'cm^-1')),
        HinderedRotor(inertia=(0.565755,'amu*angstrom^2'), symmetry=1, barrier=(13.0078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11337,'amu*angstrom^2'), symmetry=1, barrier=(25.5985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63997,0.0313927,6.9387e-07,-4.94196e-08,3.71681e-11,-74675.6,12.0463], Tmin=(10,'K'), Tmax=(587.838,'K')), NASAPolynomial(coeffs=[6.02864,0.0268837,-1.77702e-05,5.50932e-09,-6.475e-13,-75159.3,0.065464], Tmin=(587.838,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-620.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(O)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC(F)[CH]OC(F)C(F)F(6727)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {6,S} {9,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1113.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.358666,0.103835,-0.000149726,1.14204e-07,-3.5019e-11,-133816,31.3506], Tmin=(100,'K'), Tmax=(795.982,'K')), NASAPolynomial(coeffs=[13.5965,0.0337045,-1.7562e-05,3.50696e-09,-2.50081e-13,-136037,-32.7901], Tmin=(795.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1113.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H))"""),
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
    label = '[O]C(OC(F)C(F)F)C(F)F(8747)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
11 C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1274.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,176,294,457,589,547,707,1086,1160,1127,1157,1302,1442,1365,1447,3042,3152,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0696513,0.103396,-0.00012981,6.26268e-08,6.74048e-12,-153142,33.8167], Tmin=(100,'K'), Tmax=(540.044,'K')), NASAPolynomial(coeffs=[10.848,0.0438633,-2.37043e-05,4.78512e-09,-3.42382e-13,-154632,-15.0081], Tmin=(540.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1274.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(CCOJ)"""),
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
    label = 'O[CH]OC(F)C(F)F(6500)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-854.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,3025,407.5,1350,352.5,390.947,390.948,390.948,390.948],'cm^-1')),
        HinderedRotor(inertia=(0.151547,'amu*angstrom^2'), symmetry=1, barrier=(16.4366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985131,'amu*angstrom^2'), symmetry=1, barrier=(10.6846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270556,'amu*angstrom^2'), symmetry=1, barrier=(29.3441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985133,'amu*angstrom^2'), symmetry=1, barrier=(10.6846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430213,0.0821221,-0.000108733,7.27826e-08,-1.92388e-11,-102705,26.0442], Tmin=(100,'K'), Tmax=(925.998,'K')), NASAPolynomial(coeffs=[14.4683,0.021483,-1.05067e-05,2.06576e-09,-1.46989e-13,-105305,-40.6023], Tmin=(925.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-854.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-OsOsHH) + radical(OCJO)"""),
)

species(
    label = 'O[C](OC(F)C(F)F)C(F)F(8748)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {11,S} {15,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {14,S}
11 C u1 p0 c0 {6,S} {7,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1294.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,522,611,926,1093,1137,1374,1416,3112,360,370,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702957,0.112242,-0.000158339,1.16312e-07,-3.43039e-11,-155574,35.8261], Tmin=(100,'K'), Tmax=(826.415,'K')), NASAPolynomial(coeffs=[15.1699,0.0354148,-1.8895e-05,3.82362e-09,-2.75228e-13,-158197,-37.7251], Tmin=(826.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1294.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Cs_P)"""),
)

species(
    label = 'OC(O[CH]F)C(F)F(6567)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-862.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,580,1155,1237,1373,3147,258.398,258.984,259.005,259.364],'cm^-1')),
        HinderedRotor(inertia=(0.130298,'amu*angstrom^2'), symmetry=1, barrier=(6.21547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130683,'amu*angstrom^2'), symmetry=1, barrier=(6.22549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00251107,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386419,'amu*angstrom^2'), symmetry=1, barrier=(18.3631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527514,0.0835498,-0.000132835,1.14782e-07,-3.96135e-11,-103562,28.3633], Tmin=(100,'K'), Tmax=(782.578,'K')), NASAPolynomial(coeffs=[9.76922,0.029237,-1.51687e-05,2.9909e-09,-2.10207e-13,-104792,-12.5726], Tmin=(782.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-862.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CsFHHO) + radical(Csj(F1s)(O2s-Cs)(H))"""),
)

species(
    label = 'OC(O[C](F)C(F)F)C(F)F(8749)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {8,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1296.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,522,611,926,1093,1137,1374,1416,3112,395,473,707,1436,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.731168,0.113704,-0.000174427,1.44554e-07,-4.84505e-11,-155791,37.5641], Tmin=(100,'K'), Tmax=(728.497,'K')), NASAPolynomial(coeffs=[12.8318,0.039236,-2.11012e-05,4.24777e-09,-3.03298e-13,-157767,-23.5736], Tmin=(728.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1296.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sF1sH)(F1s)(O2s-Cs))"""),
)

species(
    label = 'OC(OC(F)C(F)F)[C](F)F(8750)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1292.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.59173,0.112785,-0.000163026,1.16294e-07,-2.85211e-11,-155295,36.2775], Tmin=(100,'K'), Tmax=(611.105,'K')), NASAPolynomial(coeffs=[13.1937,0.0386689,-2.06614e-05,4.14198e-09,-2.94899e-13,-157281,-25.9023], Tmin=(611.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1292.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(F1s))"""),
)

species(
    label = 'OC(OC(F)[C](F)F)C(F)F(8751)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {6,S} {11,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
11 C u1 p0 c0 {4,S} {5,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
"""),
    E0 = (-1292.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,487,638,688,1119,1325,1387,3149,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.708536,0.112903,-0.000162032,1.22202e-07,-3.71273e-11,-155258,36.6101], Tmin=(100,'K'), Tmax=(802.027,'K')), NASAPolynomial(coeffs=[14.5309,0.0368953,-1.98712e-05,4.02977e-09,-2.90126e-13,-157703,-33.5487], Tmin=(802.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1292.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCsFFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = 'OC(=CF)OC(F)C(F)F(8752)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1123.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,194,682,905,1196,1383,3221,188.554,189.007,189.198,189.793],'cm^-1')),
        HinderedRotor(inertia=(0.870227,'amu*angstrom^2'), symmetry=1, barrier=(22.322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862937,'amu*angstrom^2'), symmetry=1, barrier=(22.3268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885781,'amu*angstrom^2'), symmetry=1, barrier=(22.3253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885857,'amu*angstrom^2'), symmetry=1, barrier=(22.3317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.614868,0.101654,-0.000124999,7.46338e-08,-1.73713e-11,-134996,28.8672], Tmin=(100,'K'), Tmax=(1053.6,'K')), NASAPolynomial(coeffs=[20.6079,0.0210807,-1.02853e-05,2.04776e-09,-1.47761e-13,-139468,-74.6282], Tmin=(1053.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1123.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'OC(OC=CF)C(F)F(8753)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {8,D} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (-961.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,319.515,319.922,320.049,320.087],'cm^-1')),
        HinderedRotor(inertia=(0.256836,'amu*angstrom^2'), symmetry=1, barrier=(18.6704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257564,'amu*angstrom^2'), symmetry=1, barrier=(18.6543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256378,'amu*angstrom^2'), symmetry=1, barrier=(18.6627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257523,'amu*angstrom^2'), symmetry=1, barrier=(18.6624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.215351,0.0877637,-9.55629e-05,5.04849e-08,-1.03856e-11,-115435,30.1623], Tmin=(100,'K'), Tmax=(1189.24,'K')), NASAPolynomial(coeffs=[20.325,0.018677,-8.42425e-06,1.63705e-09,-1.17027e-13,-120320,-72.4933], Tmin=(1189.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-961.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'OC(OC(F)=CF)C(F)F(8754)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1148.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180,585.288,1510.88],'cm^-1')),
        HinderedRotor(inertia=(0.703681,'amu*angstrom^2'), symmetry=1, barrier=(16.179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701451,'amu*angstrom^2'), symmetry=1, barrier=(16.1277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701832,'amu*angstrom^2'), symmetry=1, barrier=(16.1365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702332,'amu*angstrom^2'), symmetry=1, barrier=(16.148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.258628,0.0986802,-0.000126945,8.23522e-08,-2.12108e-11,-137942,31.0527], Tmin=(100,'K'), Tmax=(946.474,'K')), NASAPolynomial(coeffs=[16.5985,0.0274375,-1.4036e-05,2.82167e-09,-2.03507e-13,-141133,-49.3457], Tmin=(946.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1148.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OC(OC=C(F)F)C(F)F(5024)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1153.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,350.341,350.341,350.341,350.341],'cm^-1')),
        HinderedRotor(inertia=(0.213103,'amu*angstrom^2'), symmetry=1, barrier=(18.5609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213102,'amu*angstrom^2'), symmetry=1, barrier=(18.5609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213102,'amu*angstrom^2'), symmetry=1, barrier=(18.5609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213103,'amu*angstrom^2'), symmetry=1, barrier=(18.5609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.360787,0.0956352,-0.000112598,6.4793e-08,-1.46029e-11,-138574,31.7661], Tmin=(100,'K'), Tmax=(1084.26,'K')), NASAPolynomial(coeffs=[19.6087,0.021964,-1.06776e-05,2.12578e-09,-1.53407e-13,-142904,-66.1904], Tmin=(1084.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1153.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    E0 = (-719.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-498.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-413.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-517.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-559.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-383.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-529.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-575.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-579.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-590.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-187.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-571.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-192.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-571.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-159.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-427.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-404.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-399.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-524.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-494.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-477.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-454.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-503.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-475.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-511.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-477.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-472.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-472.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-615.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-355.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-628.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-595.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    products = ['HF(38)', 'O=CC(F)F(990)', 'O=CC(F)F(990)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.29418e+10,'s^-1'), n=0.745914, Ea=(173.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]C(F)F-2(1828)', 'OC(O)C(F)F(6564)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.5279e-12,'m^3/(mol*s)'), n=5.02686, Ea=(42.8703,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'OC(OC(F)F)C(F)F(6011)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.12553e-07,'m^3/(mol*s)'), n=3.34134, Ea=(133.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(43)', 'OCOC(F)C(F)F(6741)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33604e-10,'m^3/(mol*s)'), n=4.72997, Ea=(122.032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HF(38)', 'OC([C]F)OC(F)C(F)F(8740)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.4487,'m^3/(mol*s)'), n=1.52711, Ea=(26.4489,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.02090651070256912, var=1.7269912187795557, Tref=1000.0, N=2, data_mean=0.0, correlation='HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R',), comment="""Estimated from node HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHF(40)', 'OC(F)OC(F)C(F)F(6014)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96165e-06,'m^3/(mol*s)'), n=3.30609, Ea=(149.698,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s',), comment="""Estimated from node CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CF2(43)', 'OC(OCF)C(F)F(6737)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.33604e-10,'m^3/(mol*s)'), n=4.72997, Ea=(127.616,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HF(38)', 'OC(OC(F)[C]F)C(F)F(8741)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.4487,'m^3/(mol*s)'), n=1.52711, Ea=(30.3063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.02090651070256912, var=1.7269912187795557, Tref=1000.0, N=2, data_mean=0.0, correlation='HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R',), comment="""Estimated from node HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R_Ext-3Br1sCCl1sHI1s-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OC(F)C(F)OC(F)C(F)F(8742)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC=C(F)F(1567)', 'OC(F)C(F)F(1614)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.78846e-05,'m^3/(mol*s)'), n=2.92667, Ea=(185.212,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OCs] for rate rule [Cd/disub_Cd/monosub;H_OCs]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['FOC(F)C(F)F(3101)', 'OC=CF(455)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.119e-05,'m^3/(mol*s)'), n=3.25, Ea=(261.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd/monosub_Cd/monosub;H_OR]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H2O(3)', 'FC(F)=COC(F)C(F)F(6332)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OF(140)', 'FC=COC(F)C(F)F(8743)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.119e-05,'m^3/(mol*s)'), n=3.25, Ea=(261.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd/monosub_Cd/monosub;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHFCF2(55)', 'OC(O)C(F)F(6564)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.57693e-05,'m^3/(mol*s)'), n=2.92667, Ea=(185.212,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OCs] for rate rule [Cd/disub_Cd/monosub;H_OCs]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHFCHF[Z](59)', 'OC(OF)C(F)F(4136)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.238e-05,'m^3/(mol*s)'), n=3.25, Ea=(261.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd/monosub_Cd/monosub;H_OR]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'OC(O[CH]C(F)F)C(F)F(8744)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'OC(OC(F)[CH]F)C(F)F(8745)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'OC([CH]F)OC(F)C(F)F(8746)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)C(F)F(1662)', 'O[CH]C(F)F(1906)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CHF2-CHF(67)', '[O]C(O)C(F)F(4954)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['OH(7)', 'FC(F)[CH]OC(F)C(F)F(6727)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[O]C(OC(F)C(F)F)C(F)F(8747)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CHF2(82)', 'O[CH]OC(F)C(F)F(6500)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -13.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', 'O[C](OC(F)C(F)F)C(F)F(8748)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CHF2(82)', 'OC(O[CH]F)C(F)F(6567)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(5)', 'OC(O[C](F)C(F)F)C(F)F(8749)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.06317e+29,'m^3/(mol*s)'), n=-7.6173, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.30158007722040336, var=26.748226087879157, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R
Ea raised from -1.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(5)', 'OC(OC(F)C(F)F)[C](F)F(8750)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl
Ea raised from -2.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(5)', 'OC(OC(F)[C](F)F)C(F)F(8751)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.06317e+29,'m^3/(mol*s)'), n=-7.6173, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.30158007722040336, var=26.748226087879157, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R
Ea raised from -2.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'OC(=CF)OC(F)C(F)F(8752)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(181.629,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F2(78)', 'OC(OC=CF)C(F)F(8753)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(6.30349,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction31',
    reactants = ['HF(38)', 'OC(OC(F)=CF)C(F)F(8754)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(192.585,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HF(38)', 'OC(OC=C(F)F)C(F)F(5024)'],
    products = ['OC(OC(F)C(F)F)C(F)F(6763)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(231.302,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #2568',
    isomers = [
        'OC(OC(F)C(F)F)C(F)F(6763)',
    ],
    reactants = [
        ('HF(38)', 'O=CC(F)F(990)', 'O=CC(F)F(990)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2568',
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

