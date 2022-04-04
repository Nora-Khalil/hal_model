species(
    label = 'CC(F)(F)C(OF)C(=O)F(8310)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-960.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,180,621.603,622.25],'cm^-1')),
        HinderedRotor(inertia=(0.119002,'amu*angstrom^2'), symmetry=1, barrier=(2.73609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698776,'amu*angstrom^2'), symmetry=1, barrier=(16.0662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698737,'amu*angstrom^2'), symmetry=1, barrier=(16.0653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698955,'amu*angstrom^2'), symmetry=1, barrier=(16.0704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169846,0.0895525,-0.000105733,6.41873e-08,-1.574e-11,-115414,30.9673], Tmin=(100,'K'), Tmax=(983.278,'K')), NASAPolynomial(coeffs=[14.6671,0.0305773,-1.57656e-05,3.18895e-09,-2.31014e-13,-118265,-38.729], Tmin=(983.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-960.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(COCsFO)"""),
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
    label = 'O=CC(=O)F(4234)',
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
    label = 'CH2CF2(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CC(F)(F)C(F)OF(5936)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-799.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.470475,'amu*angstrom^2'), symmetry=1, barrier=(10.8171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29368,'amu*angstrom^2'), symmetry=1, barrier=(29.7443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407965,'amu*angstrom^2'), symmetry=1, barrier=(9.37991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2934.34,'J/mol'), sigma=(5.23047,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.34 K, Pc=46.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57271,0.0349612,0.000139035,-5.05108e-07,4.64945e-10,-96180.9,12.2953], Tmin=(10,'K'), Tmax=(401.734,'K')), NASAPolynomial(coeffs=[6.60006,0.0425801,-3.04076e-05,1.00536e-08,-1.2449e-12,-96728.9,-3.3416], Tmin=(401.734,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-799.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""CC(F)(F)C(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CC(OF)C(=O)F(8360)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-492.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,180,934.324],'cm^-1')),
        HinderedRotor(inertia=(0.0690357,'amu*angstrom^2'), symmetry=1, barrier=(1.58727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842068,'amu*angstrom^2'), symmetry=1, barrier=(19.3608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842057,'amu*angstrom^2'), symmetry=1, barrier=(19.3606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25491,0.0522627,-4.43752e-05,1.82847e-08,-2.95877e-12,-59151.3,25.1678], Tmin=(100,'K'), Tmax=(1488.61,'K')), NASAPolynomial(coeffs=[15.4003,0.014253,-6.07466e-06,1.13194e-09,-7.8108e-14,-63362.7,-48.7031], Tmin=(1488.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-492.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(COCsFO)"""),
)

species(
    label = 'CH3F(39)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-246.136,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.0219,'amu')),
        NonlinearRotor(inertia=([3.21737,19.5168,19.5168],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([1097.73,1206.94,1206.96,1493.65,1493.66,1509.77,3011.42,3086.11,3086.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (34.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05191,-0.00340134,1.7341e-05,2.49456e-09,-1.5179e-11,-29605.1,3.91759], Tmin=(10,'K'), Tmax=(534.708,'K')), NASAPolynomial(coeffs=[0.487074,0.0143029,-7.17972e-06,1.7168e-09,-1.57925e-13,-29095.7,20.0824], Tmin=(534.708,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-246.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C([C]F)OF(8132)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u0 p1 c0 {2,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-315.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,617,898,1187,294.231,294.235],'cm^-1')),
        HinderedRotor(inertia=(0.115971,'amu*angstrom^2'), symmetry=1, barrier=(7.12467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115967,'amu*angstrom^2'), symmetry=1, barrier=(7.12453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767789,'amu*angstrom^2'), symmetry=1, barrier=(47.1685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29169,0.0644044,-9.37974e-05,7.10008e-08,-2.14837e-11,-37893.1,25.7046], Tmin=(100,'K'), Tmax=(807.806,'K')), NASAPolynomial(coeffs=[10.3781,0.0194121,-1.02538e-05,2.05521e-09,-1.46802e-13,-39361.1,-16.1931], Tmin=(807.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(COCsFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'O=C(F)C(F)OF(5789)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {5,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-655.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,180,243.882],'cm^-1')),
        HinderedRotor(inertia=(0.69142,'amu*angstrom^2'), symmetry=1, barrier=(15.8971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80703,'amu*angstrom^2'), symmetry=1, barrier=(41.5471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62699,0.0398861,-4.04284e-05,2.02407e-08,-4.0651e-12,-78828.3,11.5278], Tmin=(10,'K'), Tmax=(1137.44,'K')), NASAPolynomial(coeffs=[9.76274,0.0183084,-1.19725e-05,3.56222e-09,-3.99254e-13,-80224.1,-18.8637], Tmin=(1137.44,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-655.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(F)C(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C[C]F(325)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (80.9091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1049.2],'cm^-1')),
        HinderedRotor(inertia=(0.0264149,'amu*angstrom^2'), symmetry=1, barrier=(3.87767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50524,0.00692172,1.56037e-05,-2.26104e-08,8.49049e-12,9752.45,9.13911], Tmin=(100,'K'), Tmax=(981.317,'K')), NASAPolynomial(coeffs=[5.27397,0.00958137,-3.54752e-06,6.48805e-10,-4.59759e-14,8930.11,-1.78188], Tmin=(981.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C(F)C(OF)C(F)F(8361)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-903.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,486,617,768,1157,1926,180,496.977,1375.25],'cm^-1')),
        HinderedRotor(inertia=(0.814739,'amu*angstrom^2'), symmetry=1, barrier=(18.7325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.814295,'amu*angstrom^2'), symmetry=1, barrier=(18.7223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813587,'amu*angstrom^2'), symmetry=1, barrier=(18.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895458,0.0702475,-7.95477e-05,4.44717e-08,-9.89387e-12,-108499,28.1778], Tmin=(100,'K'), Tmax=(1086.37,'K')), NASAPolynomial(coeffs=[14.4561,0.0203169,-1.06056e-05,2.16406e-09,-1.57761e-13,-111446,-38.3678], Tmin=(1086.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-903.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO)"""),
)

species(
    label = 'CC(F)(OF)C(F)C(=O)F(8362)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-933.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,333,384,448,608,1254,1480,206,363,518,1175,1317,1481,3059,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,202.243,202.249,202.276],'cm^-1')),
        HinderedRotor(inertia=(0.963985,'amu*angstrom^2'), symmetry=1, barrier=(27.9791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.963737,'amu*angstrom^2'), symmetry=1, barrier=(27.9782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474949,'amu*angstrom^2'), symmetry=1, barrier=(13.7866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474839,'amu*angstrom^2'), symmetry=1, barrier=(13.7865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0172308,0.0953287,-0.000123057,8.41676e-08,-2.34466e-11,-112193,28.6593], Tmin=(100,'K'), Tmax=(867.744,'K')), NASAPolynomial(coeffs=[13.157,0.0347571,-1.83491e-05,3.72083e-09,-2.69023e-13,-114473,-32.8683], Tmin=(867.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-933.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(COCsFO)"""),
)

species(
    label = 'FOF(328)',
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
    label = 'CC(F)(F)C=C=O(1678)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-524.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.433647,'amu*angstrom^2'), symmetry=1, barrier=(9.97039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433288,'amu*angstrom^2'), symmetry=1, barrier=(9.96215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69845,0.0550453,-6.0104e-05,3.71272e-08,-9.69954e-12,-63063.4,19.3548], Tmin=(100,'K'), Tmax=(907.296,'K')), NASAPolynomial(coeffs=[8.14529,0.0266235,-1.31159e-05,2.60154e-09,-1.86343e-13,-64233.3,-11.1203], Tmin=(907.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-524.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'CC(F)=CC(=O)F(8363)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-564.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,323,467,575,827,1418,3010,987.5,1337.5,450,1655,255,533,799,832,1228,180],'cm^-1')),
        HinderedRotor(inertia=(0.218387,'amu*angstrom^2'), symmetry=1, barrier=(5.02116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217836,'amu*angstrom^2'), symmetry=1, barrier=(5.00847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76762,0.0535036,-6.01124e-05,4.03004e-08,-1.15982e-11,-67826.6,19.7443], Tmin=(100,'K'), Tmax=(824.712,'K')), NASAPolynomial(coeffs=[7.04586,0.0279039,-1.35525e-05,2.66413e-09,-1.89619e-13,-68697.2,-4.70303], Tmin=(824.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cd-Cd(CO)H) + group(COCFO)"""),
)

species(
    label = 'CC(F)(F)C=C(F)OOF(8364)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {4,S} {5,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-618.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,274,345,380,539,705,1166,1213,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.633566,'amu*angstrom^2'), symmetry=1, barrier=(14.5669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632693,'amu*angstrom^2'), symmetry=1, barrier=(14.5469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.203446,0.103364,-0.000170712,1.55973e-07,-5.65657e-11,-74224.2,31.6252], Tmin=(100,'K'), Tmax=(790.889,'K')), NASAPolynomial(coeffs=[8.73311,0.0416379,-2.22925e-05,4.44019e-09,-3.13394e-13,-75120.8,-6.12389], Tmin=(790.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-618.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFO)"""),
)

species(
    label = 'CC(F)(F)OC(F)=COF(8365)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {14,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-855.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,223,363,546,575,694,1179,1410,2750,2800,2850,1350,1500,750,1050,1375,1000,326,540,652,719,1357,3010,987.5,1337.5,450,1655,309.578,309.82,309.821,309.876],'cm^-1')),
        HinderedRotor(inertia=(0.143909,'amu*angstrom^2'), symmetry=1, barrier=(9.76771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354818,'amu*angstrom^2'), symmetry=1, barrier=(24.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294793,'amu*angstrom^2'), symmetry=1, barrier=(20.0295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354433,'amu*angstrom^2'), symmetry=1, barrier=(24.1537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.177963,0.0961393,-0.000115944,7.0335e-08,-1.70144e-11,-102724,28.416], Tmin=(100,'K'), Tmax=(1002.75,'K')), NASAPolynomial(coeffs=[16.8083,0.0283805,-1.45835e-05,2.94652e-09,-2.13442e-13,-106131,-53.5792], Tmin=(1002.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-855.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFFO) + group(Cs-CsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CC(F)(F)[C](OF)C([O])F(8366)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {6,S} {10,S} {14,S}
10 C u1 p0 c0 {5,S} {7,S} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-605.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,2750,2800,2850,1350,1500,750,1050,1375,1000,391,562,707,872,1109,1210,1289,3137,360,370,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.270219,'amu*angstrom^2'), symmetry=1, barrier=(6.21286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270237,'amu*angstrom^2'), symmetry=1, barrier=(6.21328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270242,'amu*angstrom^2'), symmetry=1, barrier=(6.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270239,'amu*angstrom^2'), symmetry=1, barrier=(6.21332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.773831,0.117418,-0.000204527,1.83464e-07,-6.34896e-11,-72649.2,33.1181], Tmin=(100,'K'), Tmax=(837.869,'K')), NASAPolynomial(coeffs=[11.4726,0.036375,-1.90177e-05,3.69768e-09,-2.5532e-13,-73908.8,-19.068], Tmin=(837.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-605.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJO)"""),
)

species(
    label = '[CH2]C(F)(F)C(OF)C([O])F(8367)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {7,S} {12,S}
10 C u1 p0 c0 {8,S} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-573.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,3000,3100,440,815,1455,1000,235.756,235.757,235.757],'cm^-1')),
        HinderedRotor(inertia=(0.17659,'amu*angstrom^2'), symmetry=1, barrier=(6.96506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17659,'amu*angstrom^2'), symmetry=1, barrier=(6.96505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176589,'amu*angstrom^2'), symmetry=1, barrier=(6.96505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638756,'amu*angstrom^2'), symmetry=1, barrier=(25.1938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.663067,0.111316,-0.000170033,1.32357e-07,-3.99981e-11,-68805.3,33.9561], Tmin=(100,'K'), Tmax=(692.656,'K')), NASAPolynomial(coeffs=[14.8013,0.0311193,-1.60863e-05,3.17101e-09,-2.23269e-13,-71166.1,-36.549], Tmin=(692.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsHHH) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C(OF)[C](O)F(8368)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 C u1 p0 c0 {8,S} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-605.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,395,473,707,1436,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1966,0.125894,-0.000218907,1.9107e-07,-6.44019e-11,-72640.2,35.7505], Tmin=(100,'K'), Tmax=(834.425,'K')), NASAPolynomial(coeffs=[14.5539,0.0319838,-1.6999e-05,3.31777e-09,-2.29317e-13,-74627.9,-33.5457], Tmin=(834.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-605.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsHHH) + radical(CsCsF1sO2s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CC(F)(F)C(OF)=C(O)F(8369)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,S} {14,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {5,S} {7,S} {10,D}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-847.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,274,345,380,539,705,1166,1213,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.240652,'amu*angstrom^2'), symmetry=1, barrier=(5.53307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24074,'amu*angstrom^2'), symmetry=1, barrier=(5.53508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240765,'amu*angstrom^2'), symmetry=1, barrier=(5.53567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21096,'amu*angstrom^2'), symmetry=1, barrier=(27.8424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (160.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.65418,0.115885,-0.000203919,1.88298e-07,-6.73011e-11,-101735,29.945], Tmin=(100,'K'), Tmax=(823.191,'K')), NASAPolynomial(coeffs=[9.65635,0.0409373,-2.20745e-05,4.36192e-09,-3.04637e-13,-102590,-12.6777], Tmin=(823.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-847.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'C[C](F)C(OF)C(=O)F(8370)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {1,S} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-538.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,486,617,768,1157,1926,347.934,348.052,1693.82],'cm^-1')),
        HinderedRotor(inertia=(0.0845525,'amu*angstrom^2'), symmetry=1, barrier=(7.24295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0843313,'amu*angstrom^2'), symmetry=1, barrier=(7.25594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209483,'amu*angstrom^2'), symmetry=1, barrier=(17.9928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209178,'amu*angstrom^2'), symmetry=1, barrier=(17.9809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696625,0.0790972,-9.88593e-05,6.72413e-08,-1.88247e-11,-64626.8,31.3326], Tmin=(100,'K'), Tmax=(860.262,'K')), NASAPolynomial(coeffs=[10.8749,0.0317681,-1.63292e-05,3.28062e-09,-2.36114e-13,-66377.9,-16.2393], Tmin=(860.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFH) + group(Cs-CsHHH) + group(COCsFO) + radical(CsCsCsF1s)"""),
)

species(
    label = 'CC(F)(F)C([C]=O)OF(8371)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-546.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,222,329,445,522,589,1214,1475,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,385.679,385.718],'cm^-1')),
        HinderedRotor(inertia=(0.118742,'amu*angstrom^2'), symmetry=1, barrier=(12.5422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118836,'amu*angstrom^2'), symmetry=1, barrier=(12.5427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118835,'amu*angstrom^2'), symmetry=1, barrier=(12.5438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118791,'amu*angstrom^2'), symmetry=1, barrier=(12.5436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372725,0.0834744,-0.000101879,6.32141e-08,-1.56465e-11,-65632.6,30.1284], Tmin=(100,'K'), Tmax=(981.393,'K')), NASAPolynomial(coeffs=[14.7213,0.0249931,-1.24965e-05,2.49708e-09,-1.79779e-13,-68449,-38.8262], Tmin=(981.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-546.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(F)(F)C([O])C(=O)F(8372)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-851.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,180.878,1300.94,4000],'cm^-1')),
        HinderedRotor(inertia=(0.69409,'amu*angstrom^2'), symmetry=1, barrier=(15.9585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69421,'amu*angstrom^2'), symmetry=1, barrier=(15.9612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694032,'amu*angstrom^2'), symmetry=1, barrier=(15.9572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782336,0.0755372,-8.74785e-05,5.3834e-08,-1.35365e-11,-102299,29.3759], Tmin=(100,'K'), Tmax=(955.746,'K')), NASAPolynomial(coeffs=[11.9046,0.028988,-1.44212e-05,2.87368e-09,-2.06442e-13,-104425,-23.7789], Tmin=(955.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-851.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(COCsFO) + radical(C=OCOJ)"""),
)

species(
    label = '[O]F(127)',
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
    label = 'CC(F)(F)[CH]C(=O)F(8373)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {4,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-763.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,197.588,975.582],'cm^-1')),
        HinderedRotor(inertia=(0.322911,'amu*angstrom^2'), symmetry=1, barrier=(8.47668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437006,'amu*angstrom^2'), symmetry=1, barrier=(11.8962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35448,'amu*angstrom^2'), symmetry=1, barrier=(35.8387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15513,0.067577,-7.54417e-05,4.55418e-08,-1.14056e-11,-91695,22.8299], Tmin=(100,'K'), Tmax=(951.764,'K')), NASAPolynomial(coeffs=[10.3418,0.028968,-1.45932e-05,2.92032e-09,-2.10222e-13,-93443.7,-21.0364], Tmin=(951.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(COCsFO) + radical(CCJC=O)"""),
)

species(
    label = 'O=C(F)[CH]OF(4318)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-323.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,611,648,830,1210,1753,180,783.373],'cm^-1')),
        HinderedRotor(inertia=(0.41351,'amu*angstrom^2'), symmetry=1, barrier=(9.50742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.48155,'amu*angstrom^2'), symmetry=1, barrier=(57.0558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01266,0.0470476,-7.14345e-05,5.69319e-08,-1.80482e-11,-38878.9,18.3743], Tmin=(100,'K'), Tmax=(773.894,'K')), NASAPolynomial(coeffs=[8.45195,0.0137648,-6.92359e-06,1.35885e-09,-9.56141e-14,-39875.5,-11.0409], Tmin=(773.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(OCJC=O)"""),
)

species(
    label = 'CH3-CF2(65)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-315.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.203516,'amu*angstrom^2'), symmetry=1, barrier=(4.67924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.042,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97026,0.00248557,9.36596e-05,-2.62618e-07,2.48555e-10,-37947.6,9.23011], Tmin=(10,'K'), Tmax=(267.682,'K')), NASAPolynomial(coeffs=[2.69162,0.021592,-1.34044e-05,4.02185e-09,-4.66099e-13,-37879.1,13.7136], Tmin=(267.682,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-315.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3(19)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C(F)C(OF)[C](F)F(2503)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-700.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,190,488,555,1236,1407,486,617,768,1157,1926,280.465,284.148,1409.28],'cm^-1')),
        HinderedRotor(inertia=(0.206174,'amu*angstrom^2'), symmetry=1, barrier=(11.4419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359689,'amu*angstrom^2'), symmetry=1, barrier=(20.063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353371,'amu*angstrom^2'), symmetry=1, barrier=(19.9822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.853719,0.0744802,-0.000101767,6.95983e-08,-1.89229e-11,-84163.8,29.2716], Tmin=(100,'K'), Tmax=(896.877,'K')), NASAPolynomial(coeffs=[12.9879,0.0203625,-1.12562e-05,2.31906e-09,-1.69008e-13,-86340.3,-27.9479], Tmin=(896.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-700.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + radical(CsCsF1sF1s)"""),
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
    label = 'CC(F)(F)[CH]OF(364)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-408.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,205.682,205.682],'cm^-1')),
        HinderedRotor(inertia=(0.00398481,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502912,'amu*angstrom^2'), symmetry=1, barrier=(15.0976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502908,'amu*angstrom^2'), symmetry=1, barrier=(15.0976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07254,0.06931,-9.07284e-05,6.2569e-08,-1.73961e-11,-48980.1,21.9035], Tmin=(100,'K'), Tmax=(873.645,'K')), NASAPolynomial(coeffs=[11.1073,0.0233654,-1.1844e-05,2.37327e-09,-1.70586e-13,-50733.5,-25.153], Tmin=(873.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-408.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
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
    label = 'CC(F)(F)[C](OF)C(=O)F(8374)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 C u0 p0 c0 {3,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-780.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,611,648,830,1210,1753,233.287,233.294,1869.21],'cm^-1')),
        HinderedRotor(inertia=(0.231063,'amu*angstrom^2'), symmetry=1, barrier=(8.92874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231194,'amu*angstrom^2'), symmetry=1, barrier=(8.92887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231147,'amu*angstrom^2'), symmetry=1, barrier=(8.92891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519534,'amu*angstrom^2'), symmetry=1, barrier=(20.0767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0265962,0.0969729,-0.000143003,1.06493e-07,-2.91838e-11,-93751.6,31.919], Tmin=(100,'K'), Tmax=(628.954,'K')), NASAPolynomial(coeffs=[12.0745,0.0323072,-1.72946e-05,3.47035e-09,-2.47233e-13,-95503.6,-22.4983], Tmin=(628.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-780.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(COCsFO) + radical(C2CsJO)"""),
)

species(
    label = '[CH2]C(F)(F)C(OF)C(=O)F(8375)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-748.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,486,617,768,1157,1926,3000,3100,440,815,1455,1000,180,730.102,1374.46],'cm^-1')),
        HinderedRotor(inertia=(0.828775,'amu*angstrom^2'), symmetry=1, barrier=(19.0552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264558,'amu*angstrom^2'), symmetry=1, barrier=(6.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828175,'amu*angstrom^2'), symmetry=1, barrier=(19.0414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826853,'amu*angstrom^2'), symmetry=1, barrier=(19.011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0222468,0.0928883,-0.000119457,7.72565e-08,-1.989e-11,-89902.9,33.1328], Tmin=(100,'K'), Tmax=(945.135,'K')), NASAPolynomial(coeffs=[15.7289,0.0264145,-1.39587e-05,2.84135e-09,-2.06256e-13,-92871.9,-41.7566], Tmin=(945.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-748.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(COCsFO) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'OC(F)=COF(713)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {2,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {4,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-376.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,326,540,652,719,1357,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.996132,'amu*angstrom^2'), symmetry=1, barrier=(22.903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995009,'amu*angstrom^2'), symmetry=1, barrier=(22.8772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72754,0.0499512,-6.1217e-05,3.67587e-08,-8.61434e-12,-45145.5,17.7213], Tmin=(100,'K'), Tmax=(1046.65,'K')), NASAPolynomial(coeffs=[11.9808,0.0107666,-5.06019e-06,9.89726e-10,-7.07149e-14,-47291.8,-32.2122], Tmin=(1046.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CC(F)(F)C(=O)C(=O)F(4488)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {4,D} {6,S} {9,S}
9  C u0 p0 c0 {3,S} {5,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-976.105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,348,427,405,1245,1236,1280,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,286,619,818,1246,1924,323.155,323.164],'cm^-1')),
        HinderedRotor(inertia=(0.114325,'amu*angstrom^2'), symmetry=1, barrier=(8.48059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114357,'amu*angstrom^2'), symmetry=1, barrier=(8.48226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259996,'amu*angstrom^2'), symmetry=1, barrier=(19.2825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749806,0.0780058,-0.000105608,7.8053e-08,-2.36061e-11,-117287,25.4864], Tmin=(100,'K'), Tmax=(801.127,'K')), NASAPolynomial(coeffs=[10.3806,0.0299189,-1.55704e-05,3.12664e-09,-2.24312e-13,-118830,-18.8411], Tmin=(801.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-976.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO)"""),
)

species(
    label = 'CC(F)=C(OF)C(=O)F(8376)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {9,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-586.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,323,467,575,827,1418,350,440,435,1725,255,533,799,832,1228,180,1739.17],'cm^-1')),
        HinderedRotor(inertia=(0.307104,'amu*angstrom^2'), symmetry=1, barrier=(7.06093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308505,'amu*angstrom^2'), symmetry=1, barrier=(7.09314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306208,'amu*angstrom^2'), symmetry=1, barrier=(7.04032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371105,0.0899915,-0.000153567,1.41084e-07,-5.05445e-11,-70423.1,26.4789], Tmin=(100,'K'), Tmax=(819.069,'K')), NASAPolynomial(coeffs=[8.02986,0.0343025,-1.80922e-05,3.55788e-09,-2.48349e-13,-71064.3,-5.19704], Tmin=(819.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-586.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-Cds(Cds-O2d)O2s) + group(COCFO)"""),
)

species(
    label = 'CC(F)(F)C(=C=O)OF(8337)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-506.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,274,345,380,539,705,1166,1213,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.308674,'amu*angstrom^2'), symmetry=1, barrier=(7.09702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308563,'amu*angstrom^2'), symmetry=1, barrier=(7.09446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29105,'amu*angstrom^2'), symmetry=1, barrier=(29.6838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0246819,0.0992633,-0.000173063,1.55887e-07,-5.41305e-11,-60830.4,25.7736], Tmin=(100,'K'), Tmax=(839.547,'K')), NASAPolynomial(coeffs=[9.98009,0.0316118,-1.64861e-05,3.20051e-09,-2.20769e-13,-61806.1,-16.5495], Tmin=(839.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(O2d-Cdd) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=C(F)C(OF)C(=O)F(8377)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-598.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,323,467,575,827,1418,486,617,768,1157,1926,2950,3100,1380,975,1025,1650,180,409.736,900.024],'cm^-1')),
        HinderedRotor(inertia=(0.152308,'amu*angstrom^2'), symmetry=1, barrier=(3.50185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753313,'amu*angstrom^2'), symmetry=1, barrier=(17.3201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.753322,'amu*angstrom^2'), symmetry=1, barrier=(17.3203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849713,0.0714613,-7.83649e-05,4.34699e-08,-9.68458e-12,-71849.7,28.7652], Tmin=(100,'K'), Tmax=(1080.4,'K')), NASAPolynomial(coeffs=[13.8128,0.0234667,-1.17295e-05,2.35152e-09,-1.69807e-13,-74650.7,-34.7766], Tmin=(1080.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + group(COCsFO) + group(Cds-CdsHH)"""),
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
    E0 = (-483.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-147.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (35.6932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (11.8669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (3.23534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-76.4269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-227.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (82.0983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (90.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-114.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-334.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-175.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-140.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-172.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-271.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-57.8756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-66.3569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-371.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-253.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-225.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-156.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-190.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-161.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-129.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-257.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-414.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-277.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-194.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-297.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    products = ['HF(38)', 'O=CC(=O)F(4234)', 'CH2CF2(57)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(69.8229,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'CC(F)(F)C(F)OF(5936)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(363.638,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(43)', 'CC(OF)C(=O)F(8360)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.22791e-07,'m^3/(mol*s)'), n=3.59921, Ea=(324.582,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node CC_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3F(39)', 'O=C(F)C([C]F)OF(8132)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(166.325,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)C(F)OF(5789)', 'C[C]F(325)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(170.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(S)(25)', 'O=C(F)C(OF)C(F)F(8361)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(486642,'m^3/(mol*s)'), n=0.225671, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(F)(OF)C(F)C(=O)F(8362)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FOF(328)', 'CC(F)(F)C=C=O(1678)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.75748e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd_HNd;H_OH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FOF(328)', 'CC(F)=CC(=O)F(8363)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC(F)(F)C=C(F)OOF(8364)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(95.8105,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC(F)(F)OC(F)=COF(8365)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(112.817,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(F)(F)[C](OF)C([O])F(8366)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(F)(F)C(OF)C([O])F(8367)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(F)(F)C(OF)[C](O)F(8368)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC(F)(F)C(OF)=C(O)F(8369)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(167.84,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'C[C](F)C(OF)C(=O)F(8370)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'CC(F)(F)C([C]=O)OF(8371)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'CC(F)(F)C([O])C(=O)F(8372)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]F(127)', 'CC(F)(F)[CH]C(=O)F(8373)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C(F)[CH]OF(4318)', 'CH3-CF2(65)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(6.11201,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(19)', 'O=C(F)C(OF)[C](F)F(2503)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.29446e+18,'m^3/(mol*s)'), n=-4.19701, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R
Ea raised from -12.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CFO(51)', 'CC(F)(F)[CH]OF(364)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -6.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'CC(F)(F)[C](OF)C(=O)F(8374)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', '[CH2]C(F)(F)C(OF)C(=O)F(8375)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    products = ['CH2CF2(57)', 'OC(F)=COF(713)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.87992e+10,'s^-1'), n=0.653804, Ea=(295.391,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.006151428258272596, var=1.8383731639583567, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-2C-R
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'CC(F)(F)C(=O)C(=O)F(4488)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(435.521,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', 'CC(F)=C(OF)C(=O)F(8376)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(182.452,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', 'CC(F)(F)C(=C=O)OF(8337)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(185.632,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'C=C(F)C(OF)C(=O)F(8377)'],
    products = ['CC(F)(F)C(OF)C(=O)F(8310)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(10.477,'m^3/(mol*s)'), n=1.29062, Ea=(174.028,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_5R!H->F_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_5R!H->F_Ext-4COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #2307',
    isomers = [
        'CC(F)(F)C(OF)C(=O)F(8310)',
    ],
    reactants = [
        ('HF(38)', 'O=CC(=O)F(4234)', 'CH2CF2(57)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2307',
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

