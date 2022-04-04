species(
    label = '[O]OOC([O])F(12632)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-131.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,180,180,837.507,837.548,837.563],'cm^-1')),
        HinderedRotor(inertia=(0.101362,'amu*angstrom^2'), symmetry=1, barrier=(2.33052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101367,'amu*angstrom^2'), symmetry=1, barrier=(2.33062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02817,0.0366543,-3.84973e-05,1.93692e-08,-3.69751e-12,-15774.9,22.6682], Tmin=(100,'K'), Tmax=(1392.18,'K')), NASAPolynomial(coeffs=[12.3025,0.00456309,-1.15058e-06,1.58655e-10,-9.57171e-15,-18386.5,-29.4042], Tmin=(1392.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(ROOJ)"""),
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
    label = 'FC1OO1(314)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-250.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([509,613,660,1171,1360,1414,3084,1076.4,1139.43],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2964.22,'J/mol'), sigma=(4.92574,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=463.00 K, Pc=56.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07975,-0.00670713,7.61072e-05,-1.23288e-07,6.38153e-11,-30112,8.13301], Tmin=(10,'K'), Tmax=(624.171,'K')), NASAPolynomial(coeffs=[2.49608,0.0165225,-1.11531e-05,3.48846e-09,-4.10791e-13,-30169.1,12.9858], Tmin=(624.171,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-250.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FC1OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])F(441)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-145.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,1064.74],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3649.77,'J/mol'), sigma=(5.8153,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.09 K, Pc=42.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4785,0.00278736,3.40141e-05,-4.87395e-08,1.9468e-11,-17414.9,12.3963], Tmin=(100,'K'), Tmax=(944.753,'K')), NASAPolynomial(coeffs=[8.92825,0.00190203,1.90853e-07,-1.26597e-11,-4.31401e-15,-19434.9,-18.8266], Tmin=(944.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sH)) + radical(O2sj(Cs-F1sO2sH))"""),
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
    label = '[O]OC([O])F(366)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-167.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,451,553,637,1069,1180,1265,1301,3056],'cm^-1')),
        HinderedRotor(inertia=(0.491221,'amu*angstrom^2'), symmetry=1, barrier=(11.2941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.73,'J/mol'), sigma=(6.011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.88 K, Pc=39.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73707,0.0289956,-3.52661e-05,2.30449e-08,-6.07483e-12,-20134.2,16.0103], Tmin=(100,'K'), Tmax=(921.451,'K')), NASAPolynomial(coeffs=[7.07196,0.0101777,-4.63275e-06,8.81523e-10,-6.15893e-14,-20933.1,-4.54836], Tmin=(921.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(ROOJ)"""),
)

species(
    label = '[O]OO[CH]F(365)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u1 p2 c0 {3,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (25.1845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,367.813,367.893,369.12,1888.8,4000],'cm^-1')),
        HinderedRotor(inertia=(0.165449,'amu*angstrom^2'), symmetry=1, barrier=(15.9163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479392,'amu*angstrom^2'), symmetry=1, barrier=(45.9795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3349.04,'J/mol'), sigma=(5.47852,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.11 K, Pc=46.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63192,0.0302498,-3.83392e-05,2.46662e-08,-6.21158e-12,3078.15,18.3779], Tmin=(100,'K'), Tmax=(976.55,'K')), NASAPolynomial(coeffs=[8.24696,0.00725037,-3.01159e-06,5.48905e-10,-3.74905e-14,1981.47,-8.57817], Tmin=(976.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHHO) + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'FC1OOOO1(12647)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-257.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14395,-0.00774159,0.000117941,-1.60447e-07,6.41263e-11,-30876.4,17.1347], Tmin=(100,'K'), Tmax=(930.043,'K')), NASAPolynomial(coeffs=[18.6848,-0.00632707,5.57793e-06,-9.95541e-10,5.43361e-14,-36719.1,-72.5835], Tmin=(930.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFHOO) + ring(Cyclopentane)"""),
)

species(
    label = 'FC1OOO1(414)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-211.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17662,-0.000295381,7.50007e-05,-1.03845e-07,4.09606e-11,-25378.2,12.0173], Tmin=(100,'K'), Tmax=(950.177,'K')), NASAPolynomial(coeffs=[14.3675,-0.000486106,1.2312e-06,-1.16945e-10,-4.83923e-15,-29622.9,-52.5461], Tmin=(950.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsFHOO) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(F)OOO(12648)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {6,D}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-563.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14571,0.0236377,2.83401e-05,-6.86128e-08,3.17537e-11,-67745.8,20.7463], Tmin=(100,'K'), Tmax=(932.873,'K')), NASAPolynomial(coeffs=[19.0568,-0.00499856,3.83586e-06,-6.83696e-10,3.81141e-14,-72810.2,-69.8977], Tmin=(932.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO)"""),
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
    label = '[O]OOC=O(12649)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-189.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,511.709,511.723,511.734,511.745],'cm^-1')),
        HinderedRotor(inertia=(0.000643727,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247173,'amu*angstrom^2'), symmetry=1, barrier=(45.9353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00145,0.00571033,5.63213e-05,-8.74895e-08,3.6661e-11,-22722.6,17.7459], Tmin=(100,'K'), Tmax=(925.954,'K')), NASAPolynomial(coeffs=[14.831,-0.00354083,3.51109e-06,-6.5518e-10,3.76059e-14,-26707.4,-48.1028], Tmin=(925.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-OdOsH) + radical(ROOJ)"""),
)

species(
    label = '[O]O[O](129)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 O u1 p2 c0 {1,S}
"""),
    E0 = (192.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([767.866,3898.78,3898.78],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3633.48,'J/mol'), sigma=(5.77645,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.54 K, Pc=42.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47471,0.00468948,-6.33886e-06,3.99268e-09,-7.67833e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.31,'K')), NASAPolynomial(coeffs=[2.72434,0.000671943,1.37794e-06,-3.54953e-10,2.60885e-14,24449.9,18.0452], Tmin=(1873.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'CHFO(47)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OOC(=O)F(6945)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
"""),
    E0 = (-411.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,664,788,1296,1923,589.236,589.239,589.239,589.239,589.244],'cm^-1')),
        HinderedRotor(inertia=(0.116087,'amu*angstrom^2'), symmetry=1, barrier=(28.6018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116088,'amu*angstrom^2'), symmetry=1, barrier=(28.6018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42966,0.0198098,2.44598e-05,-6.0377e-08,2.84584e-11,-49476.3,20.4862], Tmin=(100,'K'), Tmax=(920.984,'K')), NASAPolynomial(coeffs=[17.2534,-0.00628768,4.61034e-06,-8.72831e-10,5.39349e-14,-53830.4,-58.6246], Tmin=(920.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]F(388)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53591,0.00799828,-5.22385e-07,-4.56598e-09,2.08746e-12,-2348.33,9.31393], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97189,0.00388217,-1.62977e-06,3.36434e-10,-2.54185e-14,-3160.19,-3.94933], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]OO[C](O)F(12650)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {1,S} {2,S} {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-181.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06048,0.0401424,-4.84397e-05,2.82256e-08,-6.31287e-12,-21714.3,23.196], Tmin=(100,'K'), Tmax=(1107.14,'K')), NASAPolynomial(coeffs=[11.5068,0.00601401,-2.20152e-06,3.8337e-10,-2.59371e-14,-23806,-23.3384], Tmin=(1107.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(ROOJ) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O][C](F)OOO(12651)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u1 p2 c0 {6,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-77.6159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,586,761,1411,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25692,0.0313721,-1.72649e-05,-5.58502e-09,5.63541e-12,-9266.05,23.2678], Tmin=(100,'K'), Tmax=(982.887,'K')), NASAPolynomial(coeffs=[12.4198,0.0057351,-2.13416e-06,4.26895e-10,-3.32587e-14,-12023.3,-29.4504], Tmin=(982.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.6159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]OO[CH]OF(12652)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {2,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (171.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,180,180,180,655.289,655.292,655.298],'cm^-1')),
        HinderedRotor(inertia=(0.030994,'amu*angstrom^2'), symmetry=1, barrier=(9.44524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104136,'amu*angstrom^2'), symmetry=1, barrier=(3.17402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0525856,'amu*angstrom^2'), symmetry=1, barrier=(16.0237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51139,0.054079,-8.01786e-05,5.51084e-08,-1.43249e-11,20774,22.4271], Tmin=(100,'K'), Tmax=(959.798,'K')), NASAPolynomial(coeffs=[13.6581,0.00345378,-1.05522e-06,1.46676e-10,-8.00853e-15,18442.4,-35.6748], Tmin=(959.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O][CH]OOOF(12653)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {1,S} {3,S}
5 O u1 p2 c0 {6,S}
6 C u1 p0 c0 {2,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (208.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,3025,407.5,1350,352.5,503.556,504.075,506.474,507.356,508.268],'cm^-1')),
        HinderedRotor(inertia=(0.000654152,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000660757,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000656082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97613,0.0514698,-9.38144e-05,9.00777e-08,-3.32868e-11,25102.1,21.6833], Tmin=(100,'K'), Tmax=(818.266,'K')), NASAPolynomial(coeffs=[5.53012,0.0202902,-1.13485e-05,2.27015e-09,-1.59592e-13,24982.7,8.07444], Tmin=(818.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2sFO) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
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
    E0 = (2.6308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (111.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (304.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-87.9467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (70.0066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-70.0862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-13.2659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-43.4071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-72.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-79.5182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (209.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (18.9125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (43.0118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (288.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (314.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OOC([O])F(12632)'],
    products = ['O2(2)', 'FC1OO1(314)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.11041e+14,'s^-1'), n=-0.372113, Ea=(97.6903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OO_intra]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]OC([O])F(366)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[O]OO[CH]F(365)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OOC([O])F(12632)'],
    products = ['FC1OOOO1(12647)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OOC([O])F(12632)'],
    products = ['O(6)', 'FC1OOO1(414)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.67697e+10,'s^-1'), n=0.54, Ea=(165.066,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R3OO;Y_rad_intra;OO] for rate rule [R3OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OOC([O])F(12632)'],
    products = ['O=C(F)OOO(12648)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[O]OOC=O(12649)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(66.4544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]O[O](129)', 'CHFO(47)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.35412e-07,'m^3/(mol*s)'), n=3.42791, Ea=(121.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.9589805767432776, var=2.714361599831737, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]OOC(=O)F(6945)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(91.3469,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O2(2)', '[O]C([O])F(441)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.96604e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(37.3882,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]O[O](129)', '[O][CH]F(388)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OOC([O])F(12632)'],
    products = ['[O]OO[C](O)F(12650)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](F)OOO(12651)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OO[CH]OF(12652)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.4537,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]OOOF(12653)'],
    products = ['[O]OOC([O])F(12632)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(69.7596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #3297',
    isomers = [
        '[O]OOC([O])F(12632)',
    ],
    reactants = [
        ('O2(2)', 'FC1OO1(314)'),
        ('O2(2)', '[O]C([O])F(441)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3297',
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

