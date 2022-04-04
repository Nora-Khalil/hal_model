species(
    label = 'CC(F)C(OF)[C](F)F(5315)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {2,S} {3,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-603.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,250,417,511,1155,1315,1456,3119,2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407,312.196,312.403,312.408],'cm^-1')),
        HinderedRotor(inertia=(0.145415,'amu*angstrom^2'), symmetry=1, barrier=(10.1294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145996,'amu*angstrom^2'), symmetry=1, barrier=(10.1271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147847,'amu*angstrom^2'), symmetry=1, barrier=(10.1344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333971,'amu*angstrom^2'), symmetry=1, barrier=(23.1057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.124659,0.0983029,-0.00014564,1.17275e-07,-3.81461e-11,-72385.5,29.6176], Tmin=(100,'K'), Tmax=(750.96,'K')), NASAPolynomial(coeffs=[11.9196,0.0341358,-1.74435e-05,3.44469e-09,-2.43525e-13,-74194.1,-25.0371], Tmin=(750.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-603.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = 'CH2CHF(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC[C](F)F(1659)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-331.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0186107,'amu*angstrom^2'), symmetry=1, barrier=(5.09918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369651,'amu*angstrom^2'), symmetry=1, barrier=(8.49901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67676,0.0370876,8.23942e-05,-9.17766e-07,1.83952e-09,-39878.1,11.543], Tmin=(10,'K'), Tmax=(207.01,'K')), NASAPolynomial(coeffs=[5.61985,0.0252535,-1.81626e-05,6.06645e-09,-7.59349e-13,-40013.6,3.89849], Tmin=(207.01,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-331.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""FOC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C[C]F(137)',
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
        HinderedRotor(inertia=(0.0264144,'amu*angstrom^2'), symmetry=1, barrier=(3.87782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1971.36,'J/mol'), sigma=(5.118e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50526,0.00692148,1.56045e-05,-2.26115e-08,8.49098e-12,9752.45,9.13904], Tmin=(100,'K'), Tmax=(981.304,'K')), NASAPolynomial(coeffs=[5.2739,0.00958149,-3.54759e-06,6.48822e-10,-4.59773e-14,8930.14,-1.78146], Tmin=(981.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'CC(OF)[C](F)F(6765)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-377.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407,250.084,250.084],'cm^-1')),
        HinderedRotor(inertia=(0.208077,'amu*angstrom^2'), symmetry=1, barrier=(9.23468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208077,'amu*angstrom^2'), symmetry=1, barrier=(9.23468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208077,'amu*angstrom^2'), symmetry=1, barrier=(9.23468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31318,0.0567592,-7.29998e-05,5.3657e-08,-1.64962e-11,-45445.1,12.9281], Tmin=(10,'K'), Tmax=(761.87,'K')), NASAPolynomial(coeffs=[8.67525,0.0286072,-1.75733e-05,5.15692e-09,-5.81528e-13,-46262.2,-11.4824], Tmin=(761.87,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-377.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""CC(OF)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FCC(OF)[C](F)F(6807)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-547.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.275503,'amu*angstrom^2'), symmetry=1, barrier=(6.33436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21771,'amu*angstrom^2'), symmetry=1, barrier=(27.9975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126199,'amu*angstrom^2'), symmetry=1, barrier=(2.90156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.54784,0.055521,0.00025029,-2.6672e-06,6.12618e-09,-65866.2,12.8668], Tmin=(10,'K'), Tmax=(175.287,'K')), NASAPolynomial(coeffs=[6.48714,0.0382156,-2.75138e-05,9.15773e-09,-1.14116e-12,-66045.7,1.6237], Tmin=(175.287,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-547.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""FCC(OF)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC(OF)C(F)[C](F)F(9032)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-583.626,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407,289.334,289.336,289.339],'cm^-1')),
        HinderedRotor(inertia=(0.140931,'amu*angstrom^2'), symmetry=1, barrier=(8.37231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140932,'amu*angstrom^2'), symmetry=1, barrier=(8.37231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463766,'amu*angstrom^2'), symmetry=1, barrier=(27.5508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140931,'amu*angstrom^2'), symmetry=1, barrier=(8.37232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.115691,0.0938773,-0.000130818,9.53254e-08,-2.60816e-11,-70061.9,29.5007], Tmin=(100,'K'), Tmax=(646.421,'K')), NASAPolynomial(coeffs=[11.407,0.034495,-1.7359e-05,3.41097e-09,-2.4073e-13,-71740.8,-21.7418], Tmin=(646.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-583.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsHHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'CC(F)C(F)(F)[CH]OF(9033)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
9  C u1 p0 c0 {5,S} {7,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-622.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331883,0.102209,-0.00014617,1.08736e-07,-3.22431e-11,-74724.4,28.2663], Tmin=(100,'K'), Tmax=(825.105,'K')), NASAPolynomial(coeffs=[14.509,0.0302581,-1.53599e-05,3.03771e-09,-2.15569e-13,-77173.3,-40.4783], Tmin=(825.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsCsFH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJO)"""),
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
    label = 'CC(F)C([C]F)OF-2(9034)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {2,S} {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u2 p0 c0 {3,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-212.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,90,150,250,247,323,377,431,433,1065,1274,411.443,412.467,414.339,1568.99],'cm^-1')),
        HinderedRotor(inertia=(0.000976881,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000989771,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000158506,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768655,'amu*angstrom^2'), symmetry=1, barrier=(9.20119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463413,0.0842708,-0.000116185,8.74212e-08,-2.67552e-11,-25474,25.2954], Tmin=(100,'K'), Tmax=(794.301,'K')), NASAPolynomial(coeffs=[11.0252,0.031083,-1.5743e-05,3.1191e-09,-2.21826e-13,-27151.8,-23.227], Tmin=(794.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFHH) + group(Cs-CsHHH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]F(744)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = 'CC(F)[CH]OF(1522)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {4,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,578.81],'cm^-1')),
        HinderedRotor(inertia=(0.595194,'amu*angstrom^2'), symmetry=1, barrier=(13.6847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0941531,'amu*angstrom^2'), symmetry=1, barrier=(2.16476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595123,'amu*angstrom^2'), symmetry=1, barrier=(13.6831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62434,0.0559216,-6.22993e-05,3.81205e-08,-9.64894e-12,-20484,20.1072], Tmin=(100,'K'), Tmax=(944.894,'K')), NASAPolynomial(coeffs=[9.18645,0.0239087,-1.14788e-05,2.26382e-09,-1.61856e-13,-21913.1,-15.9469], Tmin=(944.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
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
    label = 'CC(F)C=C(F)F(763)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {2,S} {3,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-615.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.335246,'amu*angstrom^2'), symmetry=1, barrier=(7.70798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3354,'amu*angstrom^2'), symmetry=1, barrier=(7.71152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75427,0.0219089,0.000203259,-7.97237e-07,9.09415e-10,-74057.1,11.8947], Tmin=(10,'K'), Tmax=(304.394,'K')), NASAPolynomial(coeffs=[4.48812,0.0403525,-2.60355e-05,8.08495e-09,-9.63543e-13,-74231.9,7.0898], Tmin=(304.394,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-615.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""CC(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC=C(F)F(1650)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-380.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,1091.78],'cm^-1')),
        HinderedRotor(inertia=(0.655419,'amu*angstrom^2'), symmetry=1, barrier=(29.8096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79505,0.0174471,4.03379e-05,-1.08299e-07,6.89768e-11,-45737.6,11.4304], Tmin=(10,'K'), Tmax=(583.249,'K')), NASAPolynomial(coeffs=[5.4314,0.0228414,-1.62698e-05,5.25091e-09,-6.31965e-13,-46211.1,1.99525], Tmin=(583.249,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-380.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""FOCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3-CHF(63)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-89.9142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,334,575,1197,1424,3202],'cm^-1')),
        HinderedRotor(inertia=(0.0832224,'amu*angstrom^2'), symmetry=1, barrier=(1.91345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0516,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(2,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9096,0.00783843,1.54328e-05,-2.14834e-08,7.90131e-12,-10814.1,7.80299], Tmin=(10,'K'), Tmax=(888.916,'K')), NASAPolynomial(coeffs=[2.27498,0.019093,-1.0138e-05,2.62844e-09,-2.67692e-13,-10677.5,14.6301], Tmin=(888.916,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-89.9142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CC(F)C(OF)=C(F)F(9035)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-611.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,284,328,853,1146,1135,1297,3239,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,182,240,577,636,1210,1413,180,2078.17],'cm^-1')),
        HinderedRotor(inertia=(0.218453,'amu*angstrom^2'), symmetry=1, barrier=(5.02267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94356,'amu*angstrom^2'), symmetry=1, barrier=(44.6862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217588,'amu*angstrom^2'), symmetry=1, barrier=(5.00278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414142,0.090047,-0.000150743,1.43314e-07,-5.35858e-11,-73375.3,27.1785], Tmin=(100,'K'), Tmax=(804.26,'K')), NASAPolynomial(coeffs=[5.78711,0.0422914,-2.24471e-05,4.45031e-09,-3.13016e-13,-73559.3,6.65644], Tmin=(804.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'C[CH]C(OF)[C](F)F(9036)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {5,S} {6,S} {13,S}
8  C u1 p0 c0 {1,S} {2,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-190.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,418.765,418.845,420.14],'cm^-1')),
        HinderedRotor(inertia=(0.0742523,'amu*angstrom^2'), symmetry=1, barrier=(9.24628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0742145,'amu*angstrom^2'), symmetry=1, barrier=(9.25261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153681,'amu*angstrom^2'), symmetry=1, barrier=(19.0963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00727168,'amu*angstrom^2'), symmetry=1, barrier=(9.28289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752937,0.0763125,-9.30763e-05,6.12462e-08,-1.64361e-11,-22839.6,29.0743], Tmin=(100,'K'), Tmax=(900.072,'K')), NASAPolynomial(coeffs=[11.4086,0.0289563,-1.41538e-05,2.78817e-09,-1.98694e-13,-24757.7,-21.2112], Tmin=(900.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CCJCO) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'CC(F)C([O])[C](F)F(4772)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-507.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407,274.035,275.535,276.588],'cm^-1')),
        HinderedRotor(inertia=(0.00227223,'amu*angstrom^2'), symmetry=1, barrier=(0.120408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24844,'amu*angstrom^2'), symmetry=1, barrier=(13.1323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247737,'amu*angstrom^2'), symmetry=1, barrier=(13.1807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607007,0.080754,-0.000108321,8.03304e-08,-2.43561e-11,-60881.8,26.635], Tmin=(100,'K'), Tmax=(800.144,'K')), NASAPolynomial(coeffs=[10.5006,0.0312947,-1.56008e-05,3.07672e-09,-2.18461e-13,-62465.1,-18.8898], Tmin=(800.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'CC(F)[CH][C](F)F(770)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {7,S} {12,S}
7  C u1 p0 c0 {2,S} {3,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-365.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,180,855.155],'cm^-1')),
        HinderedRotor(inertia=(0.220282,'amu*angstrom^2'), symmetry=1, barrier=(5.06471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220994,'amu*angstrom^2'), symmetry=1, barrier=(5.08108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219681,'amu*angstrom^2'), symmetry=1, barrier=(5.05091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20289,0.0692373,-0.000107031,1.00317e-07,-3.801e-11,-43834.9,24.9156], Tmin=(100,'K'), Tmax=(776.618,'K')), NASAPolynomial(coeffs=[5.27389,0.0358292,-1.84774e-05,3.6744e-09,-2.60355e-13,-44092,8.71985], Tmin=(776.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FO[CH][C](F)F(1648)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-149.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,190,488,555,1236,1407,287.768,646.375],'cm^-1')),
        HinderedRotor(inertia=(0.194641,'amu*angstrom^2'), symmetry=1, barrier=(11.3648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2772,'amu*angstrom^2'), symmetry=1, barrier=(16.2049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65469,0.0553565,-8.8846e-05,6.9474e-08,-2.11107e-11,-17911.9,20.3653], Tmin=(100,'K'), Tmax=(812.524,'K')), NASAPolynomial(coeffs=[10.8185,0.010241,-5.55348e-06,1.12943e-09,-8.09562e-14,-19401,-21.9418], Tmin=(812.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'CC(F)[C](OF)[C](F)F(9037)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {5,S} {6,S} {9,S}
9  C u1 p0 c0 {2,S} {3,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-422.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.214584,'amu*angstrom^2'), symmetry=1, barrier=(4.93371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218812,'amu*angstrom^2'), symmetry=1, barrier=(5.03093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219542,'amu*angstrom^2'), symmetry=1, barrier=(5.04769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657163,'amu*angstrom^2'), symmetry=1, barrier=(15.1095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.277481,0.106331,-0.000187914,1.71654e-07,-6.02839e-11,-50723.1,30.5755], Tmin=(100,'K'), Tmax=(838.951,'K')), NASAPolynomial(coeffs=[9.73221,0.0350882,-1.84846e-05,3.6036e-09,-2.49117e-13,-51575,-11.0251], Tmin=(838.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(C2CsJO) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = 'F[CH]C(OF)[C](F)F(9038)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
7  C u1 p0 c0 {1,S} {6,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-353.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,334,575,1197,1424,3202,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.278423,'amu*angstrom^2'), symmetry=1, barrier=(6.40149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278314,'amu*angstrom^2'), symmetry=1, barrier=(6.39899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26017,'amu*angstrom^2'), symmetry=1, barrier=(28.9739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218452,0.0981436,-0.000194793,1.87445e-07,-6.73736e-11,-42450.8,26.2754], Tmin=(100,'K'), Tmax=(852.711,'K')), NASAPolynomial(coeffs=[7.89173,0.0295164,-1.66676e-05,3.30284e-09,-2.28554e-13,-42573.1,-2.56455], Tmin=(852.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-O2sCsH)(F1s)(H)) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'C[C](F)C(OF)[C](F)F(6770)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {1,S} {6,S} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-412.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,190,488,555,1236,1407,193.095,193.369,193.832],'cm^-1')),
        HinderedRotor(inertia=(0.227677,'amu*angstrom^2'), symmetry=1, barrier=(5.97729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00457747,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224352,'amu*angstrom^2'), symmetry=1, barrier=(5.98903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706712,'amu*angstrom^2'), symmetry=1, barrier=(18.4462,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.210627,0.104065,-0.000180819,1.64192e-07,-5.7712e-11,-49477.6,31.9249], Tmin=(100,'K'), Tmax=(829.828,'K')), NASAPolynomial(coeffs=[9.80467,0.035023,-1.84813e-05,3.61793e-09,-2.51215e-13,-50424.8,-10.2169], Tmin=(829.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-412.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CsCsCsF1s) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)C(OF)[C](F)F(2828)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-396.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,3000,3100,440,815,1455,1000,286.74,287.99,289.581],'cm^-1')),
        HinderedRotor(inertia=(0.105565,'amu*angstrom^2'), symmetry=1, barrier=(6.19393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219474,'amu*angstrom^2'), symmetry=1, barrier=(13.0658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102795,'amu*angstrom^2'), symmetry=1, barrier=(6.17896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56377,'amu*angstrom^2'), symmetry=1, barrier=(33.6912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.322216,0.103924,-0.000171064,1.47558e-07,-5.01451e-11,-47515.8,31.7221], Tmin=(100,'K'), Tmax=(801.99,'K')), NASAPolynomial(coeffs=[12.1696,0.0314464,-1.64778e-05,3.23861e-09,-2.26378e-13,-49192.3,-23.7475], Tmin=(801.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'CC(F)C(=O)[C](F)F(9039)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {4,D} {5,S} {8,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-652.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,179,346,818,1406,1524,248.174,248.269],'cm^-1')),
        HinderedRotor(inertia=(0.00273349,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169166,'amu*angstrom^2'), symmetry=1, barrier=(7.39604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583976,'amu*angstrom^2'), symmetry=1, barrier=(25.57,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00819,0.0713787,-9.45282e-05,7.10877e-08,-2.20865e-11,-78362,24.2407], Tmin=(100,'K'), Tmax=(778.819,'K')), NASAPolynomial(coeffs=[9.05552,0.0300452,-1.49149e-05,2.93463e-09,-2.08071e-13,-79615.4,-12.5709], Tmin=(778.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-652.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + radical(Csj(CO-CsO2d)(F1s)(F1s))"""),
)

species(
    label = 'CC=C(OF)[C](F)F(9040)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {4,S} {6,D} {8,S}
8  C u1 p0 c0 {1,S} {2,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-279.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,161,297,490,584,780,1358,180],'cm^-1')),
        HinderedRotor(inertia=(0.393232,'amu*angstrom^2'), symmetry=1, barrier=(9.04119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00369945,'amu*angstrom^2'), symmetry=1, barrier=(9.03006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19155,'amu*angstrom^2'), symmetry=1, barrier=(50.3881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742864,0.0809622,-0.000133235,1.22948e-07,-4.44978e-11,-33449.3,26.0832], Tmin=(100,'K'), Tmax=(821.293,'K')), NASAPolynomial(coeffs=[6.50204,0.0354098,-1.80709e-05,3.51748e-09,-2.44539e-13,-33805,3.02612], Tmin=(821.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Csj(Cd-O2sCd)(F1s)(F1s))"""),
)

species(
    label = 'C=CC(OF)[C](F)F(9029)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-262.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,334.559,334.793,334.794],'cm^-1')),
        HinderedRotor(inertia=(0.193464,'amu*angstrom^2'), symmetry=1, barrier=(15.4226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193821,'amu*angstrom^2'), symmetry=1, barrier=(15.4266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193773,'amu*angstrom^2'), symmetry=1, barrier=(15.4248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904236,0.071312,-8.19462e-05,4.84518e-08,-1.15416e-11,-31500.2,27.2997], Tmin=(100,'K'), Tmax=(1012.84,'K')), NASAPolynomial(coeffs=[12.8692,0.0240591,-1.19658e-05,2.38983e-09,-1.72106e-13,-33923.9,-30.577], Tmin=(1012.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'CC(F)C([C]F)OF(5303)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
8  C u0 p1 c0 {2,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-216.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,250,417,511,1155,1315,1456,3119,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,617,898,1187,272.983,272.987],'cm^-1')),
        HinderedRotor(inertia=(0.184126,'amu*angstrom^2'), symmetry=1, barrier=(9.73747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184161,'amu*angstrom^2'), symmetry=1, barrier=(9.7374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184183,'amu*angstrom^2'), symmetry=1, barrier=(9.73759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550885,'amu*angstrom^2'), symmetry=1, barrier=(29.1348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649914,0.0790303,-0.000101027,7.02818e-08,-1.99163e-11,-25970.9,26.3205], Tmin=(100,'K'), Tmax=(855.089,'K')), NASAPolynomial(coeffs=[11.1892,0.0297288,-1.45417e-05,2.85373e-09,-2.02515e-13,-27773.3,-22.8753], Tmin=(855.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'CC(F)[C](OF)C(F)F(5312)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-624.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,522,611,926,1093,1137,1374,1416,3112,360,370,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.361208,'amu*angstrom^2'), symmetry=1, barrier=(8.30488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361322,'amu*angstrom^2'), symmetry=1, barrier=(8.3075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361003,'amu*angstrom^2'), symmetry=1, barrier=(8.30016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806401,'amu*angstrom^2'), symmetry=1, barrier=(18.5407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.23684,0.102547,-0.000165789,1.45056e-07,-5.01746e-11,-75001,29.3832], Tmin=(100,'K'), Tmax=(808.513,'K')), NASAPolynomial(coeffs=[10.5213,0.036243,-1.8513e-05,3.61062e-09,-2.51558e-13,-76313.1,-17.5884], Tmin=(808.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-624.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(C2CsJO)"""),
)

species(
    label = 'C[C](F)C(OF)C(F)F(5314)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
8  C u0 p0 c0 {9,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {3,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-614.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,208.422,208.422,208.423],'cm^-1')),
        HinderedRotor(inertia=(0.264077,'amu*angstrom^2'), symmetry=1, barrier=(8.14043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00388079,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264075,'amu*angstrom^2'), symmetry=1, barrier=(8.14043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.729687,'amu*angstrom^2'), symmetry=1, barrier=(22.4934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.158486,0.100135,-0.000158145,1.36841e-07,-4.7273e-11,-73756,30.692], Tmin=(100,'K'), Tmax=(790.887,'K')), NASAPolynomial(coeffs=[10.5496,0.0362571,-1.85573e-05,3.6365e-09,-2.54635e-13,-75145.7,-16.5341], Tmin=(790.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH2]C(F)C(OF)C(F)F(5316)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {12,S}
9  C u1 p0 c0 {7,S} {13,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-598.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,3000,3100,440,815,1455,1000,288.947,288.948,288.952],'cm^-1')),
        HinderedRotor(inertia=(0.129958,'amu*angstrom^2'), symmetry=1, barrier=(7.69964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.277655,'amu*angstrom^2'), symmetry=1, barrier=(16.4506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129957,'amu*angstrom^2'), symmetry=1, barrier=(7.69968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601407,'amu*angstrom^2'), symmetry=1, barrier=(35.6312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.12767,0.0980983,-0.000140619,1.08193e-07,-3.35539e-11,-71800.2,29.9915], Tmin=(100,'K'), Tmax=(786.723,'K')), NASAPolynomial(coeffs=[12.6425,0.0331785,-1.68563e-05,3.33136e-09,-2.36124e-13,-73809.8,-28.5558], Tmin=(786.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'CC(F)C([O])C(F)(F)F(9041)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-944.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444886,0.0832757,-9.68699e-05,6.0166e-08,-1.52663e-11,-113421,25.9749], Tmin=(100,'K'), Tmax=(948.05,'K')), NASAPolynomial(coeffs=[12.5975,0.0320021,-1.57456e-05,3.11997e-09,-2.23391e-13,-115725,-32.006], Tmin=(948.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-944.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + group(Cs-CsHHH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[CH]C(OF)C(F)(F)F(9042)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u0 p0 c0 {9,S} {11,S} {12,S} {13,S}
9  C u1 p0 c0 {6,S} {8,S} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-627.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (145.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469324,0.0803929,-8.76598e-05,4.97488e-08,-1.14298e-11,-75373.8,28.8414], Tmin=(100,'K'), Tmax=(1046.38,'K')), NASAPolynomial(coeffs=[13.9377,0.0289067,-1.38527e-05,2.7245e-09,-1.94655e-13,-78192.4,-36.7463], Tmin=(1046.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-627.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + group(Cs-CsHHH) + radical(CCJCO)"""),
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
    E0 = (-316.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (33.7584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (216.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (106.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-49.1962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-207.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (95.5085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (98.1569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-224.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-229.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-162.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (117.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-191.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-27.1236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-4.08638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (25.6739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (17.6481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (34.9467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (50.9422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-238.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-82.5685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-100.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (91.4227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-136.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-193.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-242.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-316.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-275.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-153.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['HF(38)', 'O=C[C](F)F(787)', 'CH2CHF(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(51.0506,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FOC[C](F)F(1659)', 'C[C]F(137)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.67744e-07,'m^3/(mol*s)'), n=3.38172, Ea=(48.9826,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'CC(OF)[C](F)F(6765)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.7066e-06,'m^3/(mol*s)'), n=3.19155, Ea=(219.79,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(S)(25)', 'FCC(OF)[C](F)F(6807)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(973284,'m^3/(mol*s)'), n=0.225671, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H_Ext-7R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC(OF)C(F)[C](F)F(9032)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['CC(F)C(F)(F)[CH]OF(9033)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'CC(F)C([C]F)OF-2(9034)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C]F(744)', 'CC(F)[CH]OF(1522)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]F(127)', 'CC(F)C=C(F)F(763)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.20421e+10,'m^3/(mol*s)'), n=-1.2407, Ea=(53.1939,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08524545916686703, var=37.982834433071, Tref=1000.0, N=20, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FOC=C(F)F(1650)', 'CH3-CHF(63)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.000143481,'m^3/(mol*s)'), n=2.73966, Ea=(5.41394,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24334233861756036, var=1.6451213851449848, Tref=1000.0, N=114, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'CC(F)C(OF)=C(F)F(9035)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(168,'m^3/(mol*s)'), n=1.64, Ea=(1.12922,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C_Ext-7R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C[CH]C(OF)[C](F)F(9036)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -38.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'CC(F)C([O])[C](F)F(4772)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(7.3568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(127)', 'CC(F)[CH][C](F)F(770)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FO[CH][C](F)F(1648)', 'CH3-CHF(63)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -6.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'CC(F)[C](OF)[C](F)F(9037)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(1.35813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH3(19)', 'F[CH]C(OF)[C](F)F(9038)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -14.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'C[C](F)C(OF)[C](F)F(6770)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.261562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[CH2]C(F)C(OF)[C](F)F(2828)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -1.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'CC(F)C(=O)[C](F)F(9039)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(459.241,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'CC=C(OF)[C](F)F(9040)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(242.129,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'C=CC(OF)[C](F)F(9029)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(25.0352,'m^3/(mol*s)'), n=1.25316, Ea=(208.307,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_N-5R!H->F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F(37)', 'CC(F)C([C]F)OF(5303)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(518506,'m^3/(mol*s)'), n=0.4717, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F',), comment="""Estimated from node Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CF2(43)', 'CC(F)[CH]OF(1522)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.38543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['CC(F)[C](OF)C(F)F(5312)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.06381e+10,'s^-1'), n=0.647667, Ea=(174.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['C[C](F)C(OF)C(F)F(5314)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.37585e+08,'s^-1'), n=0.919844, Ea=(125.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(F)C(OF)C(F)F(5316)'],
    products = ['CC(F)C(OF)[C](F)F(5315)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(203.145,'s^-1'), n=2.83, Ea=(46.4821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['CC(F)C([O])C(F)(F)F(9041)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.5808,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC(F)C(OF)[C](F)F(5315)'],
    products = ['C[CH]C(OF)C(F)(F)F(9042)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(214.493,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #2964',
    isomers = [
        'CC(F)C(OF)[C](F)F(5315)',
    ],
    reactants = [
        ('HF(38)', 'O=C[C](F)F(787)', 'CH2CHF(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2964',
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

