species(
    label = 'C=C(CF)C(F)(F)OF(9176)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-694.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,265.575,277.279,283.71],'cm^-1')),
        HinderedRotor(inertia=(0.16824,'amu*angstrom^2'), symmetry=1, barrier=(8.47467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155274,'amu*angstrom^2'), symmetry=1, barrier=(8.38597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534483,'amu*angstrom^2'), symmetry=1, barrier=(27.3736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30749,0.0885658,-0.000124295,9.48202e-08,-2.94566e-11,-83362.8,26.0941], Tmin=(100,'K'), Tmax=(782.191,'K')), NASAPolynomial(coeffs=[11.2433,0.0326433,-1.70559e-05,3.42161e-09,-2.44965e-13,-85073.6,-23.9786], Tmin=(782.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-694.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
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
    label = 'C=C=CF(6101)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(CF)OF(581)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-220.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,227.156,227.247],'cm^-1')),
        HinderedRotor(inertia=(0.0961877,'amu*angstrom^2'), symmetry=1, barrier=(3.52341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148211,'amu*angstrom^2'), symmetry=1, barrier=(20.4618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66033,0.0328402,-1.66216e-05,-5.98148e-10,2.21535e-12,-26504.5,11.9746], Tmin=(10,'K'), Tmax=(1011.53,'K')), NASAPolynomial(coeffs=[7.48837,0.0241566,-1.33153e-05,3.53053e-09,-3.6402e-13,-27609.1,-8.16942], Tmin=(1011.53,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-220.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(CF)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOF(499)',
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
    label = 'C=C([C]F)CF(9182)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p1 c0 {2,S} {4,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-45.5103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,315,622,1128,180],'cm^-1')),
        HinderedRotor(inertia=(0.630444,'amu*angstrom^2'), symmetry=1, barrier=(14.4952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0403116,'amu*angstrom^2'), symmetry=1, barrier=(24.373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0713,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22591,0.041419,-3.22591e-05,1.26084e-08,-2.06327e-12,-5412.11,17.5742], Tmin=(100,'K'), Tmax=(1373.91,'K')), NASAPolynomial(coeffs=[9.05101,0.0215485,-1.0565e-05,2.08175e-09,-1.47816e-13,-7287.53,-17.521], Tmin=(1373.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.5103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(CJ2_singlet-FC)"""),
)

species(
    label = 'H2(8)',
    structure = adjacencyList("""1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-8.60349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3765.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (2.01594,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(496.376,'J/mol'), sigma=(2.8327,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212708,-2.78619e-07,3.40262e-10,-7.7602e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.09,'K')), NASAPolynomial(coeffs=[2.78813,0.000587684,1.5899e-07,-5.52698e-11,4.34281e-15,-596.125,0.112931], Tmin=(1959.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C([C]F)C(F)(F)OF(9183)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {9,S}
8  C u0 p0 c0 {7,D} {10,S} {11,S}
9  C u0 p1 c0 {4,S} {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-355.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,350,440,435,1725,2950,3100,1380,975,1025,1650,315,622,1128,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.148794,'amu*angstrom^2'), symmetry=1, barrier=(3.42106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28896,'amu*angstrom^2'), symmetry=1, barrier=(29.6357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866573,'amu*angstrom^2'), symmetry=1, barrier=(19.9242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494773,0.0840619,-0.000120161,8.82298e-08,-2.59616e-11,-42666.2,25.4154], Tmin=(100,'K'), Tmax=(828.423,'K')), NASAPolynomial(coeffs=[12.6228,0.0255017,-1.41265e-05,2.89848e-09,-2.10177e-13,-44675.6,-30.8122], Tmin=(828.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(CJ2_singlet-FC)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'C=CC(F)(F)OF(4807)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-472.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.822754,'amu*angstrom^2'), symmetry=1, barrier=(18.9167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821372,'amu*angstrom^2'), symmetry=1, barrier=(18.885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77039,0.0151198,0.000164752,-4.29255e-07,3.13787e-10,-56873.6,11.9516], Tmin=(10,'K'), Tmax=(490.439,'K')), NASAPolynomial(coeffs=[6.79229,0.0337748,-2.47412e-05,8.35297e-09,-1.04968e-12,-57690.7,-5.78357], Tmin=(490.439,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-472.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDCC(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(F)C(F)(F)OF(4967)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-650.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.493364,'amu*angstrom^2'), symmetry=1, barrier=(11.3434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815502,'amu*angstrom^2'), symmetry=1, barrier=(18.75,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.68127,0.0224957,0.000176709,-5.29719e-07,4.28239e-10,-78282.4,12.8424], Tmin=(10,'K'), Tmax=(453.826,'K')), NASAPolynomial(coeffs=[7.96789,0.0342935,-2.61585e-05,9.0192e-09,-1.14616e-12,-79182,-10.0767], Tmin=(453.826,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-650.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FCC(CF)=C(F)OF(9184)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-620.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,166,258,854,1008,1066,1142,1198,1304,1253,1397,1428,1520,3056,3148,3059,3251,350,440,435,1725,326,540,652,719,1357,180,180,2016.52],'cm^-1')),
        HinderedRotor(inertia=(0.241938,'amu*angstrom^2'), symmetry=1, barrier=(5.56263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241978,'amu*angstrom^2'), symmetry=1, barrier=(5.56356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70775,'amu*angstrom^2'), symmetry=1, barrier=(39.2645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345626,0.0910815,-0.000151962,1.42167e-07,-5.23621e-11,-74462.5,28.0673], Tmin=(100,'K'), Tmax=(805.169,'K')), NASAPolynomial(coeffs=[6.71227,0.0401815,-2.12357e-05,4.20179e-09,-2.95151e-13,-74863.1,2.61069], Tmin=(805.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-620.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFO)"""),
)

species(
    label = 'FCC(COF)=C(F)F(9185)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-619.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,312.308,313.979,314.55],'cm^-1')),
        HinderedRotor(inertia=(0.132135,'amu*angstrom^2'), symmetry=1, barrier=(9.12701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00219088,'amu*angstrom^2'), symmetry=1, barrier=(9.13208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.535106,'amu*angstrom^2'), symmetry=1, barrier=(37.3989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75393,0.0802371,-0.000101735,6.18232e-08,-9.19113e-12,-74446.1,27.0655], Tmin=(100,'K'), Tmax=(587.31,'K')), NASAPolynomial(coeffs=[9.38986,0.0343117,-1.7366e-05,3.42772e-09,-2.42778e-13,-75682.8,-11.8946], Tmin=(587.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF)"""),
)

species(
    label = '[CH2]C([CH]F)C(F)(F)OF(9186)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  C u1 p0 c0 {6,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-429.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,322.78,323.779,2404.06],'cm^-1')),
        HinderedRotor(inertia=(0.11576,'amu*angstrom^2'), symmetry=1, barrier=(8.58633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115708,'amu*angstrom^2'), symmetry=1, barrier=(8.63149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.634915,'amu*angstrom^2'), symmetry=1, barrier=(47.9151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635686,'amu*angstrom^2'), symmetry=1, barrier=(47.9046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15985,0.0925086,-0.000130741,9.24314e-08,-2.32741e-11,-51578.3,29.8373], Tmin=(100,'K'), Tmax=(648.543,'K')), NASAPolynomial(coeffs=[12.3531,0.0298107,-1.46538e-05,2.83354e-09,-1.97632e-13,-53422.9,-25.7358], Tmin=(648.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-429.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCFFO) + group(CsCsFHH) + group(Cs-CsHHH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Isobutyl)"""),
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
    label = '[CH2]C(CF)=C(F)OF(9187)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {6,D}
8  C u1 p0 c0 {6,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-293.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,326,540,652,719,1357,3000,3100,440,815,1455,1000,180,839.85],'cm^-1')),
        HinderedRotor(inertia=(0.130611,'amu*angstrom^2'), symmetry=1, barrier=(3.003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129664,'amu*angstrom^2'), symmetry=1, barrier=(2.98123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25362,'amu*angstrom^2'), symmetry=1, barrier=(51.8151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26047,0.0697339,-7.70274e-05,2.43948e-08,1.54027e-11,-35165.4,23.5138], Tmin=(100,'K'), Tmax=(531.76,'K')), NASAPolynomial(coeffs=[8.17668,0.0330668,-1.69184e-05,3.34906e-09,-2.37361e-13,-36118.1,-7.52661], Tmin=(531.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C(F)(F)OF(9188)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-367.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,256.995],'cm^-1')),
        HinderedRotor(inertia=(1.56458,'amu*angstrom^2'), symmetry=1, barrier=(35.9728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211352,'amu*angstrom^2'), symmetry=1, barrier=(8.13139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01521,'amu*angstrom^2'), symmetry=1, barrier=(23.3417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742459,0.0747496,-8.75542e-05,5.26786e-08,-1.27222e-11,-44045.3,22.4663], Tmin=(100,'K'), Tmax=(1001.97,'K')), NASAPolynomial(coeffs=[13.3932,0.0242458,-1.19469e-05,2.37247e-09,-1.70278e-13,-46580.4,-38.591], Tmin=(1001.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(CF)C([O])(F)F(9189)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-580.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.743082,'amu*angstrom^2'), symmetry=1, barrier=(17.0849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865659,'amu*angstrom^2'), symmetry=1, barrier=(19.9032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18502,0.0627566,-6.05145e-05,2.98641e-08,-5.99367e-12,-69736.8,24.2687], Tmin=(100,'K'), Tmax=(1182.75,'K')), NASAPolynomial(coeffs=[12.5951,0.0241677,-1.15744e-05,2.27824e-09,-1.62735e-13,-72435.9,-32.6934], Tmin=(1182.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(O2sj(Cs-F1sF1sCd))"""),
)

species(
    label = '[O]F(126)',
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
    label = '[CH2]C(CF)=C(F)F(9190)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u1 p0 c0 {5,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-452.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413,345.562],'cm^-1')),
        HinderedRotor(inertia=(0.00143423,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144134,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59277,0.0373765,-1.457e-05,-8.10142e-09,6.06854e-12,-54442.4,13.7997], Tmin=(10,'K'), Tmax=(895.001,'K')), NASAPolynomial(coeffs=[6.94047,0.0308839,-1.78826e-05,4.93899e-09,-5.27335e-13,-55380.9,-3.87465], Tmin=(895.001,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-452.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[CH2]C(CF)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FO[C](F)F(4688)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
"""),
    E0 = (-317.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.523179,'amu*angstrom^2'), symmetry=1, barrier=(12.0289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66912,0.0316874,-4.80467e-05,3.64259e-08,-1.08808e-11,-38142.7,14.3743], Tmin=(100,'K'), Tmax=(821.789,'K')), NASAPolynomial(coeffs=[7.60337,0.00767068,-4.20998e-06,8.64397e-10,-6.26025e-14,-38953.7,-8.46215], Tmin=(821.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsF1sF1sO2s)"""),
)

species(
    label = 'C=[C]CF(6104)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {4,D} {7,S} {8,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (60.6083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.577779,'amu*angstrom^2'), symmetry=1, barrier=(13.2843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66063,0.0349028,-0.000121885,2.76679e-07,-2.18448e-10,7289.62,9.12208], Tmin=(10,'K'), Tmax=(435.813,'K')), NASAPolynomial(coeffs=[2.65835,0.0237529,-1.34702e-05,3.69674e-09,-3.94886e-13,7570.23,15.3422], Tmin=(435.813,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(60.6083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CD[C]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2F(46)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-42.5685,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(33.0141,'amu')),
        NonlinearRotor(inertia=([1.91548,16.2277,17.9803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([576.418,1180.5,1217.62,1485.55,3118.23,3268.88],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03338,-0.00262849,2.74227e-05,-3.89096e-08,1.85259e-11,-5119.82,5.20374], Tmin=(10,'K'), Tmax=(594.366,'K')), NASAPolynomial(coeffs=[2.59024,0.00857266,-4.60348e-06,1.22743e-09,-1.29255e-13,-4974.57,11.194], Tmin=(594.366,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-42.5685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[CH2]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=[C]C(F)(F)OF(4958)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6 C u0 p0 c0 {7,D} {8,S} {9,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-221.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,140,424,499,621,667,843,876,1082,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.776179,'amu*angstrom^2'), symmetry=1, barrier=(17.8459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.775986,'amu*angstrom^2'), symmetry=1, barrier=(17.8414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74843,0.0174867,0.000153749,-4.46933e-07,3.55345e-10,-26640.7,13.4727], Tmin=(10,'K'), Tmax=(458.33,'K')), NASAPolynomial(coeffs=[7.16121,0.0294524,-2.20494e-05,7.52378e-09,-9.50815e-13,-27392,-5.11333], Tmin=(458.33,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-221.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CD[C]C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(=CF)C(F)(F)OF(9191)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {7,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-558.959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180,1393.89],'cm^-1')),
        HinderedRotor(inertia=(0.958135,'amu*angstrom^2'), symmetry=1, barrier=(22.0294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13793,'amu*angstrom^2'), symmetry=1, barrier=(49.1551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.450391,'amu*angstrom^2'), symmetry=1, barrier=(10.3554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32198,0.0860641,-0.000114438,7.79241e-08,-2.11522e-11,-67099.2,25.5283], Tmin=(100,'K'), Tmax=(898.396,'K')), NASAPolynomial(coeffs=[13.9794,0.0252546,-1.29055e-05,2.57839e-09,-1.84873e-13,-69553,-38.8971], Tmin=(898.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-558.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(CF)C(F)(F)OF(9192)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-447.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.4005,'amu*angstrom^2'), symmetry=1, barrier=(9.20828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26753,'amu*angstrom^2'), symmetry=1, barrier=(29.143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400445,'amu*angstrom^2'), symmetry=1, barrier=(9.20702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2908,0.0909616,-0.000134042,9.66258e-08,-2.41044e-11,-53645.5,26.581], Tmin=(100,'K'), Tmax=(615.698,'K')), NASAPolynomial(coeffs=[11.8142,0.0292397,-1.56896e-05,3.14384e-09,-2.23516e-13,-55313.6,-25.4469], Tmin=(615.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-447.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(CF)C(F)(F)OF(9193)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
9  C u0 p1 c0 {6,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-376.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,528,1116,1182,1331,1402,1494,3075,3110,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.067,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00399851,0.0957613,-0.000142408,1.12702e-07,-3.59003e-11,-45183.5,28.0403], Tmin=(100,'K'), Tmax=(766.738,'K')), NASAPolynomial(coeffs=[12.3198,0.0314672,-1.66226e-05,3.32963e-09,-2.37446e-13,-47073.3,-28.141], Tmin=(766.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(CsCFFO) + group(CsCsFHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=C(CF)C(=O)F(7956)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-538.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,255,533,799,832,1228,2950,3100,1380,975,1025,1650,180,784.252],'cm^-1')),
        HinderedRotor(inertia=(0.196261,'amu*angstrom^2'), symmetry=1, barrier=(4.51243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203739,'amu*angstrom^2'), symmetry=1, barrier=(4.68436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96379,0.0495426,-4.78915e-05,2.70935e-08,-6.76482e-12,-64745.1,20.1411], Tmin=(100,'K'), Tmax=(924.146,'K')), NASAPolynomial(coeffs=[6.74962,0.0288284,-1.42704e-05,2.84009e-09,-2.03922e-13,-65629.7,-2.57026], Tmin=(924.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cd-CdCs(CO)) + group(COCFO) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(6730)',
    structure = adjacencyList("""1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.00704622,2.78305e-05,-3.99443e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.706,'K')), NASAPolynomial(coeffs=[6.79957,0.00959978,-3.02068e-06,5.37826e-10,-3.92605e-14,19772.3,-12.7583], Tmin=(949.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    E0 = (-411.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (13.8234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (207.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-55.2372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-30.2239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (136.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-123.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-200.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-183.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (4.01664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-70.5457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-284.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-126.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-33.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-40.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-116.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-11.5359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-138.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-243.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-179.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(CF)C(F)(F)OF(9176)'],
    products = ['HF(38)', 'CF2O(49)', 'C=C=CF(6101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(59.0397,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'C=C(CF)OF(581)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(214.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FOF(499)', 'C=C([C]F)CF(9182)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(13.5911,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2(8)', 'C=C([C]F)C(F)(F)OF(9183)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(85.3813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHF(40)', 'C=CC(F)(F)OF(4807)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.44768e+28,'m^3/(mol*s)'), n=-6.4458, Ea=(80.2034,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.6591608693425184, var=5.4995193120720405, Tref=1000.0, N=19, data_mean=0.0, correlation='CH',), comment="""Estimated from node CH"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(S)(25)', 'C=C(F)C(F)(F)OF(4967)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(144.315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FCC(CF)=C(F)OF(9184)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(272.912,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FCC(COF)=C(F)F(9185)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(195.358,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH]F)C(F)(F)OF(9186)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]C(CF)=C(F)OF(9187)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0.519206,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2]C(=C)C(F)(F)OF(9188)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.07724e+07,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C=C(CF)C([O])(F)F(9189)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]F(126)', '[CH2]C(CF)=C(F)F(9190)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['FO[C](F)F(4688)', 'C=[C]CF(6104)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -27.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2F(46)', 'C=[C]C(F)(F)OF(4958)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -25.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH2]C(=CF)C(F)(F)OF(9191)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(7.22644,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[CH]=C(CF)C(F)(F)OF(9192)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -5.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(CF)C(F)(F)OF(9193)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(14.721,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F2(78)', 'C=C(CF)C(=O)F(7956)'],
    products = ['C=C(CF)C(F)(F)OF(9176)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(80.4712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(CF)C(F)(F)OF(9176)'],
    products = ['F2(78)', 'CF2O(49)', 'C=C=C(6730)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(290.583,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

network(
    label = 'PDepNetwork #3137',
    isomers = [
        'C=C(CF)C(F)(F)OF(9176)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(49)', 'C=C=CF(6101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3137',
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

