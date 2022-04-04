species(
    label = '[CH2]C(=CF)C1[CH]C1F(10532)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {8,D}
6  C u1 p0 c0 {3,S} {4,S} {11,S}
7  C u1 p0 c0 {5,S} {13,S} {14,S}
8  C u0 p0 c0 {2,S} {5,D} {12,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (92.2345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180,882.621,882.633,882.641,882.644,882.648,882.65],'cm^-1')),
        HinderedRotor(inertia=(0.0864859,'amu*angstrom^2'), symmetry=1, barrier=(1.98848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00359598,'amu*angstrom^2'), symmetry=1, barrier=(1.98792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06915,0.0553723,-2.9532e-05,1.92061e-09,2.18171e-12,11206.3,27.6765], Tmin=(100,'K'), Tmax=(1160.75,'K')), NASAPolynomial(coeffs=[13.3929,0.0271702,-1.15227e-05,2.16528e-09,-1.51427e-13,7384.29,-37.7553], Tmin=(1160.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.2345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Allyl_P)"""),
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
    label = 'FC1C=C1(6880)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {2,S} {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (53.2336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,725.22,725.471,725.523,725.698],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2885.09,'J/mol'), sigma=(4.85707,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=450.64 K, Pc=57.13 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.11981,-0.0104554,0.000119533,-1.9991e-07,1.08069e-10,6403.99,8.0387], Tmin=(10,'K'), Tmax=(590.766,'K')), NASAPolynomial(coeffs=[1.59625,0.0247584,-1.59029e-05,4.86605e-09,-5.67736e-13,6385.84,16.208], Tmin=(590.766,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(53.2336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1CDC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=[C]CC1[CH]C1F(10534)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (210.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,615,860,1140,1343,3152,1685,370,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3378.06,'J/mol'), sigma=(5.68685,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.65 K, Pc=41.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07295,0.0539072,-2.18765e-05,-8.94911e-09,6.65143e-12,25466.1,29.5458], Tmin=(100,'K'), Tmax=(1066.67,'K')), NASAPolynomial(coeffs=[14.1651,0.0253314,-1.05475e-05,2.00502e-09,-1.42841e-13,21505.8,-39.933], Tmin=(1066.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[CH2]C1=C[CH]C(F)C1F(10582)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {7,D} {8,S}
6  C u1 p0 c0 {4,S} {7,S} {11,S}
7  C u0 p0 c0 {5,D} {6,S} {12,S}
8  C u1 p0 c0 {5,S} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-125.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29002,0.0449571,9.99594e-06,-4.3839e-08,1.9058e-11,-14998.3,19.8283], Tmin=(100,'K'), Tmax=(1020.46,'K')), NASAPolynomial(coeffs=[15.0062,0.0250474,-1.05019e-05,2.06251e-09,-1.51809e-13,-19560.4,-55.2587], Tmin=(1020.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-allyl) + radical(Allyl_P) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC=[C]C1[CH]C1F(9373)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {2,S} {7,D} {11,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (217.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,615,860,1140,1343,3152,1685,370,180,180,491.667,1351.75,1352.35,1352.59,1352.83,1352.85],'cm^-1')),
        HinderedRotor(inertia=(0.183166,'amu*angstrom^2'), symmetry=1, barrier=(4.21135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0669,0.0402377,-2.43736e-05,6.90871e-09,-7.79827e-13,26246.2,23.5328], Tmin=(100,'K'), Tmax=(1993.53,'K')), NASAPolynomial(coeffs=[12.8579,0.0185856,-8.0819e-06,1.46051e-09,-9.65909e-14,21943.8,-35.9723], Tmin=(1993.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cds_S)"""),
)

species(
    label = 'FC=C1CC2C(F)C12(10541)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
6  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {2,S} {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-169.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58281,0.036202,3.34302e-05,-6.84352e-08,2.8006e-11,-20333.8,19.2831], Tmin=(100,'K'), Tmax=(993.9,'K')), NASAPolynomial(coeffs=[14.9892,0.0233665,-9.25599e-06,1.82263e-09,-1.36544e-13,-25029.7,-55.5304], Tmin=(993.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(CdCFH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'CC(=CF)C1=CC1F(10583)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
4  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {3,S} {5,D} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-112.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491684,0.068022,-5.55684e-05,2.0217e-08,-2.12487e-12,-13357.7,23.0124], Tmin=(100,'K'), Tmax=(1105.19,'K')), NASAPolynomial(coeffs=[16.6527,0.0214319,-8.48786e-06,1.56095e-09,-1.08823e-13,-17656.8,-59.8602], Tmin=(1105.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cd-Cs(F)-Cd)"""),
)

species(
    label = 'CC(=CF)C1C=C1F(10540)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-48.6749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26888,0.0633214,-5.74519e-05,2.96385e-08,-6.55431e-12,-5758.6,21.8565], Tmin=(100,'K'), Tmax=(1051.56,'K')), NASAPolynomial(coeffs=[9.03022,0.0337982,-1.53383e-05,2.93932e-09,-2.06764e-13,-7390.9,-15.9777], Tmin=(1051.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.6749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'FC1C[C]1C1[CH]C1F(10584)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u1 p0 c0 {3,S} {6,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (194.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52509,0.0461832,-1.23666e-05,-9.18745e-09,4.58429e-12,23547.1,26.7031], Tmin=(100,'K'), Tmax=(1197.55,'K')), NASAPolynomial(coeffs=[10.7221,0.0313252,-1.36233e-05,2.57198e-09,-1.79439e-13,20207,-24.0738], Tmin=(1197.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsCsHH) + ring(Cs(C)-Cs-Cs(F)) + ring(Cs(F)-Cs-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = '[CH2][C]1C(F)C2C(F)C12(10585)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
7  C u1 p0 c0 {4,S} {6,S} {8,S}
8  C u1 p0 c0 {7,S} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (76.3452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40157,0.0527378,-3.21561e-05,9.45898e-09,-1.12105e-12,9279.07,21.2455], Tmin=(100,'K'), Tmax=(1899.87,'K')), NASAPolynomial(coeffs=[14.3903,0.0253909,-1.05647e-05,1.88247e-09,-1.24066e-13,4343.75,-49.7534], Tmin=(1899.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.3452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]F)=CC=CF(8370)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,D} {6,S} {7,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {8,D} {9,S}
6  C u1 p0 c0 {1,S} {3,S} {11,S}
7  C u1 p0 c0 {3,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-84.6254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,234,589,736,816,1240,3237,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,302.855],'cm^-1')),
        HinderedRotor(inertia=(0.00185525,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767982,'amu*angstrom^2'), symmetry=1, barrier=(49.225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762915,'amu*angstrom^2'), symmetry=1, barrier=(49.2235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711089,0.0624412,-3.77056e-05,1.39698e-09,4.67836e-12,-10051,27.576], Tmin=(100,'K'), Tmax=(1013.83,'K')), NASAPolynomial(coeffs=[15.6846,0.0231252,-8.77347e-06,1.59776e-09,-1.12021e-13,-14102.7,-49.877], Tmin=(1013.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.6254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[CH2]C1([CH]F)C2C(F)C21(10586)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
7  C u1 p0 c0 {3,S} {13,S} {14,S}
8  C u1 p0 c0 {2,S} {3,S} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (126.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605224,0.0641699,-4.50451e-05,1.17304e-08,-1.33949e-13,15289.5,23.7119], Tmin=(100,'K'), Tmax=(1190.01,'K')), NASAPolynomial(coeffs=[17.0966,0.0226058,-1.01349e-05,1.9664e-09,-1.401e-13,10382.5,-62.8442], Tmin=(1190.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsHHH) + group(CsCsFHH) + polycyclic(s2_3_3_ane) + radical(Neopentyl) + radical(CsCsF1sH)"""),
)

species(
    label = 'C=[C][CH]F(6102)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (196.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.900552,'amu*angstrom^2'), symmetry=1, barrier=(20.7055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79303e-05,1.96928e-08,-6.1986e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.867,'K')), NASAPolynomial(coeffs=[4.81186,0.016118,-7.47822e-06,1.46101e-09,-1.03887e-13,23340.4,3.53847], Tmin=(747.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = '[CH2]C(=CF)C1=CC1F(10587)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u0 p0 c0 {3,S} {4,D} {10,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {11,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (39.3223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,1000,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,468.817,469.053,469.06,469.08,469.145,469.541],'cm^-1')),
        HinderedRotor(inertia=(0.193496,'amu*angstrom^2'), symmetry=1, barrier=(30.1365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192782,'amu*angstrom^2'), symmetry=1, barrier=(30.1437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564873,0.0647114,-4.44143e-05,3.88649e-09,4.96093e-12,4862.67,23.0271], Tmin=(100,'K'), Tmax=(996.787,'K')), NASAPolynomial(coeffs=[18.218,0.0165767,-6.14737e-06,1.14488e-09,-8.28589e-14,215.409,-67.7401], Tmin=(996.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.3223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cd-Cs(F)-Cd) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=CF)C1C=C1(10588)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {7,D}
4  C u0 p0 c0 {2,S} {5,D} {9,S}
5  C u0 p0 c0 {2,S} {4,D} {10,S}
6  C u1 p0 c0 {3,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {3,D} {11,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (271.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,350,440,435,1725,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.920274,0.0591942,-4.90765e-05,2.09372e-08,-3.56311e-12,32773.2,21.8371], Tmin=(100,'K'), Tmax=(1412.05,'K')), NASAPolynomial(coeffs=[14.742,0.0200402,-7.48336e-06,1.29976e-09,-8.63177e-14,28869.8,-49.6137], Tmin=(1412.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cyclopropene) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CF)C1C=C1F(10589)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {8,D}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  C u1 p0 c0 {4,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {4,D} {11,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (102.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,350,440,435,1725,323,467,575,827,1418,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180,180,1007.07,1009.39,1011.97,1012.22],'cm^-1')),
        HinderedRotor(inertia=(0.00628963,'amu*angstrom^2'), symmetry=1, barrier=(4.54319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198489,'amu*angstrom^2'), symmetry=1, barrier=(4.56366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08395,0.0630993,-5.7194e-05,2.73638e-08,-5.35563e-12,12472.7,22.7929], Tmin=(100,'K'), Tmax=(1211.13,'K')), NASAPolynomial(coeffs=[12.2857,0.0261031,-1.13735e-05,2.1418e-09,-1.49298e-13,9759.37,-33.3947], Tmin=(1211.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cd(F)-Cd) + radical(Allyl_P)"""),
)

species(
    label = 'FC1[CH][CH]1(8666)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {4,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (308.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,3150,900,1100,1800.48,1800.48,1800.48,1800.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60776,0.0050847,2.30871e-05,-2.54932e-08,7.79552e-12,37164.6,14.3884], Tmin=(100,'K'), Tmax=(1144.07,'K')), NASAPolynomial(coeffs=[3.87852,0.0162024,-7.30699e-06,1.43499e-09,-1.02849e-13,36313.1,9.595], Tmin=(1144.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(cyclopropane)"""),
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
    label = '[CH2]C(=CF)[C]1C=C1(10590)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  C u1 p0 c0 {3,S} {4,S} {5,S}
3  C u0 p0 c0 {2,S} {6,S} {7,D}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (403.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,3150,900,1100,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180,784.632,784.674,784.731,784.761,784.813,784.884,784.941],'cm^-1')),
        HinderedRotor(inertia=(0.00723615,'amu*angstrom^2'), symmetry=1, barrier=(3.16174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13757,'amu*angstrom^2'), symmetry=1, barrier=(3.16301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19418,0.0538617,-4.39298e-05,1.87195e-08,-3.19881e-12,48636.4,19.106], Tmin=(100,'K'), Tmax=(1403.87,'K')), NASAPolynomial(coeffs=[13.3139,0.019329,-7.0322e-06,1.19752e-09,-7.84903e-14,45233.5,-43.4762], Tmin=(1403.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cyclopropene) + radical(Allyl_T) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CF)[C]1CC1F(8463)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {6,S} {13,S} {14,S}
8  C u0 p0 c0 {2,S} {6,D} {12,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1.6975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23506,0.0488883,-4.67296e-06,-2.77272e-08,1.35968e-11,-94.1745,23.6724], Tmin=(100,'K'), Tmax=(1012.55,'K')), NASAPolynomial(coeffs=[13.7669,0.0260166,-1.02468e-05,1.92064e-09,-1.37215e-13,-3997.36,-43.685], Tmin=(1012.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.6975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(Allyl_T) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CF)C1C[C]1F(10591)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {8,D}
6  C u1 p0 c0 {1,S} {3,S} {4,S}
7  C u1 p0 c0 {5,S} {13,S} {14,S}
8  C u0 p0 c0 {2,S} {5,D} {12,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (56.7958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824766,0.0600611,-4.28396e-05,1.51218e-08,-2.13488e-12,6953.45,28.9261], Tmin=(100,'K'), Tmax=(1664.35,'K')), NASAPolynomial(coeffs=[16.3743,0.0226908,-9.16001e-06,1.63141e-09,-1.0854e-13,1777.4,-54.0131], Tmin=(1664.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.7958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(Allyl_P)"""),
)

species(
    label = 'CC([CH]F)=C1[CH]C1F(10592)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
4  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {8,S}
7  C u1 p0 c0 {3,S} {5,S} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (37.3055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950463,0.0592826,-4.05176e-05,1.33905e-08,-1.76924e-12,4602.96,25.6574], Tmin=(100,'K'), Tmax=(1754.51,'K')), NASAPolynomial(coeffs=[16.5487,0.0237213,-1.01149e-05,1.83835e-09,-1.23185e-13,-870.513,-58.3642], Tmin=(1754.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.3055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cd(Cd)-Cs-Cs(F)) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'CC(=[C]F)C1[CH]C1F(10593)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
5  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {4,S} {14,S}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (191.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,167,640,1190,229.122,229.993,1292.22,1292.28,1292.28,1292.28,1292.34],'cm^-1')),
        HinderedRotor(inertia=(0.154044,'amu*angstrom^2'), symmetry=1, barrier=(5.7726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153631,'amu*angstrom^2'), symmetry=1, barrier=(5.77454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51488,0.0561794,-3.91718e-05,1.39637e-08,-2.09957e-12,23081.6,26.2685], Tmin=(100,'K'), Tmax=(1471.45,'K')), NASAPolynomial(coeffs=[10.4771,0.0318166,-1.43365e-05,2.71173e-09,-1.87859e-13,20444.1,-20.4307], Tmin=(1471.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'CC(=CF)C1[CH][C]1F(10594)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {1,S} {3,S} {7,S}
7  C u1 p0 c0 {3,S} {6,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (131.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,212,367,445,1450,194,682,905,1196,1383,3221,209.929,210.643,1485.99,1486,1486.04,1486.05,1486.08],'cm^-1')),
        HinderedRotor(inertia=(0.0909486,'amu*angstrom^2'), symmetry=1, barrier=(2.8792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916769,'amu*angstrom^2'), symmetry=1, barrier=(2.87891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8752,0.0530606,-3.60177e-05,1.33996e-08,-2.30474e-12,15851.9,27.1995], Tmin=(100,'K'), Tmax=(1198.49,'K')), NASAPolynomial(coeffs=[6.11291,0.0389172,-1.83162e-05,3.55312e-09,-2.50806e-13,14836.1,5.9877], Tmin=(1198.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(cyclopropane)"""),
)

species(
    label = 'C=C([C]F)C1CC1F(10595)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  C u2 p0 c0 {2,S} {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (101.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,250,417,511,1155,1315,1456,3119,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940137,0.0588804,-4.07476e-05,1.35647e-08,-1.79712e-12,12347.2,25.9408], Tmin=(100,'K'), Tmax=(1759.19,'K')), NASAPolynomial(coeffs=[17.0338,0.0222875,-9.54656e-06,1.7408e-09,-1.16839e-13,6684.8,-60.7927], Tmin=(1759.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'FC=C(CF)C1[CH][CH]1(10596)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {3,S} {7,S} {12,S}
7  C u1 p0 c0 {3,S} {6,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (159.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,194,682,905,1196,1383,3221,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37464,0.0536106,-3.4188e-05,1.04691e-08,-1.28915e-12,19317.2,29.0215], Tmin=(100,'K'), Tmax=(1839.44,'K')), NASAPolynomial(coeffs=[14.6509,0.0247404,-1.06453e-05,1.9365e-09,-1.29478e-13,14433.1,-43.1202], Tmin=(1839.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cyclopropane) + radical(cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = '[CH]C(=C)C1C(F)C1F(10597)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  C u2 p0 c0 {6,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (90.0611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,194,306,365,469,408,614,1107,1203,1236,1394,1392,1520,3058,3180,350,440,435,1725,2950,3100,1380,975,1025,1650,417.667,417.681,417.72,417.72,417.735,417.735,417.747,417.752],'cm^-1')),
        HinderedRotor(inertia=(0.000966244,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000966123,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850272,0.0587014,-2.27856e-05,-8.63029e-09,6.30637e-12,10954.2,26.5967], Tmin=(100,'K'), Tmax=(1083.32,'K')), NASAPolynomial(coeffs=[13.5161,0.031863,-1.32175e-05,2.46248e-09,-1.72263e-13,7040.61,-40.9193], Tmin=(1083.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.0611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(AllylJ2_triplet) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]=C(CF)C1[CH]C1F(10598)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {4,S} {13,S}
8  C u1 p0 c0 {6,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (204.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3120,650,792.5,1650,237.771,237.771,237.773,1387.59,1387.59,1387.59,1387.59],'cm^-1')),
        HinderedRotor(inertia=(0.0049319,'amu*angstrom^2'), symmetry=1, barrier=(6.73851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684115,'amu*angstrom^2'), symmetry=1, barrier=(27.4457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45918,0.0559352,-3.7449e-05,1.22062e-08,-1.63169e-12,24640.8,27.1783], Tmin=(100,'K'), Tmax=(1663.62,'K')), NASAPolynomial(coeffs=[12.6522,0.0290228,-1.31838e-05,2.48238e-09,-1.7046e-13,20916.6,-32.5191], Tmin=(1663.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(F)C1[CH]C1F(10533)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {12,S}
7  C u0 p0 c0 {8,D} {13,S} {14,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (149.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,2950,3100,1380,975,1025,1650,1685,370,182.432,182.433,182.454,182.481,182.536,1476.17,1476.17],'cm^-1')),
        HinderedRotor(inertia=(0.000387505,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0153282,'amu*angstrom^2'), symmetry=1, barrier=(23.702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3318.08,'J/mol'), sigma=(5.64812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.28 K, Pc=41.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35584,0.0625565,-5.13572e-05,2.25151e-08,-4.23792e-12,18018.1,24.8012], Tmin=(100,'K'), Tmax=(1201.04,'K')), NASAPolynomial(coeffs=[9.69,0.0348001,-1.6692e-05,3.27347e-09,-2.32731e-13,16016.2,-16.9331], Tmin=(1201.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(C-F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-HH))"""),
)

species(
    label = 'F[CH]C1=C[CH]C(F)C1(10599)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {7,D} {8,S}
6  C u1 p0 c0 {4,S} {7,S} {12,S}
7  C u0 p0 c0 {5,D} {6,S} {13,S}
8  C u1 p0 c0 {2,S} {5,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-131.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46296,0.0399937,2.33397e-05,-5.70115e-08,2.35563e-11,-15742.8,21.429], Tmin=(100,'K'), Tmax=(1010.72,'K')), NASAPolynomial(coeffs=[14.6957,0.0251964,-1.04605e-05,2.06216e-09,-1.52777e-13,-20336.8,-52.0464], Tmin=(1010.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-allyl) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[CH]F(230)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'C=[C]C1[CH]C1F(8780)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {2,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (409.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2950,3100,1380,975,1025,1650,1685,370,180,180,180,1233.24,1233.28,1233.59,1233.59],'cm^-1')),
        HinderedRotor(inertia=(0.00285399,'amu*angstrom^2'), symmetry=1, barrier=(3.07976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34261,0.0305819,-3.16973e-06,-1.09787e-08,4.47581e-12,49306.7,21.688], Tmin=(100,'K'), Tmax=(1183.47,'K')), NASAPolynomial(coeffs=[8.05562,0.0237444,-1.03112e-05,1.94899e-09,-1.3616e-13,47081,-10.5263], Tmin=(1183.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane) + radical(Cds_S)"""),
)

species(
    label = 'C=C1C(F)C2C(F)C12(10561)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-178.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47083,0.0412808,1.48327e-05,-4.54313e-08,1.87212e-11,-21335.4,17.1312], Tmin=(100,'K'), Tmax=(1039.93,'K')), NASAPolynomial(coeffs=[14.0877,0.0261271,-1.14515e-05,2.28084e-09,-1.68151e-13,-25764.3,-52.9094], Tmin=(1039.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCCFH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'C=C(CF)C1=CC1F(10600)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {3,S} {5,D} {12,S}
8  C u0 p0 c0 {6,D} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-95.8863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.384607,0.068934,-5.90873e-05,2.50464e-08,-4.17996e-12,-11393.5,24.2024], Tmin=(100,'K'), Tmax=(1446.65,'K')), NASAPolynomial(coeffs=[18.2689,0.0194841,-7.81392e-06,1.41799e-09,-9.66696e-14,-16567.9,-68.6825], Tmin=(1446.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.8863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cd-Cs(F)-Cd)"""),
)

species(
    label = 'C=C(CF)C1C=C1F(10539)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {12,S}
8  C u0 p0 c0 {5,D} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-32.3844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37537,0.0618923,-5.36052e-05,2.60591e-08,-5.49396e-12,-3803.93,22.2685], Tmin=(100,'K'), Tmax=(1085.88,'K')), NASAPolynomial(coeffs=[8.7014,0.0349058,-1.63268e-05,3.17239e-09,-2.24783e-13,-5394.97,-13.679], Tmin=(1085.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.3844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'F[CH][C]1CC2C(F)C12(10601)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
6  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
7  C u1 p0 c0 {4,S} {6,S} {8,S}
8  C u1 p0 c0 {2,S} {7,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (87.6117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61809,0.050394,-2.75393e-05,6.39963e-09,-5.34668e-13,10624.7,21.1637], Tmin=(100,'K'), Tmax=(2020,'K')), NASAPolynomial(coeffs=[19.2217,0.0207424,-9.38732e-06,1.68497e-09,-1.09103e-13,2450.54,-78.7701], Tmin=(2020,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.6117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFHH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = '[CH2]C(CF)=C1[CH]C1F(10602)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {8,S}
7  C u1 p0 c0 {3,S} {5,S} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (43.0352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.760703,0.0617967,-4.43318e-05,1.54286e-08,-2.13723e-12,5300.41,25.4701], Tmin=(100,'K'), Tmax=(1694.23,'K')), NASAPolynomial(coeffs=[17.4071,0.0224951,-9.53576e-06,1.73662e-09,-1.16834e-13,-340.153,-63.6156], Tmin=(1694.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cd(Cd)-Cs-Cs(F)) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(CF)C1[CH][C]1F(10603)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {2,S} {3,S} {7,S}
7  C u1 p0 c0 {3,S} {6,S} {12,S}
8  C u0 p0 c0 {5,D} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (147.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,212,367,445,1450,2950,3100,1380,975,1025,1650,180,180,1220.97,1221.53,1221.56,1222.03,1222.07,1222.33],'cm^-1')),
        HinderedRotor(inertia=(0.151168,'amu*angstrom^2'), symmetry=1, barrier=(3.47565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150344,'amu*angstrom^2'), symmetry=1, barrier=(3.45671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1424,0.0502934,-2.93943e-05,7.93486e-09,-8.81391e-13,17797.1,26.9874], Tmin=(100,'K'), Tmax=(1817.28,'K')), NASAPolynomial(coeffs=[9.47523,0.0341532,-1.60721e-05,3.04761e-09,-2.09062e-13,15131.9,-12.7695], Tmin=(1817.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(cyclopropane)"""),
)

species(
    label = '[CH]C(=CF)C1CC1F(10604)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {6,D} {13,S}
8  C u2 p0 c0 {6,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (85.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924145,0.0570414,-1.91952e-05,-1.14889e-08,7.11786e-12,10403.8,26.6735], Tmin=(100,'K'), Tmax=(1080.86,'K')), NASAPolynomial(coeffs=[13.1112,0.032425,-1.34618e-05,2.50942e-09,-1.75615e-13,6572.75,-38.605], Tmin=(1080.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(C(F)F)C1[CH][CH]1(10605)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {8,D}
6  C u1 p0 c0 {3,S} {7,S} {11,S}
7  C u1 p0 c0 {3,S} {6,S} {12,S}
8  C u0 p0 c0 {5,D} {13,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (133.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24909,0.0536073,-3.37466e-05,1.0053e-08,-1.18817e-12,16116.3,29.7507], Tmin=(100,'K'), Tmax=(1946.08,'K')), NASAPolynomial(coeffs=[16.7362,0.0217747,-9.21066e-06,1.64774e-09,-1.08401e-13,10088.5,-55.2775], Tmin=(1946.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(cyclopropane) + radical(cyclopropane)"""),
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
    E0 = (-30.2422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (258.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (42.5434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (476.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-21.9579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (33.1579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-5.26896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (132.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (94.9989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (24.1141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (34.4625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (129.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (133.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (266.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (192.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (195.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (382.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (246.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (131.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (131.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (172.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (224.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (91.1575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (23.5188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (247.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (196.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (257.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (121.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (40.4633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (501.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-21.9579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (33.1579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-5.26896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (94.9989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (425.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (178.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (130.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (53.0773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (232.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['C=C=CF(6101)', 'FC1C=C1(6880)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]CC1[CH]C1F(10534)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2]C1=C[CH]C(F)C1F(10582)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(72.7856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(18)', 'FC=[C]C1[CH]C1F(9373)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['FC=C1CC2C(F)C12(10541)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['CC(=CF)C1=CC1F(10583)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['CC(=CF)C1C=C1F(10540)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['FC1C[C]1C1[CH]C1F(10584)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.57691e+12,'s^-1'), n=0.0709135, Ea=(162.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3_D;doublebond_intra;radadd_intra_cs2H] + [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2][C]1C(F)C2C(F)C12(10585)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH]F)=CC=CF(8370)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2]C1([CH]F)C2C(F)C21(10586)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C][CH]F(6102)', 'FC1C=C1(6880)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(2.52929,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[CH2]C(=CF)C1=CC1F(10587)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(5.08432,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[CH2]C(=CF)C1C=C1(10588)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(44.3293,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH2]C(=CF)C1C=C1F(10589)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.596503,'m^3/(mol*s)'), n=2.29663, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08546600576587018, var=0.6163705976545122, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C=CF(6101)', 'FC1[CH][CH]1(8666)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00100541,'m^3/(mol*s)'), n=2.87982, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][CH]F(6102)', 'FC1[CH][CH]1(8666)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[CH2]C(=CF)[C]1C=C1(10590)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(247.032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2]C(=CF)[C]1CC1F(8463)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2]C(=CF)C1C[C]1F(10591)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['CC([CH]F)=C1[CH]C1F(10592)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(=[C]F)C1[CH]C1F(10593)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC(=CF)C1[CH][C]1F(10594)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(32400,'s^-1'), n=2.04, Ea=(82.4248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([C]F)C1CC1F(10595)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC=C(CF)C1[CH][CH]1(10596)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(209.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(=C)C1C(F)C1F(10597)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.12998e-06,'s^-1'), n=5.16802, Ea=(228.655,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(CF)C1[CH]C1F(10598)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(175.382,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(F)C1[CH]C1F(10533)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['F[CH]C1=C[CH]C(F)C1(10599)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(70.7055,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]F(230)', 'C=[C]C1[CH]C1F(8780)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['C=C1C(F)C2C(F)C12(10561)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['C=C(CF)C1=CC1F(10600)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['C=C(CF)C1C=C1F(10539)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['F[CH][C]1CC2C(F)C12(10601)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CHF(40)', 'C=[C]C1[CH]C1F(8780)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH2]C(CF)=C1[CH]C1F(10602)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C(CF)C1[CH][C]1F(10603)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    products = ['[CH]C(=CF)C1CC1F(10604)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(41.5432,'s^-1'), n=2.90518, Ea=(83.3195,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 2.8284271247461903
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C(C(F)F)C1[CH][CH]1(10605)'],
    products = ['[CH2]C(=CF)C1[CH]C1F(10532)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.00220343,'s^-1'), n=4.50663, Ea=(222.302,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

network(
    label = 'PDepNetwork #3577',
    isomers = [
        '[CH2]C(=CF)C1[CH]C1F(10532)',
    ],
    reactants = [
        ('C=C=CF(6101)', 'FC1C=C1(6880)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3577',
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

