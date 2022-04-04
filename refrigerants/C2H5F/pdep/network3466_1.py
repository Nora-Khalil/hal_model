species(
    label = '[CH2]C(F)(F)C1[CH]C1F(10035)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
7  C u1 p0 c0 {4,S} {5,S} {11,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-230.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,180,180,180,736.741,737.44,738.314,739.086],'cm^-1')),
        HinderedRotor(inertia=(0.00591965,'amu*angstrom^2'), symmetry=1, barrier=(2.28737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00591886,'amu*angstrom^2'), symmetry=1, barrier=(2.28034,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02775,0.0651388,-5.76351e-05,2.57605e-08,-4.69066e-12,-27631.6,23.7468], Tmin=(100,'K'), Tmax=(1290.23,'K')), NASAPolynomial(coeffs=[13.5442,0.0263354,-1.2523e-05,2.45103e-09,-1.74138e-13,-30861.4,-39.8271], Tmin=(1290.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs-Cs(C-FF)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CH2CF2(56)',
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
    label = 'FC1C=C1(7087)',
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
    label = 'F[C](F)CC1[CH]C1F(10036)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {13,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-157.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180,899.459,899.544,899.722,899.755,899.818,899.906],'cm^-1')),
        HinderedRotor(inertia=(0.173192,'amu*angstrom^2'), symmetry=1, barrier=(3.98203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173181,'amu*angstrom^2'), symmetry=1, barrier=(3.98177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3163.46,'J/mol'), sigma=(5.42643,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=494.12 K, Pc=44.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25745,0.0593097,-4.76195e-05,1.94276e-08,-3.25576e-12,-18898.6,26.429], Tmin=(100,'K'), Tmax=(1382,'K')), NASAPolynomial(coeffs=[12.5385,0.0266585,-1.21805e-05,2.33218e-09,-1.63263e-13,-22016.7,-31.6455], Tmin=(1382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C1[CH]C1F(10100)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u1 p0 c0 {4,S} {5,S} {10,S}
7  C u1 p0 c0 {2,S} {3,S} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-159.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,180,695.202,695.218,695.222,695.226,695.228,695.263],'cm^-1')),
        HinderedRotor(inertia=(0.00859749,'amu*angstrom^2'), symmetry=1, barrier=(2.94944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53496,0.0469392,-3.54121e-05,1.04421e-08,-3.56677e-13,-19081.1,22.1408], Tmin=(100,'K'), Tmax=(1130.08,'K')), NASAPolynomial(coeffs=[13.5855,0.0145549,-6.05805e-06,1.16671e-09,-8.36921e-14,-22460.5,-40.371], Tmin=(1130.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Csj(Cs-CsCsH)(F1s)(F1s)_1959_ring)"""),
)

species(
    label = 'FC1C2CC(F)(F)C12(10086)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {4,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
7  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-513.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32073,0.0474566,-1.01548e-05,-1.85865e-08,9.27793e-12,-61629.4,17.0973], Tmin=(100,'K'), Tmax=(1087.54,'K')), NASAPolynomial(coeffs=[14.2254,0.023996,-1.09031e-05,2.16656e-09,-1.57896e-13,-65855.8,-52.7696], Tmin=(1087.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-513.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'CC(F)(F)C1=CC1F(10101)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-422.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29231,0.0621942,-5.47378e-05,2.49154e-08,-4.71939e-12,-50693.9,20.6923], Tmin=(100,'K'), Tmax=(1223.65,'K')), NASAPolynomial(coeffs=[11.4263,0.0290674,-1.41299e-05,2.7916e-09,-1.9938e-13,-53174,-30.2437], Tmin=(1223.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd-Cs(F)-Cd)"""),
)

species(
    label = 'CC(F)(F)C1C=C1F(10038)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-426.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18512,0.0615108,-5.26987e-05,2.32985e-08,-4.23205e-12,-51192.7,21.6749], Tmin=(100,'K'), Tmax=(1287.66,'K')), NASAPolynomial(coeffs=[12.3645,0.0267829,-1.2244e-05,2.35358e-09,-1.65574e-13,-54071.8,-35.0855], Tmin=(1287.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-426.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(C-FF))"""),
)

species(
    label = '[CH2]C(F)(F)[CH]C=CF(9223)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u1 p0 c0 {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {6,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-323.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180,2625.91],'cm^-1')),
        HinderedRotor(inertia=(0.0178721,'amu*angstrom^2'), symmetry=1, barrier=(9.93515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34533,'amu*angstrom^2'), symmetry=1, barrier=(30.9318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310242,'amu*angstrom^2'), symmetry=1, barrier=(30.9897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948448,0.0699033,-7.22044e-05,3.94689e-08,-8.85008e-12,-38759.4,24.3231], Tmin=(100,'K'), Tmax=(1062.75,'K')), NASAPolynomial(coeffs=[11.9964,0.0283204,-1.35127e-05,2.65125e-09,-1.89091e-13,-41107.7,-29.6494], Tmin=(1062.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Allyl_S) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2][C](F)F(4220)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 C u1 p0 c0 {3,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([190,488,555,1236,1407,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25819,0.0186259,-1.56412e-05,8.22293e-09,-2.04575e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.673,'K')), NASAPolynomial(coeffs=[4.4701,0.0131647,-6.41279e-06,1.29206e-09,-9.37475e-14,-13305.9,7.85286], Tmin=(887.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sF1s) + radical(Cs_P)"""),
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
    label = '[CH2]C(F)(F)C1=CC1F(10102)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {4,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-210.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,248,333,466,604,684,796,1061,1199,2950,1000,3000,3100,440,815,1455,1000,180,180,180,1095.89,1096.98],'cm^-1')),
        HinderedRotor(inertia=(0.206896,'amu*angstrom^2'), symmetry=1, barrier=(4.75695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.933898,'amu*angstrom^2'), symmetry=1, barrier=(21.4721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21229,0.0670382,-7.51706e-05,4.58548e-08,-1.16843e-11,-25270,22.2957], Tmin=(100,'K'), Tmax=(932.651,'K')), NASAPolynomial(coeffs=[9.88232,0.0298531,-1.5364e-05,3.10385e-09,-2.24562e-13,-26887.1,-18.9275], Tmin=(932.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd-Cs(F)-Cd) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1C=C1(10103)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-5.91975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22702,0.0541195,-4.31398e-05,1.71388e-08,-2.71884e-12,-606.526,23.0726], Tmin=(100,'K'), Tmax=(1497.03,'K')), NASAPolynomial(coeffs=[14.3846,0.0189626,-7.91293e-06,1.45125e-09,-9.90375e-14,-4545.96,-45.714], Tmin=(1497.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.91975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(C-FF)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C=C1F(10104)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-214.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,215,315,519,588,595,1205,1248,323,467,575,827,1418,3000,3100,440,815,1455,1000,180,180,180,458.34,1491.13,3348.01],'cm^-1')),
        HinderedRotor(inertia=(0.0216567,'amu*angstrom^2'), symmetry=1, barrier=(34.0383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160918,'amu*angstrom^2'), symmetry=1, barrier=(25.3149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1118,0.0639425,-6.32015e-05,3.22136e-08,-6.65467e-12,-25684.9,23.5769], Tmin=(100,'K'), Tmax=(1155.64,'K')), NASAPolynomial(coeffs=[12.7068,0.0238088,-1.11086e-05,2.16221e-09,-1.53627e-13,-28364.8,-34.0396], Tmin=(1155.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(C-FF)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'C=C(F)C1[CH]C1F(9919)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-36.9173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180,1286.82,1286.83,1286.84,1286.88,1286.93],'cm^-1')),
        HinderedRotor(inertia=(0.17084,'amu*angstrom^2'), symmetry=1, barrier=(3.92795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81561,0.0409054,-1.55558e-05,-4.07366e-09,2.98396e-12,-4355.64,23.2048], Tmin=(100,'K'), Tmax=(1194.8,'K')), NASAPolynomial(coeffs=[10.5485,0.024307,-1.05837e-05,2.005e-09,-1.40322e-13,-7344.51,-24.2556], Tmin=(1194.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.9173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane)"""),
)

species(
    label = 'FC1[CH][CH]1(9903)',
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
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,3150,900,1100,1800.38,1800.4,1800.51,1800.62],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60777,0.00508457,2.30875e-05,-2.54937e-08,7.79572e-12,37164.6,14.3883], Tmin=(100,'K'), Tmax=(1144.06,'K')), NASAPolynomial(coeffs=[3.87844,0.0162025,-7.30706e-06,1.43501e-09,-1.0285e-13,36313.1,9.59549], Tmin=(1144.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(cyclopropane)"""),
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
    label = '[CH2]C(F)(F)[C]1C=C1(10105)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {6,D} {8,S}
6  C u0 p0 c0 {4,S} {5,D} {9,S}
7  C u1 p0 c0 {3,S} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (126.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,3150,900,1100,3000,3100,440,815,1455,1000,180,1150.33,1150.39,1150.46,1150.49,1150.5,1150.68,1150.97],'cm^-1')),
        HinderedRotor(inertia=(2.49514,'amu*angstrom^2'), symmetry=1, barrier=(57.3682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284492,'amu*angstrom^2'), symmetry=1, barrier=(6.54104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50107,0.0487851,-3.79864e-05,1.49133e-08,-2.35173e-12,15256.7,20.341], Tmin=(100,'K'), Tmax=(1500.01,'K')), NASAPolynomial(coeffs=[12.9436,0.0182713,-7.47239e-06,1.35138e-09,-9.13992e-14,11823.9,-39.5023], Tmin=(1500.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(C-FF)) + radical(Allyl_T) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)=C1[CH]C1F(10106)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {3,S} {4,S} {9,S}
6  C u0 p0 c0 {2,S} {4,D} {7,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (40.6355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2950,1000,271,519,563,612,1379,3000,3100,440,815,1455,1000,722.212,725.776,726.61,730.432,733.8],'cm^-1')),
        HinderedRotor(inertia=(0.000315611,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.486,0.0477113,-3.19631e-05,8.47106e-09,-3.56326e-13,4984.08,20.1522], Tmin=(100,'K'), Tmax=(1251.9,'K')), NASAPolynomial(coeffs=[13.1475,0.0192884,-8.49615e-06,1.61311e-09,-1.12847e-13,1371.74,-41.4938], Tmin=(1251.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cd(Cd-F)-Cs-Cs) + radical(Allyl_S) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)[C]1CC1F(9245)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-233.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58959,0.0614646,-5.1677e-05,8.29939e-10,2.29174e-11,-27985.9,23.9476], Tmin=(100,'K'), Tmax=(517.731,'K')), NASAPolynomial(coeffs=[6.18328,0.0380637,-1.89064e-05,3.73682e-09,-2.66171e-13,-28623.6,3.2449], Tmin=(517.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C[C]1F(10107)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-189.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,215,315,519,588,595,1205,1248,212,367,445,1450,3000,3100,440,815,1455,1000,180,180,844.604,845.484,847.751,849.583,851.778,852.231],'cm^-1')),
        HinderedRotor(inertia=(0.00665028,'amu*angstrom^2'), symmetry=1, barrier=(3.41357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0068317,'amu*angstrom^2'), symmetry=1, barrier=(3.47983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752207,0.0641636,-5.42642e-05,2.27297e-08,-3.78431e-12,-22615.4,26.1368], Tmin=(100,'K'), Tmax=(1435.34,'K')), NASAPolynomial(coeffs=[16.3259,0.0207623,-8.90727e-06,1.66272e-09,-1.14938e-13,-27086.1,-54.6257], Tmin=(1435.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CC(F)(F)[C]1[CH]C1F(10108)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {8,S}
8  C u1 p0 c0 {5,S} {7,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-257.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04669,0.0718431,-0.000100845,9.01486e-08,-3.37981e-11,-30843.4,23.8695], Tmin=(100,'K'), Tmax=(750.239,'K')), NASAPolynomial(coeffs=[5.58346,0.0395483,-2.00683e-05,3.9678e-09,-2.81083e-13,-31296,4.80629], Tmin=(750.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs-Cs(C-FF)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'CC(F)(F)C1[CH][C]1F(10109)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-166.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,551.91,551.911,551.911,551.913,551.914,551.914,551.914],'cm^-1')),
        HinderedRotor(inertia=(0.000553429,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000553425,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651125,0.0647838,-5.51985e-05,2.32354e-08,-3.86575e-12,-19848.9,25.9671], Tmin=(100,'K'), Tmax=(1444.69,'K')), NASAPolynomial(coeffs=[17.0479,0.0193859,-8.06337e-06,1.48475e-09,-1.01922e-13,-24586.6,-59.1701], Tmin=(1444.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'FCC(F)(F)C1[CH][CH]1(10110)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
7  C u1 p0 c0 {4,S} {8,S} {12,S}
8  C u1 p0 c0 {4,S} {7,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-182.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,222,329,445,522,589,1214,1475,528,1116,1182,1331,1402,1494,3075,3110,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20324,0.0476799,-8.50841e-06,-6.27959e-08,6.4098e-11,-21900.4,23.832], Tmin=(100,'K'), Tmax=(447.12,'K')), NASAPolynomial(coeffs=[3.93084,0.0418455,-2.12114e-05,4.26933e-09,-3.08779e-13,-22151.1,15.8125], Tmin=(447.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + ring(Cs-Cs(C-FF)-Cs) + radical(cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = '[CH2][C](F)C1C(F)C1F(10111)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {5,S} {11,S}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u1 p0 c0 {7,S} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-185.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,194,306,365,469,408,614,1107,1203,1236,1394,1392,1520,3058,3180,212,367,445,1450,3000,3100,440,815,1455,1000,180,1066.69,1066.71,1066.91,1067.19],'cm^-1')),
        HinderedRotor(inertia=(0.0904441,'amu*angstrom^2'), symmetry=1, barrier=(2.07949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900207,'amu*angstrom^2'), symmetry=1, barrier=(2.06975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.396,0.0628405,-6.26054e-05,3.6571e-08,-9.33712e-12,-22221.3,25.5925], Tmin=(100,'K'), Tmax=(911.57,'K')), NASAPolynomial(coeffs=[7.67355,0.0352935,-1.72751e-05,3.41813e-09,-2.44586e-13,-23365.8,-4.1116], Tmin=(911.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC[C](F)C1[CH]C1F(10112)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {13,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-186.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,333.732,333.993,334.023,334.096,334.153,1055.59,1680.05,1680.07],'cm^-1')),
        HinderedRotor(inertia=(0.00150886,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479394,'amu*angstrom^2'), symmetry=1, barrier=(37.927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4676,0.0610918,-5.53133e-05,2.79797e-08,-6.16378e-12,-22380.6,24.3232], Tmin=(100,'K'), Tmax=(1039.27,'K')), NASAPolynomial(coeffs=[8.37769,0.0344958,-1.69267e-05,3.35565e-09,-2.40368e-13,-23816.9,-9.28015], Tmin=(1039.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(C-F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(CsCsCsF1s)"""),
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
    E0 = (-134.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (183.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (317.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-126.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-71.3498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-109.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-1.34087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (40.0675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (96.7857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (197.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (93.2685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (173.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (43.1291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (295.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (167.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (87.0891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (19.6401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (85.1224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-40.7455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (31.3747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (116.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (41.8928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['CH2CF2(56)', 'FC1C=C1(7087)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)CC1[CH]C1F(10036)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(18)', 'F[C](F)C1[CH]C1F(10100)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['FC1C2CC(F)(F)C12(10086)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['CC(F)(F)C1=CC1F(10101)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['CC(F)(F)C1C=C1F(10038)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)[CH]C=CF(9223)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C](F)F(4220)', 'FC1C=C1(7087)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.10266e-12,'m^3/(mol*s)'), n=5.04048, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5477765617574784, var=1.0081445727216138, Tref=1000.0, N=71, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH2]C(F)(F)C1=CC1F(10102)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]C(F)(F)C1C=C1(10103)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(34.7322,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[CH2]C(F)(F)C1C=C1F(10104)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C=C(F)C1[CH]C1F(9919)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(41.962,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CF2(56)', 'FC1[CH][CH]1(9903)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00100541,'m^3/(mol*s)'), n=2.87982, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](F)F(4220)', 'FC1[CH][CH]1(9903)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', '[CH2]C(F)(F)[C]1C=C1(10105)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(227.024,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', '[CH2]C(F)=C1[CH]C1F(10106)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(231.685,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['[CH2]C(F)(F)[C]1CC1F(9245)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2] for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(F)(F)C1C[C]1F(10107)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.802e+08,'s^-1'), n=1.453, Ea=(178.297,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S_cy3;C_rad_out_noH;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    products = ['CC(F)(F)[C]1[CH]C1F(10108)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4313e+06,'s^-1'), n=1.72764, Ea=(94.0044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC(F)(F)C1[CH][C]1F(10109)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.52029e+07,'s^-1'), n=1.31708, Ea=(101.589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FCC(F)(F)C1[CH][CH]1(10110)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(219.402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C](F)C1C(F)C1F(10111)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(205.848,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC[C](F)C1[CH]C1F(10112)'],
    products = ['[CH2]C(F)(F)C1[CH]C1F(10035)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(132.814,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #3466',
    isomers = [
        '[CH2]C(F)(F)C1[CH]C1F(10035)',
    ],
    reactants = [
        ('CH2CF2(56)', 'FC1C=C1(7087)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3466',
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

