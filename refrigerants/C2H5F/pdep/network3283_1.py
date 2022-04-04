species(
    label = '[CH2]C(F)(F)C[C]=CF(8554)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u0 p0 c0 {3,S} {8,D} {13,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-217.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,615,860,1140,1343,3152,1685,370,301.541,307.687,307.713],'cm^-1')),
        HinderedRotor(inertia=(0.126767,'amu*angstrom^2'), symmetry=1, barrier=(8.31901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00502292,'amu*angstrom^2'), symmetry=1, barrier=(8.31424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266005,'amu*angstrom^2'), symmetry=1, barrier=(17.9009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869218,0.0736591,-8.5759e-05,5.42682e-08,-1.41171e-11,-26055.1,27.6851], Tmin=(100,'K'), Tmax=(923.576,'K')), NASAPolynomial(coeffs=[10.992,0.0298179,-1.45564e-05,2.87254e-09,-2.05117e-13,-27925,-20.3467], Tmin=(923.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
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
    label = 'C=C=CF(5941)',
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
    label = 'FC=[C]CC[C](F)F(8557)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u0 p0 c0 {3,S} {8,D} {13,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-195.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,190,488,555,1236,1407,615,860,1140,1343,3152,1685,370,286.799,288.161,288.913,1531.89],'cm^-1')),
        HinderedRotor(inertia=(0.0020534,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156696,'amu*angstrom^2'), symmetry=1, barrier=(9.42348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162198,'amu*angstrom^2'), symmetry=1, barrier=(9.4197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3147.27,'J/mol'), sigma=(5.30768,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=491.60 K, Pc=47.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14003,0.0685016,-8.12231e-05,5.69191e-08,-1.69358e-11,-23374.3,28.5103], Tmin=(100,'K'), Tmax=(802.969,'K')), NASAPolynomial(coeffs=[8.08571,0.0339038,-1.65965e-05,3.26624e-09,-2.32355e-13,-24489.8,-3.47495], Tmin=(802.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[CH2]C(=CF)C([CH2])(F)F(8552)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {5,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-312.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693414,0.0783194,-9.67325e-05,6.51695e-08,-1.80024e-11,-37418.7,24.9261], Tmin=(100,'K'), Tmax=(873.424,'K')), NASAPolynomial(coeffs=[11.0426,0.0309215,-1.53292e-05,3.03354e-09,-2.16562e-13,-39226.5,-23.6018], Tmin=(873.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Csj(Cs-F1sF1sCd)(H)(H)) + radical(Allyl_P)"""),
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
    label = 'FC=[C]C[C](F)F(4639)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-178.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,615,860,1140,1343,3152,1685,370,291.713,298.097,1337.84],'cm^-1')),
        HinderedRotor(inertia=(0.09286,'amu*angstrom^2'), symmetry=1, barrier=(5.493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947839,'amu*angstrom^2'), symmetry=1, barrier=(5.51882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92501,0.0506398,-6.06767e-05,4.3892e-08,-1.36833e-11,-21352.6,23.6183], Tmin=(100,'K'), Tmax=(762.881,'K')), NASAPolynomial(coeffs=[6.49197,0.0266932,-1.35906e-05,2.74312e-09,-1.98193e-13,-22049.4,2.82153], Tmin=(762.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(F1s)) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[C]=CF(1334)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (404.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([194,682,905,1196,1383,3221],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.0070472,-1.17409e-06,-1.98501e-09,8.12281e-13,48709.9,8.54797], Tmin=(100,'K'), Tmax=(1284.59,'K')), NASAPolynomial(coeffs=[5.40185,0.00468,-2.11337e-06,4.24439e-10,-3.06832e-14,47991.3,-1.49755], Tmin=(1284.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])(F)F(4485)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 C u1 p0 c0 {3,S} {8,S} {9,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-147.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00588817,'amu*angstrom^2'), symmetry=1, barrier=(7.16366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556778,'amu*angstrom^2'), symmetry=1, barrier=(12.8014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0606,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09088,0.042004,-4.29653e-05,2.30839e-08,-4.98348e-12,-17634.5,17.7144], Tmin=(100,'K'), Tmax=(1116.45,'K')), NASAPolynomial(coeffs=[9.7005,0.0147404,-6.33535e-06,1.21104e-09,-8.56179e-14,-19333.6,-19.8359], Tmin=(1116.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'FC=C1CC(F)(F)C1(9203)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-520.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20764,0.0561294,-3.97281e-05,1.35258e-08,-1.84395e-12,-62439.6,20.8113], Tmin=(100,'K'), Tmax=(1695.74,'K')), NASAPolynomial(coeffs=[15.383,0.0226916,-1.015e-05,1.8974e-09,-1.29579e-13,-67247.1,-55.0631], Tmin=(1695.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-520.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'CC(F)(F)C=C=CF(9220)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {12,S}
7  C u0 p0 c0 {3,S} {8,D} {13,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-477.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16355,0.0681406,-7.63685e-05,4.92009e-08,-1.34632e-11,-57306.5,22.4195], Tmin=(100,'K'), Tmax=(867.457,'K')), NASAPolynomial(coeffs=[8.65353,0.0336039,-1.66496e-05,3.30643e-09,-2.36839e-13,-58606,-12.6506], Tmin=(867.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-477.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C(F)C[C]=CF(8522)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {7,D} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-73.6978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.14557,'amu*angstrom^2'), symmetry=1, barrier=(3.34695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145559,'amu*angstrom^2'), symmetry=1, barrier=(3.34669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5748,0.0547926,-5.07732e-05,2.56555e-08,-5.40731e-12,-8777.69,22.9586], Tmin=(100,'K'), Tmax=(1116.58,'K')), NASAPolynomial(coeffs=[9.74025,0.0255409,-1.14768e-05,2.19307e-09,-1.54119e-13,-10601.2,-17.3354], Tmin=(1116.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.6978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]F(5942)',
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
        HinderedRotor(inertia=(0.900551,'amu*angstrom^2'), symmetry=1, barrier=(20.7054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79302e-05,1.96926e-08,-6.19852e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.897,'K')), NASAPolynomial(coeffs=[4.81186,0.0161179,-7.47821e-06,1.46101e-09,-1.03887e-13,23340.4,3.53844], Tmin=(747.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = '[CH2]C(F)(F)C=C=CF(9221)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {8,D} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {8,D} {12,S}
8  C u0 p0 c0 {5,D} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-265.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,113,247,382,1207,3490,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(0.254544,'amu*angstrom^2'), symmetry=1, barrier=(5.85247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00756105,'amu*angstrom^2'), symmetry=1, barrier=(5.89992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791831,0.0769537,-0.000113169,9.50139e-08,-3.27667e-11,-31870.7,25.0336], Tmin=(100,'K'), Tmax=(730.546,'K')), NASAPolynomial(coeffs=[8.78364,0.0314483,-1.61463e-05,3.20093e-09,-2.26902e-13,-32991.8,-10.694], Tmin=(730.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = 'C#CCC([CH2])(F)F(6933)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {7,T}
7  C u0 p0 c0 {6,T} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-106.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,182.411],'cm^-1')),
        HinderedRotor(inertia=(0.946456,'amu*angstrom^2'), symmetry=1, barrier=(22.3548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321393,'amu*angstrom^2'), symmetry=1, barrier=(7.58999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.97835,'amu*angstrom^2'), symmetry=1, barrier=(70.3538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918948,0.0699448,-8.72006e-05,5.86025e-08,-1.57635e-11,-12685.2,22.5214], Tmin=(100,'K'), Tmax=(907.772,'K')), NASAPolynomial(coeffs=[11.5386,0.02315,-9.87611e-06,1.8148e-09,-1.24017e-13,-14613.2,-27.6844], Tmin=(907.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)CC#CF(9222)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {8,T}
8  C u0 p0 c0 {3,S} {7,T}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-206.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,239,401,1367,304.452,2177.71],'cm^-1')),
        HinderedRotor(inertia=(0.138853,'amu*angstrom^2'), symmetry=1, barrier=(9.11937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138788,'amu*angstrom^2'), symmetry=1, barrier=(9.11874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779415,'amu*angstrom^2'), symmetry=1, barrier=(51.0748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02705,0.0738773,-9.1357e-05,4.48079e-08,2.41092e-12,-24694.4,24.5196], Tmin=(100,'K'), Tmax=(574.083,'K')), NASAPolynomial(coeffs=[9.86039,0.0282066,-1.35097e-05,2.57997e-09,-1.7842e-13,-25970.2,-15.4725], Tmin=(574.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)=C[C]=CF(9102)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {7,S} {8,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,D} {11,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-0.603809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,615,860,1140,1343,3152,1685,370,180,1342.4],'cm^-1')),
        HinderedRotor(inertia=(3.62265,'amu*angstrom^2'), symmetry=1, barrier=(83.2919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.6181,'amu*angstrom^2'), symmetry=1, barrier=(83.1873,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32549,0.0549621,-5.29543e-05,2.70331e-08,-5.51666e-12,26.9404,23.1365], Tmin=(100,'K'), Tmax=(1187.76,'K')), NASAPolynomial(coeffs=[12.1451,0.0185247,-6.93795e-06,1.20483e-09,-8.02873e-14,-2543.27,-30.9237], Tmin=(1187.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.603809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948448,0.0699033,-7.22044e-05,3.94689e-08,-8.85008e-12,-38759.4,24.3231], Tmin=(100,'K'), Tmax=(1062.75,'K')), NASAPolynomial(coeffs=[11.9964,0.0283204,-1.35127e-05,2.65125e-09,-1.89091e-13,-41107.7,-29.6494], Tmin=(1062.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Allyl_S) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)CC=[C]F(9224)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-214.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,215,315,519,588,595,1205,1248,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,167,640,1190,256.071,256.188],'cm^-1')),
        HinderedRotor(inertia=(0.00257411,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179818,'amu*angstrom^2'), symmetry=1, barrier=(8.37784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286772,'amu*angstrom^2'), symmetry=1, barrier=(13.3382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658102,0.0787521,-0.000100466,6.97601e-08,-1.97209e-11,-25641,26.9697], Tmin=(100,'K'), Tmax=(857.317,'K')), NASAPolynomial(coeffs=[11.2034,0.0295504,-1.438e-05,2.81758e-09,-1.99808e-13,-27449.1,-22.2819], Tmin=(857.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'CC(F)(F)[CH][C]=CF(9225)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {8,S} {12,S}
7  C u0 p0 c0 {3,S} {8,D} {13,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-288.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22716,0.0650316,-6.31697e-05,3.32051e-08,-7.33439e-12,-34601.1,23.6526], Tmin=(100,'K'), Tmax=(1061.37,'K')), NASAPolynomial(coeffs=[10.265,0.0309703,-1.50319e-05,2.96856e-09,-2.12307e-13,-36519.5,-20.488], Tmin=(1061.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Allyl_S) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = 'CC(F)(F)C[C]=[C]F(9226)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {5,S} {8,D}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-179.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,167,640,1190,365.734,365.751],'cm^-1')),
        HinderedRotor(inertia=(0.0866282,'amu*angstrom^2'), symmetry=1, barrier=(8.2229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866266,'amu*angstrom^2'), symmetry=1, barrier=(8.2229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866188,'amu*angstrom^2'), symmetry=1, barrier=(8.2229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866425,0.074721,-9.4435e-05,6.75344e-08,-2.00194e-11,-21479.6,26.551], Tmin=(100,'K'), Tmax=(813.374,'K')), NASAPolynomial(coeffs=[9.53133,0.0321088,-1.58507e-05,3.12426e-09,-2.22196e-13,-22889.1,-13.4623], Tmin=(813.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cdj(Cs-CsHH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=[C]C[C](F)CF(9227)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {8,D} {13,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-181.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,615,860,1140,1343,3152,1685,370,310.436,310.44,310.491,2523.6],'cm^-1')),
        HinderedRotor(inertia=(0.154774,'amu*angstrom^2'), symmetry=1, barrier=(10.5847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154754,'amu*angstrom^2'), symmetry=1, barrier=(10.5849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467888,'amu*angstrom^2'), symmetry=1, barrier=(31.9978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02447,0.0713055,-9.04604e-05,6.84239e-08,-2.18601e-11,-21707.5,29.0505], Tmin=(100,'K'), Tmax=(752.395,'K')), NASAPolynomial(coeffs=[7.95411,0.0344656,-1.70161e-05,3.34875e-09,-2.37705e-13,-22750.3,-2.40966], Tmin=(752.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + radical(CsCsCsF1s) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[CH2][C](F)CC(F)=CF(9228)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {6,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-230.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504897,0.0854468,-0.000136074,1.21478e-07,-4.28998e-11,-27635.3,29.0476], Tmin=(100,'K'), Tmax=(817.673,'K')), NASAPolynomial(coeffs=[7.84277,0.0352007,-1.7575e-05,3.40005e-09,-2.35902e-13,-28355.6,-1.94302], Tmin=(817.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH]=[C]CC(F)(F)CF(9229)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
7  C u1 p0 c0 {5,S} {8,D}
8  C u1 p0 c0 {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-168.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.662879,'amu*angstrom^2'), symmetry=1, barrier=(15.2409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66257,'amu*angstrom^2'), symmetry=1, barrier=(15.2338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662701,'amu*angstrom^2'), symmetry=1, barrier=(15.2368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580275,0.0811381,-0.000109931,8.21796e-08,-2.50338e-11,-20154.5,26.5486], Tmin=(100,'K'), Tmax=(797.761,'K')), NASAPolynomial(coeffs=[10.6326,0.0307325,-1.51506e-05,2.96995e-09,-2.10025e-13,-21758.3,-19.6765], Tmin=(797.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)CC([CH2])(F)F(9230)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u1 p0 c0 {6,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-222.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541705,0.0814441,-0.000105912,7.41073e-08,-2.10019e-11,-26667,26.7641], Tmin=(100,'K'), Tmax=(856.959,'K')), NASAPolynomial(coeffs=[11.7605,0.0290778,-1.42497e-05,2.79817e-09,-1.98591e-13,-28589.8,-25.6281], Tmin=(856.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsHH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    E0 = (-118.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (149.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (51.3948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (302.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (356.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-110.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-40.3959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (152.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-24.3019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (5.57624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (46.5354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (128.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (104.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (197.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (74.9587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (39.0363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (54.4887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-24.6385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-47.5509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (54.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (102.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (104.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (32.8961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['CH2CF2(56)', 'C=C=CF(5941)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC=[C]CC[C](F)F(8557)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['[CH2]C(=CF)C([CH2])(F)F(8552)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(18)', 'FC=[C]C[C](F)F(4639)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CF(1334)', '[CH2]C([CH2])(F)F(4485)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['FC=C1CC(F)(F)C1(9203)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['CC(F)(F)C=C=CF(9220)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'C=C(F)C[C]=CF(8522)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(54.3623,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2CF2(56)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.000816002,'m^3/(mol*s)'), n=2.63352, Ea=(42.2953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.23368250386652978, var=0.35313233945400097, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing_Ext-1R!H-R_N-Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing_Ext-1R!H-R_N-Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C](F)F(4220)', 'C=C=CF(5941)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.00136,'m^3/(mol*s)'), n=2.41, Ea=(6.26933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_4R!H->C_Ext-2C-R_Sp-2C=1C_N-Sp-6R!H-2C_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_4R!H->C_Ext-2C-R_Sp-2C=1C_N-Sp-6R!H-2C_Ext-3C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[CH2]C(F)(F)C=C=CF(9221)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(1.73391,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C#CCC([CH2])(F)F(6933)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(63.3095,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[CH2]C(F)(F)CC#CF(9222)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](F)F(4220)', 'C=[C][CH]F(5942)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(11.5628,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', '[CH2]C(F)=C[C]=CF(9102)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(257.781,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['[CH2]C(F)(F)[CH]C=CF(9223)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(F)(F)CC=[C]F(9224)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['CC(F)(F)[CH][C]=CF(9225)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.86259e+06,'s^-1'), n=1.72764, Ea=(94.0044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC(F)(F)C[C]=[C]F(9226)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC=[C]C[C](F)CF(9227)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(136.48,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['[CH2][C](F)CC(F)=CF(9228)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(221.49,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC(F)(F)CF(9229)'],
    products = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(174.262,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(F)(F)C[C]=CF(8554)'],
    products = ['[CH]=C(F)CC([CH2])(F)F(9230)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(151.539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #3283',
    isomers = [
        '[CH2]C(F)(F)C[C]=CF(8554)',
    ],
    reactants = [
        ('CH2CF2(56)', 'C=C=CF(5941)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3283',
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

