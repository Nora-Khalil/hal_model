species(
    label = 'O=C[C](F)OC([CH]F)=C(F)F(5252)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {9,S} {10,D}
8  C u1 p0 c0 {1,S} {5,S} {11,S}
9  C u1 p0 c0 {2,S} {7,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {7,D}
11 C u0 p0 c0 {6,D} {8,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-746.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,280,501,1494,1531,234,589,736,816,1240,3237,182,240,577,636,1210,1413,2782.5,750,1395,475,1775,1000,319.445,319.446,319.446],'cm^-1')),
        HinderedRotor(inertia=(0.580721,'amu*angstrom^2'), symmetry=1, barrier=(42.0522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580721,'amu*angstrom^2'), symmetry=1, barrier=(42.0521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580721,'amu*angstrom^2'), symmetry=1, barrier=(42.0521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0222735,'amu*angstrom^2'), symmetry=1, barrier=(42.0522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0331784,0.0956802,-0.000129353,9.42172e-08,-2.80929e-11,-89616.6,31.74], Tmin=(100,'K'), Tmax=(811.473,'K')), NASAPolynomial(coeffs=[12.0902,0.0362472,-1.94908e-05,3.95943e-09,-2.85904e-13,-91573.4,-23.9092], Tmin=(811.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCOF1sO2s) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
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
    label = 'FC=C=C(F)F(1325)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C[C]F(4375)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u2 p0 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (22.7053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,163,1167],'cm^-1')),
        HinderedRotor(inertia=(0.0337629,'amu*angstrom^2'), symmetry=1, barrier=(13.6228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34132,0.0165481,-1.72264e-05,1.26789e-08,-4.5452e-12,2752.64,10.5202], Tmin=(100,'K'), Tmax=(638.16,'K')), NASAPolynomial(coeffs=[4.07864,0.0119262,-6.36171e-06,1.32811e-09,-9.81565e-14,2658.54,7.2943], Tmin=(638.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = 'O=C([CH]F)[C](F)F(3211)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 C u1 p0 c0 {1,S} {5,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-459.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,235,1215,1347,1486,3221,179,346,818,1406,1524,387.551,389.302],'cm^-1')),
        HinderedRotor(inertia=(0.00112358,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.483289,'amu*angstrom^2'), symmetry=1, barrier=(51.7608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0822,0.0453167,-5.36925e-05,3.38513e-08,-8.69315e-12,-55175.3,19.6668], Tmin=(100,'K'), Tmax=(937.879,'K')), NASAPolynomial(coeffs=[8.71354,0.0170346,-8.45973e-06,1.69905e-09,-1.22736e-13,-56419.1,-11.9003], Tmin=(937.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(CO-CsO2d)(F1s)(H)) + radical(Csj(CO-CsO2d)(F1s)(F1s))"""),
)

species(
    label = '[CH]F(804)',
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
    label = 'O=C[C](F)O[C]=C(F)F(5470)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {7,D}
6  C u1 p0 c0 {1,S} {4,S} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-450.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,562,600,623,1070,1265,1685,370,245.849,245.849,245.849,245.849],'cm^-1')),
        HinderedRotor(inertia=(0.943986,'amu*angstrom^2'), symmetry=1, barrier=(40.4882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.943985,'amu*angstrom^2'), symmetry=1, barrier=(40.4882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853104,'amu*angstrom^2'), symmetry=1, barrier=(36.5902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93215,0.071143,-8.94606e-05,5.68771e-08,-1.44211e-11,-54055.2,25.4077], Tmin=(100,'K'), Tmax=(958.769,'K')), NASAPolynomial(coeffs=[12.9863,0.0208534,-1.07835e-05,2.17094e-09,-1.56658e-13,-56366.7,-32.2395], Tmin=(958.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCOF1sO2s) + radical(Cdj(Cd-F1sF1s)(O2s-Cs))"""),
)

species(
    label = 'O=CC1(F)OC(=C(F)F)C1F(5266)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-984.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254543,0.0821829,-7.31426e-05,2.23892e-08,8.63964e-13,-118235,25.5511], Tmin=(100,'K'), Tmax=(978.586,'K')), NASAPolynomial(coeffs=[22.5004,0.0132901,-4.51156e-06,8.19596e-10,-5.98287e-14,-123844,-89.636], Tmin=(978.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-984.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCCFO) + group(CsCCFH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(2methyleneoxetane)"""),
)

species(
    label = 'O=C=C(F)OC(CF)=C(F)F(5275)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-800.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.086242,0.105154,-0.000191085,1.86542e-07,-6.98125e-11,-96164.5,34.2882], Tmin=(100,'K'), Tmax=(823.358,'K')), NASAPolynomial(coeffs=[5.46888,0.0461424,-2.52353e-05,5.01526e-09,-3.51455e-13,-95993.8,15.1595], Tmin=(823.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-800.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C[C](F)O[C]1C(F)C1(F)F(5471)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {5,S} {7,S} {8,S}
10 C u1 p0 c0 {4,S} {5,S} {11,S}
11 C u0 p0 c0 {6,D} {10,S} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-662.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.115212,0.103386,-0.000178335,1.67492e-07,-6.10703e-11,-79498.8,31.8972], Tmin=(100,'K'), Tmax=(823.271,'K')), NASAPolynomial(coeffs=[7.11075,0.043152,-2.28102e-05,4.48108e-09,-3.12429e-13,-79837.2,3.61265], Tmin=(823.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-662.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(C2CsJOCs) + radical(CsCOF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH]C(OC1(F)[CH]O1)=C(F)F(5472)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {5,S} {10,S} {11,D}
9  C u1 p0 c0 {6,S} {7,S} {12,S}
10 C u1 p0 c0 {2,S} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-600.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0104916,0.0918601,-0.000105833,6.05421e-08,-1.3832e-11,-72042.3,31.0217], Tmin=(100,'K'), Tmax=(1057.01,'K')), NASAPolynomial(coeffs=[16.9857,0.0275427,-1.45615e-05,2.97673e-09,-2.17028e-13,-75635.3,-51.9176], Tmin=(1057.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-600.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + ring(Cs(F)(O2)-O2s-Cs) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
)

species(
    label = 'O=CC1(F)O[C]([CH]F)C1(F)F(5473)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {9,S} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-737.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.361952,0.109389,-0.000194421,1.81643e-07,-6.54904e-11,-88522.2,33.3836], Tmin=(100,'K'), Tmax=(825.056,'K')), NASAPolynomial(coeffs=[8.62069,0.0401136,-2.17043e-05,4.29165e-09,-2.99729e-13,-89128.8,-2.91881], Tmin=(825.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCCFO) + group(CsCsFHH) + group(Cds-OdCsH) + ring(O2s-Cs-Cs-Cs(F)) + radical(C2CsJOCs) + radical(Csj(Cs-O2sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]OC(F)C(=C(F)F)O1(5399)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {6,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-770.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.81619,0.0906674,-8.48441e-05,3.28295e-08,-4.02678e-12,-92497.6,27.1697], Tmin=(100,'K'), Tmax=(1172.54,'K')), NASAPolynomial(coeffs=[27.4223,0.0121209,-7.11517e-06,1.57242e-09,-1.20734e-13,-100342,-118.773], Tmin=(1172.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-770.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCFF) + ring(Cyclohexanone) + radical(CCsJOCs) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=CC1(F)OC1([CH]F)[C](F)F(5474)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u1 p0 c0 {2,S} {7,S} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {7,S}
11 C u0 p0 c0 {6,D} {8,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-694.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.249336,0.104078,-0.000153404,1.09679e-07,-2.67813e-11,-83375.8,32.8567], Tmin=(100,'K'), Tmax=(617.358,'K')), NASAPolynomial(coeffs=[13.1786,0.0325823,-1.7368e-05,3.4639e-09,-2.45446e-13,-85329.3,-27.8423], Tmin=(617.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-694.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(CsCCFO) + group(CsCsFHH) + group(CsCsFFH) + group(Cds-OdCsH) + ring(O2s-Cs(F)-Cs(C)) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
)

species(
    label = '[O]C1[C](F)OC(=C(F)F)C1F(5451)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {11,D}
10 C u1 p0 c0 {2,S} {5,S} {8,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-715.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.160496,0.0785604,-6.48442e-05,1.70489e-08,1.36205e-12,-85950.3,29.3495], Tmin=(100,'K'), Tmax=(1038.84,'K')), NASAPolynomial(coeffs=[22.6438,0.014337,-6.16369e-06,1.24437e-09,-9.36257e-14,-91960.9,-87.6619], Tmin=(1038.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCFF) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[CH][C]=C(F)F(2842)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-200.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.793],'cm^-1')),
        HinderedRotor(inertia=(0.35565,'amu*angstrom^2'), symmetry=1, barrier=(8.17708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377452,-4.40203e-05,2.68135e-08,-6.58286e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.929,'K')), NASAPolynomial(coeffs=[8.46193,0.0130101,-6.31222e-06,1.26456e-09,-9.14162e-14,-25275.2,-12.1012], Tmin=(983.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = '[O][CH]C(=O)F(398)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.101,381.695],'cm^-1')),
        HinderedRotor(inertia=(0.482775,'amu*angstrom^2'), symmetry=1, barrier=(49.7784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80088,0.0290439,-3.02832e-05,1.66615e-08,-3.81631e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.22,'K')), NASAPolynomial(coeffs=[6.95396,0.0128879,-6.71487e-06,1.38086e-09,-1.01081e-13,-26608.7,-6.32755], Tmin=(1028.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
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
    label = 'O=C=C(F)OC([CH]F)=C(F)F(5475)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {1,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-650.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,197,221,431,657,2120,512.5,787.5,180.442,606.885,609.162,611.71],'cm^-1')),
        HinderedRotor(inertia=(0.00867212,'amu*angstrom^2'), symmetry=1, barrier=(2.3082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00816231,'amu*angstrom^2'), symmetry=1, barrier=(2.21442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00832295,'amu*angstrom^2'), symmetry=1, barrier=(2.21705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.113278,0.104626,-0.000190198,1.82531e-07,-6.73967e-11,-78094.4,34.9476], Tmin=(100,'K'), Tmax=(819.752,'K')), NASAPolynomial(coeffs=[7.0341,0.0415211,-2.30711e-05,4.60598e-09,-3.23462e-13,-78317.7,7.6715], Tmin=(819.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
)

species(
    label = 'O=C[C]F-2(1215)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=[C]C(F)OC([CH]F)=C(F)F(5476)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-726.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([355,410,600,1181,1341,1420,3056,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,1855,455,950,320.553,320.553,320.553],'cm^-1')),
        HinderedRotor(inertia=(0.104774,'amu*angstrom^2'), symmetry=1, barrier=(7.63975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553888,'amu*angstrom^2'), symmetry=1, barrier=(40.3877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230779,'amu*angstrom^2'), symmetry=1, barrier=(16.8276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553888,'amu*angstrom^2'), symmetry=1, barrier=(40.3877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.774879,0.11354,-0.000178837,1.4396e-07,-4.58782e-11,-87264.9,34.0837], Tmin=(100,'K'), Tmax=(770.799,'K')), NASAPolynomial(coeffs=[15.3328,0.0299423,-1.61383e-05,3.2275e-09,-2.28933e-13,-89747.8,-39.4318], Tmin=(770.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cd-O2sCd)(F1s)(H)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C][C](F)OC(CF)=C(F)F(5477)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 C u1 p0 c0 {4,S} {5,S} {11,S}
11 C u1 p0 c0 {6,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-737.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,280,501,1494,1531,1855,455,950,200.742,200.752,200.822,2597.96],'cm^-1')),
        HinderedRotor(inertia=(0.504435,'amu*angstrom^2'), symmetry=1, barrier=(14.4234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41768,'amu*angstrom^2'), symmetry=1, barrier=(40.5651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41662,'amu*angstrom^2'), symmetry=1, barrier=(40.5643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41681,'amu*angstrom^2'), symmetry=1, barrier=(40.5655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.56615,0.111431,-0.000189454,1.68535e-07,-5.85436e-11,-88598.1,32.9776], Tmin=(100,'K'), Tmax=(817.785,'K')), NASAPolynomial(coeffs=[11.3617,0.0357984,-1.90126e-05,3.73661e-09,-2.60477e-13,-89970.8,-18.6333], Tmin=(817.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)(F)OC([C]F)=CF(5478)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {5,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 C u2 p0 c0 {4,S} {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-695.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,350,440,435,1725,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.419429,0.106289,-0.000161025,1.2704e-07,-4.01988e-11,-83503.7,31.725], Tmin=(100,'K'), Tmax=(772.136,'K')), NASAPolynomial(coeffs=[13.7745,0.0327584,-1.81811e-05,3.70892e-09,-2.67136e-13,-85695.7,-33.0821], Tmin=(772.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'O=C[C](F)OC(=[C]F)C(F)F(5479)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {11,D}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {13,S}
11 C u1 p0 c0 {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-660.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,280,501,1494,1531,2782.5,750,1395,475,1775,1000,167,640,1190,195.243,195.244,195.244,1591.46],'cm^-1')),
        HinderedRotor(inertia=(0.00442238,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55563,'amu*angstrom^2'), symmetry=1, barrier=(42.083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55601,'amu*angstrom^2'), symmetry=1, barrier=(42.0823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55579,'amu*angstrom^2'), symmetry=1, barrier=(42.0838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.359385,0.107186,-0.00018077,1.64752e-07,-5.9256e-11,-79253.8,33.0972], Tmin=(100,'K'), Tmax=(795.582,'K')), NASAPolynomial(coeffs=[9.56277,0.0400868,-2.18059e-05,4.35127e-09,-3.06908e-13,-80287.9,-9.07886], Tmin=(795.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-660.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCOF1sO2s) + radical(Cdj(Cd-CsO2s)(F1s))"""),
)

species(
    label = 'F[C]F(138)',
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
    label = 'O=C[C](F)O[C]=CF(5037)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {8,S}
4  O u0 p2 c0 {6,D}
5  C u1 p0 c0 {1,S} {3,S} {6,S}
6  C u0 p0 c0 {4,D} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {8,D} {10,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-259.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,275.827,275.845,275.871,275.997],'cm^-1')),
        HinderedRotor(inertia=(1.01229,'amu*angstrom^2'), symmetry=1, barrier=(54.6592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.541843,'amu*angstrom^2'), symmetry=1, barrier=(29.2631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.541748,'amu*angstrom^2'), symmetry=1, barrier=(29.2676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20654,0.0636765,-7.49502e-05,4.49118e-08,-1.0758e-11,-31096.7,23.2643], Tmin=(100,'K'), Tmax=(1011.99,'K')), NASAPolynomial(coeffs=[12.3452,0.0196498,-9.69239e-06,1.92198e-09,-1.37826e-13,-33351.1,-30.6058], Tmin=(1011.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCOF1sO2s) + radical(Cdj(Cd-F1sH)(O2s-Cs))"""),
)

species(
    label = 'O=CC1(F)OC(=CF)C1(F)F(5268)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1012.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.923672,0.0914441,-9.92769e-05,5.0804e-08,-9.80256e-12,-121558,28.1339], Tmin=(100,'K'), Tmax=(1386.94,'K')), NASAPolynomial(coeffs=[25.3891,0.00840626,-1.73673e-06,2.01762e-10,-1.13254e-14,-128169,-104.938], Tmin=(1386.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1012.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCCFO) + group(CsCCFF) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(2methyleneoxetane)"""),
)

species(
    label = 'O=C=C(F)OC(=CF)C(F)F(5276)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {6,D} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-825.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147828,0.109013,-0.000206544,2.0595e-07,-7.76263e-11,-99159.9,34.0521], Tmin=(100,'K'), Tmax=(832.961,'K')), NASAPolynomial(coeffs=[4.23171,0.0485916,-2.68041e-05,5.32234e-09,-3.71858e-13,-98523,21.9264], Tmin=(832.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-825.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=CC1(F)O[C]([C](F)F)C1F(5480)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
9  C u1 p0 c0 {5,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {13,S}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-718.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0653016,0.0992262,-0.000162384,1.45787e-07,-5.20715e-11,-86218.9,33.6713], Tmin=(100,'K'), Tmax=(787.96,'K')), NASAPolynomial(coeffs=[9.40391,0.0378375,-2.01661e-05,4.01002e-09,-2.82843e-13,-87297.6,-7.1319], Tmin=(787.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-718.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCCFO) + group(CsCsFFH) + group(Cds-OdCsH) + ring(O2s-Cs-Cs-Cs(F)) + radical(C2CsJOCs) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C1O[C](F)[CH]OC1(F)F(5303)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u1 p0 c0 {5,S} {10,S} {12,S}
10 C u1 p0 c0 {3,S} {6,S} {9,S}
11 C u0 p0 c0 {4,S} {8,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-804.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03347,0.0955712,-9.80576e-05,4.52327e-08,-7.94602e-12,-96540.2,27.8493], Tmin=(100,'K'), Tmax=(1400.82,'K')), NASAPolynomial(coeffs=[29.2014,0.00923546,-5.60809e-06,1.23439e-09,-9.3685e-14,-105011,-128.207], Tmin=(1400.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-804.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFFO) + group(Cds-CdsCsOs) + group(CdCFH) + ring(Cyclohexanone) + radical(CCsJOCs) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1[C](F)OC(=CF)C1(F)F(5351)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {11,D}
10 C u1 p0 c0 {3,S} {5,S} {8,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-742.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.31389,0.0820167,-7.19411e-05,2.22174e-08,7.52783e-14,-89122.3,30.3664], Tmin=(100,'K'), Tmax=(1033.18,'K')), NASAPolynomial(coeffs=[23.5064,0.0130269,-5.50761e-06,1.11384e-09,-8.43171e-14,-95284.4,-91.3307], Tmin=(1033.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFF) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCFH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
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
    label = 'O=CC(F)OC([C]F)=C(F)F(5481)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {8,D}
11 C u2 p0 c0 {4,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-665.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,350,440,435,1725,2782.5,750,1395,475,1775,1000,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288652,0.103537,-0.00015361,1.20373e-07,-3.81148e-11,-79849.3,31.7266], Tmin=(100,'K'), Tmax=(769.534,'K')), NASAPolynomial(coeffs=[12.9767,0.0345828,-1.91985e-05,3.92565e-09,-2.83462e-13,-81890.9,-28.7958], Tmin=(769.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-665.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'O=[C][C](F)OC(=CF)C(F)F(5482)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {13,S}
10 C u1 p0 c0 {4,S} {5,S} {11,S}
11 C u1 p0 c0 {6,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-762.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.634407,0.115376,-0.000205257,1.88457e-07,-6.66093e-11,-91593.2,32.765], Tmin=(100,'K'), Tmax=(833.164,'K')), NASAPolynomial(coeffs=[10.1396,0.0382204,-2.05651e-05,4.03972e-09,-2.80543e-13,-92505.9,-11.9502], Tmin=(833.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-762.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[CH]C(OC(F)(F)C=O)=C(F)F(5483)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {5,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u2 p0 c0 {8,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-732.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,350,440,435,1725,2782.5,750,1395,475,1775,1000,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299018,0.101305,-0.000128416,8.57003e-08,-2.31741e-11,-87899,32.7304], Tmin=(100,'K'), Tmax=(895.523,'K')), NASAPolynomial(coeffs=[14.437,0.0354822,-1.81604e-05,3.61887e-09,-2.59074e-13,-90538.2,-36.7356], Tmin=(895.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-732.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O[C](F)C=O)C(F)(F)F(5484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {11,D}
9  C u1 p0 c0 {4,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 C u1 p0 c0 {8,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-756.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.256016,0.0966622,-0.000122251,7.67383e-08,-1.89775e-11,-90839.4,31.0034], Tmin=(100,'K'), Tmax=(988.226,'K')), NASAPolynomial(coeffs=[17.6784,0.0240698,-1.2065e-05,2.40601e-09,-1.73017e-13,-94384,-55.3075], Tmin=(988.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CdOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCOF1sO2s) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C(F)C(F)(F)C(=O)[CH]F(5485)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u1 p0 c0 {4,S} {8,S} {12,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-798.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.291683,0.101622,-0.000146065,1.08374e-07,-3.20509e-11,-95891.3,32.1938], Tmin=(100,'K'), Tmax=(826.881,'K')), NASAPolynomial(coeffs=[14.5878,0.0296421,-1.54899e-05,3.0971e-09,-2.21132e-13,-98352,-36.7623], Tmin=(826.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-798.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Csj(CO-CsO2d)(F1s)(H))"""),
)

species(
    label = '[O]C(C(=O)F)C([CH]F)=C(F)F(5249)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 C u1 p0 c0 {2,S} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-742.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269151,0.0881154,-0.000112605,7.54746e-08,-2.04419e-11,-89137.5,36.7913], Tmin=(100,'K'), Tmax=(894.936,'K')), NASAPolynomial(coeffs=[13.2379,0.0301492,-1.54461e-05,3.09615e-09,-2.22575e-13,-91458.7,-24.3354], Tmin=(894.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = '[CH]=C(F)OC([CH]F)=C(F)F(5486)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {1,S} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-357.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,293,496,537,1218,3120,650,792.5,1650,180,1727.73,1728.67],'cm^-1')),
        HinderedRotor(inertia=(2.14841,'amu*angstrom^2'), symmetry=1, barrier=(49.3961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14725,'amu*angstrom^2'), symmetry=1, barrier=(49.3695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113967,'amu*angstrom^2'), symmetry=1, barrier=(2.62033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0171202,0.100596,-0.000179592,1.70888e-07,-6.29182e-11,-42886,33.0934], Tmin=(100,'K'), Tmax=(815.201,'K')), NASAPolynomial(coeffs=[7.25725,0.0400336,-2.20876e-05,4.4079e-09,-3.09857e-13,-43234.5,4.74566], Tmin=(815.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cds-CdsHH) + radical(Csj(Cd-O2sCd)(F1s)(H)) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'FC1=COC(F)C(=C(F)F)O1(5413)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {9,D}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-967.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.257638,0.0681376,-3.66831e-05,-2.2892e-09,4.86577e-12,-116261,23.1484], Tmin=(100,'K'), Tmax=(1164.52,'K')), NASAPolynomial(coeffs=[20.3944,0.0242536,-1.27244e-05,2.63913e-09,-1.94806e-13,-122665,-84.4278], Tmin=(1164.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-967.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFO) + group(CdCFF) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclohexane)"""),
)

species(
    label = 'F[CH]C1([C](F)F)OC=C(F)O1(5487)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
8  C u1 p0 c0 {2,S} {7,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u0 p0 c0 {1,S} {5,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-726.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09863,0.0838111,-2.6429e-05,-6.18264e-08,3.92681e-11,-87203.3,30.5245], Tmin=(100,'K'), Tmax=(927.549,'K')), NASAPolynomial(coeffs=[37.7394,-0.0105386,7.87537e-06,-1.47357e-09,8.91671e-14,-97554.3,-170.885], Tmin=(927.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsOs) + group(CsCsFHH) + group(CsCsFFH) + group(CdCFO) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CsCsF1sH) + radical(CsCsF1sF1s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'O=C=[C]F(580)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (33.6712,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(58.9933,'amu')),
        NonlinearRotor(inertia=([3.58812,117.154,120.742],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([280.764,365.894,563.98,866.433,1429.66,2068.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94745,0.00353157,4.23087e-05,-1.10022e-07,8.19717e-11,4052.13,8.19715], Tmin=(10,'K'), Tmax=(472.411,'K')), NASAPolynomial(coeffs=[4.41104,0.00925722,-6.51505e-06,2.12227e-09,-2.60232e-13,3900.64,5.16846], Tmin=(472.411,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(33.6712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C([CH]F)C(F)F(3959)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-623.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,375,552.5,462.5,1710,235,1215,1347,1486,3221,520.144,520.166],'cm^-1')),
        HinderedRotor(inertia=(0.0765802,'amu*angstrom^2'), symmetry=1, barrier=(14.7033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249467,'amu*angstrom^2'), symmetry=1, barrier=(47.8942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70021,0.0289351,5.81793e-06,-4.08324e-08,2.36932e-11,-75027.4,13.7673], Tmin=(10,'K'), Tmax=(730.67,'K')), NASAPolynomial(coeffs=[6.5202,0.0268284,-1.72251e-05,5.16265e-09,-5.87804e-13,-75795.4,-1.3878], Tmin=(730.67,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-623.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""ODC([CH]F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C=[C]OC([CH]F)=C(F)F(5488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {4,S} {7,S} {8,D}
7  C u1 p0 c0 {1,S} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u1 p0 c0 {4,S} {10,D}
10 C u0 p0 c0 {5,D} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-222.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,1685,370,2120,512.5,787.5,180,180,1936.02],'cm^-1')),
        HinderedRotor(inertia=(0.00243313,'amu*angstrom^2'), symmetry=1, barrier=(6.46847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06141,'amu*angstrom^2'), symmetry=1, barrier=(47.3959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0624,'amu*angstrom^2'), symmetry=1, barrier=(47.4186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234215,0.0960926,-0.000177667,1.69133e-07,-6.1296e-11,-26635.4,34.5867], Tmin=(100,'K'), Tmax=(838.472,'K')), NASAPolynomial(coeffs=[7.20653,0.0351684,-1.91889e-05,3.781e-09,-2.62623e-13,-26832.2,7.97619], Tmin=(838.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cd-O2sCd)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=C(F)OC([CH]F)=C(F)F(5489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {11,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {1,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-562.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,293,496,537,1218,1685,370,180,180,1585.87,1585.99],'cm^-1')),
        HinderedRotor(inertia=(0.123664,'amu*angstrom^2'), symmetry=1, barrier=(2.84329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123474,'amu*angstrom^2'), symmetry=1, barrier=(2.83891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11397,'amu*angstrom^2'), symmetry=1, barrier=(48.6042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11385,'amu*angstrom^2'), symmetry=1, barrier=(48.6016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.537245,0.113927,-0.000204525,1.9165e-07,-6.92535e-11,-67487.5,38.5665], Tmin=(100,'K'), Tmax=(822.831,'K')), NASAPolynomial(coeffs=[8.92833,0.0410102,-2.25563e-05,4.48168e-09,-3.13718e-13,-68134.5,0.280443], Tmin=(822.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-562.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cd-O2sCd)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = 'F[CH]C(O[C]=COF)=C(F)F(5490)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {4,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {1,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-246.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,3010,987.5,1337.5,450,1655,1685,370,272.131,281.628,1918.79,1919.05],'cm^-1')),
        HinderedRotor(inertia=(0.190788,'amu*angstrom^2'), symmetry=1, barrier=(10.2752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706939,'amu*angstrom^2'), symmetry=1, barrier=(37.4617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21806,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.74923,'amu*angstrom^2'), symmetry=1, barrier=(37.4255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.575092,0.110889,-0.000185615,1.63216e-07,-5.64334e-11,-29450.1,39.1558], Tmin=(100,'K'), Tmax=(804.229,'K')), NASAPolynomial(coeffs=[11.8983,0.0351403,-1.87622e-05,3.7063e-09,-2.59586e-13,-31013,-15.5471], Tmin=(804.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2sCF) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cd-O2sCd)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC(=C(F)F)C(F)F(5491)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-630.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,182,240,577,636,1210,1413,3010,987.5,1337.5,450,1655,1685,370,338.993,338.993,338.994,338.994,2743.68],'cm^-1')),
        HinderedRotor(inertia=(0.461437,'amu*angstrom^2'), symmetry=1, barrier=(37.6291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0678153,'amu*angstrom^2'), symmetry=1, barrier=(5.53016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461438,'amu*angstrom^2'), symmetry=1, barrier=(37.6291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369429,0.106264,-0.000176703,1.56828e-07,-5.49424e-11,-75626.6,37.4196], Tmin=(100,'K'), Tmax=(798.583,'K')), NASAPolynomial(coeffs=[10.8387,0.0363149,-1.93785e-05,3.83837e-09,-2.69586e-13,-76976.4,-11.3752], Tmin=(798.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-630.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'F[C]C(=CF)OC(F)=COF(5492)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {5,S} {10,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u0 p0 c0 {2,S} {7,D} {13,S}
11 C u2 p0 c0 {4,S} {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-231.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,326,540,652,719,1357,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.466494,0.116258,-0.000221254,2.18033e-07,-8.14556e-11,-27706.7,35.625], Tmin=(100,'K'), Tmax=(829.472,'K')), NASAPolynomial(coeffs=[5.77463,0.0474122,-2.66816e-05,5.33122e-09,-3.73606e-13,-27409.1,14.7174], Tmin=(829.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2sCF) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C=C(F)C(F)C(=O)[C](F)F(5493)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {7,S} {11,D}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-788.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.408997,0.104588,-0.00015627,1.20995e-07,-3.7295e-11,-94651.1,33.6968], Tmin=(100,'K'), Tmax=(794.806,'K')), NASAPolynomial(coeffs=[14.3844,0.0301361,-1.57589e-05,3.13572e-09,-2.22636e-13,-97002.7,-34.2751], Tmin=(794.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-788.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Csj(CO-CsO2d)(F1s)(F1s))"""),
)

species(
    label = 'FC=C1OC(F)=COC1(F)F(5317)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {11,D}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 C u0 p0 c0 {4,S} {8,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1001.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0848401,0.0725683,-4.84837e-05,8.61846e-09,1.44855e-12,-120305,24.3582], Tmin=(100,'K'), Tmax=(1207.3,'K')), NASAPolynomial(coeffs=[21.2988,0.0227211,-1.19452e-05,2.46437e-09,-1.80776e-13,-126917,-88.1524], Tmin=(1207.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1001.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFO) + group(CdCFH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclohexane)"""),
)

species(
    label = '[O]CC(=O)F(396)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-363.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,180,475.353],'cm^-1')),
        HinderedRotor(inertia=(0.0192623,'amu*angstrom^2'), symmetry=1, barrier=(3.18239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86108,0.0121233,3.70793e-05,-8.49933e-08,5.10894e-11,-43750.9,11.7734], Tmin=(10,'K'), Tmax=(581.186,'K')), NASAPolynomial(coeffs=[3.98792,0.0215968,-1.40743e-05,4.31471e-09,-5.0294e-13,-43940.4,9.72707], Tmin=(581.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-363.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]=C=C(F)F(2804)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {4,D} {6,D}
6 C u1 p0 c0 {3,S} {5,D}
"""),
    E0 = (-151.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([94,120,354,641,825,1294,540,610,2055,137,207,812],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81697,0.0142144,8.17672e-05,-2.81639e-07,2.55173e-10,-18237.7,11.4722], Tmin=(10,'K'), Tmax=(406.163,'K')), NASAPolynomial(coeffs=[5.51158,0.0188349,-1.39947e-05,4.71475e-09,-5.89764e-13,-18551.1,2.65979], Tmin=(406.163,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-151.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""F[C]DCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(CF)[C](F)F(3844)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-607.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,375,552.5,462.5,1710,179,346,818,1406,1524,410.471,410.497],'cm^-1')),
        HinderedRotor(inertia=(0.05975,'amu*angstrom^2'), symmetry=1, barrier=(7.14244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301348,'amu*angstrom^2'), symmetry=1, barrier=(36.0284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.58756,0.0356021,-1.68932e-05,-1.76866e-08,1.74073e-11,-73004.1,13.2839], Tmin=(10,'K'), Tmax=(630.174,'K')), NASAPolynomial(coeffs=[6.36177,0.026905,-1.74046e-05,5.29611e-09,-6.13262e-13,-73530.7,-0.223027], Tmin=(630.174,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-607.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""ODC(CF)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC=C(F)OC([C]F)=C(F)F(5494)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {5,S} {10,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {7,D}
11 C u2 p0 c0 {4,S} {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-581.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,326,540,652,719,1357,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.639721,0.114193,-0.000199037,1.8101e-07,-6.42616e-11,-69814.7,34.8608], Tmin=(100,'K'), Tmax=(806.393,'K')), NASAPolynomial(coeffs=[10.952,0.0375243,-2.0767e-05,4.1512e-09,-2.9223e-13,-71061,-14.7039], Tmin=(806.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C=[C]OC(=CF)C(F)(F)F(5495)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {12,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-684.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,350,440,435,1725,194,682,905,1196,1383,3221,3010,987.5,1337.5,450,1655,1685,370,389.844,392.217,392.452,394.146,2352.99],'cm^-1')),
        HinderedRotor(inertia=(0.113807,'amu*angstrom^2'), symmetry=1, barrier=(12.2805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116177,'amu*angstrom^2'), symmetry=1, barrier=(12.1099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373462,'amu*angstrom^2'), symmetry=1, barrier=(39.6737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.118134,0.096291,-0.000130537,9.05095e-08,-2.49528e-11,-82228.8,35.7753], Tmin=(100,'K'), Tmax=(886.124,'K')), NASAPolynomial(coeffs=[15.1431,0.0273993,-1.39161e-05,2.76823e-09,-1.97821e-13,-84933.4,-36.0056], Tmin=(886.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(OC(F)=COF)=C(F)F(5496)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {5,S} {10,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {7,D}
11 C u2 p0 c0 {7,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-268.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,326,540,652,719,1357,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.463147,0.112699,-0.000193895,1.84064e-07,-6.79307e-11,-32097,37.0476], Tmin=(100,'K'), Tmax=(820.806,'K')), NASAPolynomial(coeffs=[6.41629,0.0501926,-2.67034e-05,5.25304e-09,-3.66649e-13,-32250,11.1638], Tmin=(820.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(AllylJ2_triplet)"""),
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
    E0 = (-257.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (52.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (253.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-249.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-232.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-26.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-74.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-130.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-187.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-202.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-195.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-159.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-90.4296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (50.1232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (73.6064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (65.0976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (177.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-101.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-198.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-33.4756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (15.2299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (263.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-249.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-232.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-130.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-187.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-195.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (25.6626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-132.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-193.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-52.7556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-34.9463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-121.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-95.1975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (374.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-249.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-178.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-1.2721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (240.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (89.1607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (295.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (3.63028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (290.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-116.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-249.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (48.5755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (10.0851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-67.8776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-27.0997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (259.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[C]F(4375)', 'O=C([CH]F)[C](F)F(3211)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]F(804)', 'O=C[C](F)O[C]=C(F)F(5470)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC1(F)OC(=C(F)F)C1F(5266)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=C=C(F)OC(CF)=C(F)F(5275)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=C[C](F)O[C]1C(F)C1(F)F(5471)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['F[CH]C(OC1(F)[CH]O1)=C(F)F(5472)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC1(F)O[C]([CH]F)C1(F)F(5473)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['F[C]1[CH]OC(F)C(=C(F)F)O1(5399)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC1(F)OC1([CH]F)[C](F)F(5474)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.86885e+09,'s^-1'), n=0.611527, Ea=(54.7386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]C1[C](F)OC(=C(F)F)C1F(5451)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(61.9573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=CC(=O)F(4234)', 'F[CH][C]=C(F)F(2842)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(35.9403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][CH]C(=O)F(398)', 'FC=C=C(F)F(1325)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=C=C(F)OC([CH]F)=C(F)F(5475)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(398)', 'F[CH][C]=C(F)F(2842)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[C]F-2(1215)', 'O=C([CH]F)[C](F)F(3211)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHF(40)', 'O=C[C](F)O[C]=C(F)F(5470)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]C(F)OC([CH]F)=C(F)F(5476)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_noH] for rate rule [R2H_S;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C][C](F)OC(CF)=C(F)F(5477)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6306.05,'s^-1'), n=2.09377, Ea=(50.431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_1H] for rate rule [R5HJ_1;CO_rad_out;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=CC(F)(F)OC([C]F)=CF(5478)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(173.319,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C[C](F)OC(=[C]F)C(F)F(5479)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(186.649,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]F(138)', 'O=C[C](F)O[C]=CF(5037)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC1(F)OC(=CF)C1(F)F(5268)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=C=C(F)OC(=CF)C(F)F(5276)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=CC1(F)O[C]([C](F)F)C1F(5480)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['FC=C1O[C](F)[CH]OC1(F)F(5303)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]C1[C](F)OC(=CF)C1(F)F(5351)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(61.9573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CF2(43)', 'O=C[C](F)O[C]=CF(5037)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=CC(F)OC([C]F)=C(F)F(5481)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=[C][C](F)OC(=CF)C(F)F(5482)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(57774.9,'s^-1'), n=1.76812, Ea=(63.6725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;XH_out] for rate rule [R5HJ_3;C_rad_out_noH;CO_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(OC(F)(F)C=O)=C(F)F(5483)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(190.568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[CH]=C(O[C](F)C=O)C(F)(F)F(5484)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(222.545,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]C=C(F)C(F)(F)C(=O)[CH]F(5485)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(136.243,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]C(C(=O)F)C([CH]F)=C(F)F(5249)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(162.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(6)', '[CH]=C(F)OC([CH]F)=C(F)F(5486)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['FC1=COC(F)C(=C(F)F)O1(5413)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_1H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['F[CH]C1([C](F)F)OC=C(F)O1(5487)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=C=[C]F(580)', 'O=C([CH]F)C(F)F(3959)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.33333e+07,'s^-1'), n=1.2, Ea=(256.219,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O"""),
)

reaction(
    label = 'reaction39',
    reactants = ['HF(38)', 'O=C=[C]OC([CH]F)=C(F)F(5488)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(255.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O[C]=C(F)OC([CH]F)=C(F)F(5489)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['F[CH]C(O[C]=COF)=C(F)F(5490)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(52.9997,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]C=[C]OC(=C(F)F)C(F)F(5491)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(144.903,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction43',
    reactants = ['F[C]C(=CF)OC(F)=COF(5492)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(33.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]C=C(F)C(F)C(=O)[C](F)F(5493)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(141.077,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['FC=C1OC(F)=COC1(F)F(5317)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['[O]CC(=O)F(396)', 'F[C]=C=C(F)F(2804)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.33333e+07,'s^-1'), n=1.2, Ea=(306.067,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O"""),
)

reaction(
    label = 'reaction47',
    reactants = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    products = ['O=C=[C]F(580)', 'O=C(CF)[C](F)F(3844)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.33333e+07,'s^-1'), n=1.2, Ea=(267.576,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O"""),
)

reaction(
    label = 'reaction48',
    reactants = ['OC=C(F)OC([C]F)=C(F)F(5494)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(7.546e+09,'s^-1'), n=0.732, Ea=(25.1375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Cd_rad_out_single;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_single;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O]C=[C]OC(=CF)C(F)(F)F(5495)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(169.032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]C(OC(F)=COF)=C(F)F(5496)'],
    products = ['O=C[C](F)OC([CH]F)=C(F)F(5252)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(39.2429,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1531',
    isomers = [
        'O=C[C](F)OC([CH]F)=C(F)F(5252)',
    ],
    reactants = [
        ('O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1531',
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

