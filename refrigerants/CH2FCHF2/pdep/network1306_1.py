species(
    label = 'OC(F)([CH]F)OC(F)=CF(4048)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {11,S}
10 C u0 p0 c0 {4,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-944.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,326,540,652,719,1357,334,575,1197,1424,3202,194,682,905,1196,1383,3221,233.581,233.972,234.04,234.154],'cm^-1')),
        HinderedRotor(inertia=(0.538023,'amu*angstrom^2'), symmetry=1, barrier=(20.8943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265326,'amu*angstrom^2'), symmetry=1, barrier=(10.3069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265174,'amu*angstrom^2'), symmetry=1, barrier=(10.3077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868514,'amu*angstrom^2'), symmetry=1, barrier=(33.6838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.782827,0.112505,-0.00017361,1.33345e-07,-4.01461e-11,-113480,30.8243], Tmin=(100,'K'), Tmax=(817.789,'K')), NASAPolynomial(coeffs=[17.0818,0.025116,-1.33034e-05,2.64851e-09,-1.87669e-13,-116401,-51.7668], Tmin=(817.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-944.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'O=C(F)[CH]F(215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-442.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(1.5365,'amu*angstrom^2'), symmetry=1, barrier=(35.3272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91969,0.0052474,7.61145e-05,-1.81228e-07,1.26791e-10,-53179.7,10.2212], Tmin=(10,'K'), Tmax=(489.215,'K')), NASAPolynomial(coeffs=[4.0714,0.0192239,-1.3397e-05,4.3331e-09,-5.27008e-13,-53376.6,7.73666], Tmin=(489.215,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-442.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)CF(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(F)C(O)(F)[CH]F(4958)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u1 p0 c0 {3,S} {7,S} {12,S}
10 C u0 p0 c0 {4,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1043.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.636765,0.113863,-0.000199872,1.79823e-07,-6.26834e-11,-125293,31.8183], Tmin=(100,'K'), Tmax=(822.438,'K')), NASAPolynomial(coeffs=[11.6638,0.0342663,-1.86403e-05,3.68687e-09,-2.57361e-13,-126648,-21.055], Tmin=(822.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1043.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + radical(Csj(Cs-F1sO2sCs)(F1s)(H))"""),
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
    label = '[CH]C(O)(F)OC(F)=CF(4959)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {10,S}
9  C u2 p0 c0 {6,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-523.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,326,540,652,719,1357,194,682,905,1196,1383,3221,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194925,0.0969953,-0.000137361,9.54143e-08,-2.59111e-11,-62819.2,27.855], Tmin=(100,'K'), Tmax=(905.197,'K')), NASAPolynomial(coeffs=[17.0999,0.0205731,-1.07262e-05,2.15184e-09,-1.54392e-13,-65950.4,-53.8605], Tmin=(905.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(137)',
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
    label = 'O[C](F)OC(F)=CF(4047)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {6,D} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-721.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,482,586,761,1411,194,682,905,1196,1383,3221,180,307.254,1254.52],'cm^-1')),
        HinderedRotor(inertia=(0.612981,'amu*angstrom^2'), symmetry=1, barrier=(14.0936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612358,'amu*angstrom^2'), symmetry=1, barrier=(14.0793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61896,'amu*angstrom^2'), symmetry=1, barrier=(37.2231,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.942889,0.0681912,-8.43181e-05,5.10338e-08,-1.20967e-11,-86684,24.6253], Tmin=(100,'K'), Tmax=(1032.61,'K')), NASAPolynomial(coeffs=[14.5709,0.0154002,-7.6316e-06,1.52344e-09,-1.09889e-13,-89498.4,-41.5589], Tmin=(1032.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-721.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFHOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sO2sO2s)"""),
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
    label = 'OC(F)([C]F)OC(F)=CF-2(4960)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {11,S}
10 C u2 p0 c0 {4,S} {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-775.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,194,682,905,1196,1383,3221,90,150,250,247,323,377,431,433,1065,1274,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.502727,0.105838,-0.000159496,1.17787e-07,-3.40417e-11,-93095.2,28.8825], Tmin=(100,'K'), Tmax=(851.173,'K')), NASAPolynomial(coeffs=[17.2682,0.0223223,-1.23118e-05,2.50317e-09,-1.80046e-13,-96120.3,-53.9878], Tmin=(851.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-775.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCCl_triplet)"""),
)

species(
    label = 'OC1(F)O[C](F)C(F)C1F(4961)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {7,S} {10,S} {12,S}
10 C u1 p0 c0 {4,S} {5,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1013.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.225421,0.0794386,-8.40432e-05,4.18315e-08,-6.50095e-12,-121733,25.2879], Tmin=(100,'K'), Tmax=(877.898,'K')), NASAPolynomial(coeffs=[16.4092,0.0197756,-6.15198e-06,9.45453e-10,-5.87948e-14,-125117,-53.7715], Tmin=(877.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1013.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFOO) + group(CsCFHO) + ring(Tetrahydrofuran) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'OC1(F)OC(F)([CH]F)C1F(4962)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {9,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
10 C u1 p0 c0 {4,S} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-963.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206704,0.0944746,-0.000129818,7.92182e-08,-9.69033e-12,-115796,27.5995], Tmin=(100,'K'), Tmax=(577.987,'K')), NASAPolynomial(coeffs=[11.6744,0.0331284,-1.73688e-05,3.44698e-09,-2.43641e-13,-117422,-24.0422], Tmin=(577.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-963.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCFOO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(O2s-Cs-Cs-Cs(F)) + radical(Csj(Cs-F1sO2sCs)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'OC(=CF)OC(F)=CF(4055)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u0 p0 c0 {2,S} {6,D} {10,S}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-704.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,326,540,652,719,1357,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.173795,'amu*angstrom^2'), symmetry=1, barrier=(3.99589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17242,'amu*angstrom^2'), symmetry=1, barrier=(3.96428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174149,'amu*angstrom^2'), symmetry=1, barrier=(4.00403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.196339,0.0945141,-0.000164675,1.50905e-07,-5.34578e-11,-84592.8,28.1179], Tmin=(100,'K'), Tmax=(831.331,'K')), NASAPolynomial(coeffs=[8.66448,0.0333786,-1.75752e-05,3.43765e-09,-2.38588e-13,-85296.2,-6.93386], Tmin=(831.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OC(F)=CF(878)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {2,S} {4,D} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-508.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.537415,'amu*angstrom^2'), symmetry=1, barrier=(12.3562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3292.1,'J/mol'), sigma=(5.16209,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.22 K, Pc=54.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86684,0.00827899,9.49292e-05,-2.30798e-07,1.56874e-10,-61152.3,9.15162], Tmin=(10,'K'), Tmax=(532.988,'K')), NASAPolynomial(coeffs=[6.38193,0.0179589,-1.26771e-05,4.31723e-09,-5.56992e-13,-61826,-5.20451], Tmin=(532.988,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-508.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC=C(F)OC(F)=CF(4040)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  C u0 p0 c0 {4,S} {7,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-734.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,390,483,597,572,732,631,807,1275,1439,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,180,180,1710.16],'cm^-1')),
        HinderedRotor(inertia=(0.174376,'amu*angstrom^2'), symmetry=1, barrier=(4.00925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177576,'amu*angstrom^2'), symmetry=1, barrier=(4.08282,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690355,0.0840957,-0.000149076,1.42733e-07,-5.27155e-11,-88212.3,26.5822], Tmin=(100,'K'), Tmax=(823.039,'K')), NASAPolynomial(coeffs=[5.94843,0.035688,-1.92017e-05,3.79706e-09,-2.65528e-13,-88303.8,6.94134], Tmin=(823.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-734.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O[C]([CH]F)OC(F)=CF(4963)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u1 p0 c0 {4,S} {5,S} {8,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {7,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-523.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,326,540,652,719,1357,334,575,1197,1424,3202,194,682,905,1196,1383,3221,280.279,280.28,280.28,1385.88],'cm^-1')),
        HinderedRotor(inertia=(0.171294,'amu*angstrom^2'), symmetry=1, barrier=(9.54879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171293,'amu*angstrom^2'), symmetry=1, barrier=(9.5488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171293,'amu*angstrom^2'), symmetry=1, barrier=(9.54879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39561,'amu*angstrom^2'), symmetry=1, barrier=(22.0534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.444272,0.106567,-0.000182649,1.5652e-07,-5.20827e-11,-62756.9,31.4203], Tmin=(100,'K'), Tmax=(819.453,'K')), NASAPolynomial(coeffs=[14.0183,0.0250625,-1.34874e-05,2.65293e-09,-1.84518e-13,-64761,-33.2384], Tmin=(819.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_P) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
)

species(
    label = 'OC(F)([CH]F)O[C]=CF(4964)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {6,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-510.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,615,860,1140,1343,3152,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.462247,0.103592,-0.000158297,1.17999e-07,-3.41764e-11,-61264.7,29.4858], Tmin=(100,'K'), Tmax=(852.267,'K')), NASAPolynomial(coeffs=[17.497,0.0193034,-9.94978e-06,1.95883e-09,-1.38e-13,-64325.9,-54.2864], Tmin=(852.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + radical(Cdj(Cd-F1sH)(O2s-Cs))"""),
)

species(
    label = '[CH]=C(F)OC(O)(F)[CH]F(4965)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-515.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,293,496,537,1218,3120,650,792.5,1650,244.584,244.585,244.6,244.624],'cm^-1')),
        HinderedRotor(inertia=(0.347076,'amu*angstrom^2'), symmetry=1, barrier=(14.7312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691455,'amu*angstrom^2'), symmetry=1, barrier=(29.3611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346943,'amu*angstrom^2'), symmetry=1, barrier=(14.7313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514586,'amu*angstrom^2'), symmetry=1, barrier=(21.8429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.484107,0.104268,-0.000158098,1.1676e-07,-3.35388e-11,-61853.3,29.093], Tmin=(100,'K'), Tmax=(858.723,'K')), NASAPolynomial(coeffs=[17.6734,0.0196888,-1.03563e-05,2.06114e-09,-1.46314e-13,-64971.8,-55.7409], Tmin=(858.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'O[C](F)[CH]F(881)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u1 p0 c0 {2,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-268.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,395,473,707,1436,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.799567,'amu*angstrom^2'), symmetry=1, barrier=(18.3836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798411,'amu*angstrom^2'), symmetry=1, barrier=(18.357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08184,0.0438232,-5.88366e-05,3.85022e-08,-9.80532e-12,-32280.3,17.0306], Tmin=(100,'K'), Tmax=(965.263,'K')), NASAPolynomial(coeffs=[10.537,0.00878523,-4.38802e-06,8.96773e-10,-6.5599e-14,-33912.6,-23.4618], Tmin=(965.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'CHFCF[Z](72)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-58.4765,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([18.6764,93.2568,111.933],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([203.406,408.59,750.559,801.84,1025.92,1150.34,1361.98,1782.01,3238.53],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9644,0.00217073,4.76141e-05,-9.63411e-08,5.92621e-11,-7031.25,9.1333], Tmin=(10,'K'), Tmax=(529.917,'K')), NASAPolynomial(coeffs=[3.26037,0.0150322,-1.01557e-05,3.21347e-09,-3.84694e-13,-7062.61,11.0829], Tmin=(529.917,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-58.4765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""F[C]DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)(F)[CH]F(4869)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6 C u1 p0 c0 {2,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-425.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(1.10704,'amu*angstrom^2'), symmetry=1, barrier=(25.453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10658,'amu*angstrom^2'), symmetry=1, barrier=(25.4424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53718,0.0522583,-6.27184e-05,3.55519e-08,-7.74699e-12,-51043.2,19.8523], Tmin=(100,'K'), Tmax=(1131.07,'K')), NASAPolynomial(coeffs=[14.1248,0.00774253,-3.68276e-06,7.55673e-10,-5.60173e-14,-53890.8,-42.4263], Tmin=(1131.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'F[CH][C](F)OC(F)=CF(4966)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u1 p0 c0 {2,S} {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u1 p0 c0 {4,S} {6,S} {11,S}
9  C u0 p0 c0 {3,S} {7,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-544.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([395,473,707,1436,326,540,652,719,1357,334,575,1197,1424,3202,194,682,905,1196,1383,3221,255.805,255.818,255.821,255.835],'cm^-1')),
        HinderedRotor(inertia=(0.345711,'amu*angstrom^2'), symmetry=1, barrier=(16.0527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187284,'amu*angstrom^2'), symmetry=1, barrier=(8.69469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837209,'amu*angstrom^2'), symmetry=1, barrier=(38.8785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0207582,0.0960651,-0.000152591,1.19426e-07,-3.58004e-11,-65321.7,27.2525], Tmin=(100,'K'), Tmax=(691.007,'K')), NASAPolynomial(coeffs=[14.2205,0.0224802,-1.2074e-05,2.39887e-09,-1.69056e-13,-67501.3,-37.7192], Tmin=(691.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-544.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[O]C(F)([CH]F)OC(F)=CF(4774)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {11,S}
10 C u0 p0 c0 {4,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-686.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,326,540,652,719,1357,334,575,1197,1424,3202,194,682,905,1196,1383,3221,234.442,234.666,234.709,234.731],'cm^-1')),
        HinderedRotor(inertia=(0.448235,'amu*angstrom^2'), symmetry=1, barrier=(17.4951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184053,'amu*angstrom^2'), symmetry=1, barrier=(7.16997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02638,'amu*angstrom^2'), symmetry=1, barrier=(40.0376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.445751,0.106514,-0.000168864,1.32667e-07,-4.00207e-11,-82382.1,30.4373], Tmin=(100,'K'), Tmax=(684.974,'K')), NASAPolynomial(coeffs=[14.985,0.0261334,-1.41482e-05,2.82379e-09,-1.9965e-13,-84724.3,-39.835], Tmin=(684.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-686.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'OC(F)([CH]F)OC(F)=[C]F(4967)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u1 p0 c0 {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-675.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,293,496,537,1218,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.894566,0.116287,-0.00019754,1.65699e-07,-5.41244e-11,-81058.2,32.2495], Tmin=(100,'K'), Tmax=(799.378,'K')), NASAPolynomial(coeffs=[16.3203,0.0241288,-1.33176e-05,2.64459e-09,-1.851e-13,-83618.1,-45.7448], Tmin=(799.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-675.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    label = 'O=C([CH]F)OC(F)=CF(4895)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,D}
6  C u0 p0 c0 {4,S} {5,D} {8,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u1 p0 c0 {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-698.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,235,1215,1347,1486,3221,194,682,905,1196,1383,3221,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406084,0.0871489,-0.000143191,1.26703e-07,-4.47414e-11,-83878,26.4986], Tmin=(100,'K'), Tmax=(771.448,'K')), NASAPolynomial(coeffs=[9.76068,0.0303378,-1.65758e-05,3.32689e-09,-2.35964e-13,-85074.1,-14.6024], Tmin=(771.448,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(CO-O2sO2d)(F1s)(H))"""),
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
    label = 'C#COC(O)(F)[CH]F(4968)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6  C u1 p0 c0 {2,S} {5,S} {9,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-383.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17664,'amu*angstrom^2'), symmetry=1, barrier=(27.0532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1763,'amu*angstrom^2'), symmetry=1, barrier=(27.0454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17685,'amu*angstrom^2'), symmetry=1, barrier=(27.058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17623,'amu*angstrom^2'), symmetry=1, barrier=(27.0437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.195043,0.0919628,-0.000127823,8.31612e-08,-2.06747e-11,-45972.6,23.9825], Tmin=(100,'K'), Tmax=(997.517,'K')), NASAPolynomial(coeffs=[20.2072,0.0101488,-4.79446e-06,9.3604e-10,-6.68536e-14,-50042.9,-74.3956], Tmin=(997.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(Ct-CtH) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'OC(F)([CH]F)OC#CF(4969)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {3,S} {8,T}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-483.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,334,575,1197,1424,3202,2175,525,239,401,1367,226.048,226.048,226.048,226.048],'cm^-1')),
        HinderedRotor(inertia=(0.574627,'amu*angstrom^2'), symmetry=1, barrier=(20.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574627,'amu*angstrom^2'), symmetry=1, barrier=(20.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574628,'amu*angstrom^2'), symmetry=1, barrier=(20.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574627,'amu*angstrom^2'), symmetry=1, barrier=(20.836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.278112,0.0991931,-0.00015038,1.08768e-07,-3.04508e-11,-57973.5,26.6116], Tmin=(100,'K'), Tmax=(881.866,'K')), NASAPolynomial(coeffs=[18.0602,0.0160144,-8.89891e-06,1.81311e-09,-1.30578e-13,-61207.9,-59.5548], Tmin=(881.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-483.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(CtCF) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H))"""),
)

species(
    label = 'O=C=CF(2952)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-172.285,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(60.0011,'amu')),
        NonlinearRotor(inertia=([9.01649,110.348,119.365],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([238.894,460.044,539.472,686.431,1048.64,1208.61,1432.92,2235.18,3235.83],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95581,0.00278106,4.95039e-05,-1.085e-07,7.1303e-11,-20718.2,7.95963], Tmin=(10,'K'), Tmax=(510.864,'K')), NASAPolynomial(coeffs=[3.75632,0.0135037,-8.87766e-06,2.78796e-09,-3.34853e-13,-20817.4,7.61808], Tmin=(510.864,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-172.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ODCDCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)([C]F)OC(F)=CF(4037)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {11,S}
10 C u0 p1 c0 {4,S} {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-794.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,326,540,652,719,1357,194,682,905,1196,1383,3221,617,898,1187,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05134,'amu*angstrom^2'), symmetry=1, barrier=(24.1725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05144,'amu*angstrom^2'), symmetry=1, barrier=(24.1746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05066,'amu*angstrom^2'), symmetry=1, barrier=(24.1567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05093,'amu*angstrom^2'), symmetry=1, barrier=(24.163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.453068,0.104527,-0.000159,1.19076e-07,-3.48862e-11,-95352.9,30.0254], Tmin=(100,'K'), Tmax=(840.208,'K')), NASAPolynomial(coeffs=[16.8494,0.0221568,-1.19506e-05,2.40292e-09,-1.71543e-13,-98260.6,-50.4371], Tmin=(840.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-794.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]C(F)(CF)OC(F)=CF(3315)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-887.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,528,1116,1182,1331,1402,1494,3075,3110,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.799291,'amu*angstrom^2'), symmetry=1, barrier=(18.3773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799993,'amu*angstrom^2'), symmetry=1, barrier=(18.3934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44171,'amu*angstrom^2'), symmetry=1, barrier=(33.1478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3798.73,'J/mol'), sigma=(5.99743,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.35 K, Pc=39.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133218,0.0970708,-0.000130171,8.87561e-08,-2.41032e-11,-106566,28.7602], Tmin=(100,'K'), Tmax=(897.966,'K')), NASAPolynomial(coeffs=[15.4002,0.0278767,-1.45861e-05,2.94312e-09,-2.12162e-13,-109356,-44.5078], Tmin=(897.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-887.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = 'OC(F)(CF)OC(F)=[C]F(4049)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-876.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,528,1116,1182,1331,1402,1494,3075,3110,293,496,537,1218,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.475818,0.105298,-0.000151794,1.09615e-07,-3.12406e-11,-105247,30.2085], Tmin=(100,'K'), Tmax=(860.146,'K')), NASAPolynomial(coeffs=[16.5133,0.0262919,-1.40161e-05,2.82883e-09,-2.0321e-13,-108169,-49.1946], Tmin=(860.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-876.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'O[C](OC(F)=CF)C(F)F(4970)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-942.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.271664,0.0999917,-0.000139358,9.71009e-08,-2.67106e-11,-113257,31.7883], Tmin=(100,'K'), Tmax=(890.124,'K')), NASAPolynomial(coeffs=[16.375,0.0251856,-1.3298e-05,2.68732e-09,-1.93653e-13,-116221,-46.5848], Tmin=(890.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-942.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_P)"""),
)

species(
    label = 'OC(F)(O[C]=CF)C(F)F(4971)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {7,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {10,D} {12,S}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-927.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,235,523,627,1123,1142,1372,1406,3097,615,860,1140,1343,3152,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.663166,0.106335,-0.000147996,9.98756e-08,-2.62077e-11,-111445,30.5372], Tmin=(100,'K'), Tmax=(938.404,'K')), NASAPolynomial(coeffs=[19.5417,0.0202069,-1.03181e-05,2.06127e-09,-1.47911e-13,-115237,-65.6541], Tmin=(938.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-F1sH)(O2s-Cs))"""),
)

species(
    label = '[CH]=C(F)OC(O)(F)C(F)F(4972)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {9,D} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-932.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,361,584,565,722,1474,235,523,627,1123,1142,1372,1406,3097,293,496,537,1218,3120,650,792.5,1650,186.634,186.68,186.689,186.773],'cm^-1')),
        HinderedRotor(inertia=(0.843881,'amu*angstrom^2'), symmetry=1, barrier=(20.8773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843943,'amu*angstrom^2'), symmetry=1, barrier=(20.8768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.844193,'amu*angstrom^2'), symmetry=1, barrier=(20.878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843594,'amu*angstrom^2'), symmetry=1, barrier=(20.877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.690008,0.107083,-0.000148108,9.91242e-08,-2.5815e-11,-112033,30.1614], Tmin=(100,'K'), Tmax=(944.764,'K')), NASAPolynomial(coeffs=[19.7538,0.0205274,-1.06855e-05,2.154e-09,-1.55409e-13,-115896,-67.3065], Tmin=(944.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-932.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
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
    E0 = (-420.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-409.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-31.0873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-87.1611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-143.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-445.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-470.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-180.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-511.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-286.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-30.6074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-18.2348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-23.1346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-285.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-64.0732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-96.3137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-54.8698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-43.9959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-372.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (32.1164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-63.7462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-407.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-163.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-154.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-382.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-412.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-336.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-317.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-320.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['O=C(F)[CH]F(215)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(105.19,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['O=C(F)C(F)C(O)(F)[CH]F(4958)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(116.002,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F(37)', '[CH]C(O)(F)OC(F)=CF(4959)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(137)', 'O[C](F)OC(F)=CF(4047)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(5)', 'OC(F)([C]F)OC(F)=CF-2(4960)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['OC1(F)O[C](F)C(F)C1F(4961)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.92372e+09,'s^-1'), n=0.369605, Ea=(79.4285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['OC1(F)OC(F)([CH]F)C1F(4962)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(54.9173,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'OC(=CF)OC(F)=CF(4055)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(31.7702,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C(F)[CH]F(215)', 'OC(F)=CF(878)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(19.9602,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(7)', 'FC=C(F)OC(F)=CF(4040)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1130.79,'m^3/(mol*s)'), n=1.28408, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.06277471317878429, var=1.0052381843076323, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'O[C]([CH]F)OC(F)=CF(4963)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -40.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'OC(F)([CH]F)O[C]=CF(4964)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -43.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[CH]=C(F)OC(O)(F)[CH]F(4965)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -42.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C(F)[CH]F(215)', 'O[C](F)[CH]F(881)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(6.03212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHFCF[Z](72)', '[O]C(O)(F)[CH]F(4869)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(7)', 'F[CH][C](F)OC(F)=CF(4966)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.54e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[O]C(F)([CH]F)OC(F)=CF(4774)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'OC(F)([CH]F)OC(F)=[C]F(4967)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -11.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'O=C([CH]F)OC(F)=CF(4895)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(187.013,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', 'C#COC(O)(F)[CH]F(4968)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(4.86643,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'OC(F)([CH]F)OC#CF(4969)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(281.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['HF(38)', 'O=C=CF(2952)', 'O=C(F)[CH]F(215)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.14547e+10,'s^-1'), n=0.631481, Ea=(117.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CHF(40)', 'O[C](F)OC(F)=CF(4047)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', 'OC(F)([C]F)OC(F)=CF(4037)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(8.51791,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(CF)OC(F)=CF(3315)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(18000,'s^-1'), n=2.287, Ea=(84.6842,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OC(F)(CF)OC(F)=[C]F(4049)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4708.23,'s^-1'), n=2.0907, Ea=(44.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_1H] for rate rule [R5H_DSSS;Cd_rad_out_single;Cs_H_out_1H]
Euclidian distance = 2.8284271247461903
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    products = ['O[C](OC(F)=CF)C(F)F(4970)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(188.527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OC(F)(O[C]=CF)C(F)F(4971)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.39758e+07,'s^-1'), n=1.42748, Ea=(190.72,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(F)OC(O)(F)C(F)F(4972)'],
    products = ['OC(F)([CH]F)OC(F)=CF(4048)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(192.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1306',
    isomers = [
        'OC(F)([CH]F)OC(F)=CF(4048)',
    ],
    reactants = [
        ('O=C(F)[CH]F(215)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1306',
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

