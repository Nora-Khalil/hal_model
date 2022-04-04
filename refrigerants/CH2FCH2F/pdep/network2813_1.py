species(
    label = 'F[CH]C(O[CH]C(F)F)=C(F)F(7556)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {10,S} {11,D}
9  C u1 p0 c0 {6,S} {7,S} {13,S}
10 C u1 p0 c0 {3,S} {8,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-806.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,350,440,435,1725,3025,407.5,1350,352.5,234,589,736,816,1240,3237,182,240,577,636,1210,1413,289.661,290.031,290.566,290.705],'cm^-1')),
        HinderedRotor(inertia=(0.215084,'amu*angstrom^2'), symmetry=1, barrier=(12.861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523079,'amu*angstrom^2'), symmetry=1, barrier=(31.3449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523561,'amu*angstrom^2'), symmetry=1, barrier=(31.354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275515,'amu*angstrom^2'), symmetry=1, barrier=(16.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02176,0.118965,-0.000176072,1.31193e-07,-3.86591e-11,-96843.3,34.1714], Tmin=(100,'K'), Tmax=(831.925,'K')), NASAPolynomial(coeffs=[17.3791,0.030489,-1.65424e-05,3.35054e-09,-2.40703e-13,-99904.9,-51.2159], Tmin=(831.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-806.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOC(O)) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
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
    label = 'FC=C=C(F)F(5948)',
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
    label = '[CH]C(F)F(506)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-68.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,200.055,1350.79,1616.28],'cm^-1')),
        HinderedRotor(inertia=(0.00421209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96601,0.019505,-1.18279e-05,2.09012e-09,3.19809e-13,-8190.14,12.8831], Tmin=(100,'K'), Tmax=(1211.05,'K')), NASAPolynomial(coeffs=[7.87164,0.00800482,-3.40878e-06,6.62037e-10,-4.73286e-14,-9723.2,-13.1469], Tmin=(1211.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=C([CH]F)[C](F)F(7925)',
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
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,235,1215,1347,1486,3221,179,346,818,1406,1524,388.405,388.787],'cm^-1')),
        HinderedRotor(inertia=(0.00112742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48239,'amu*angstrom^2'), symmetry=1, barrier=(51.7615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08221,0.0453166,-5.36923e-05,3.38512e-08,-8.69309e-12,-55175.3,19.6668], Tmin=(100,'K'), Tmax=(937.894,'K')), NASAPolynomial(coeffs=[8.71355,0.0170346,-8.45972e-06,1.69905e-09,-1.22736e-13,-56419.2,-11.9004], Tmin=(937.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(CO-CsO2d)(F1s)(H)) + radical(Csj(CO-CsO2d)(F1s)(F1s))"""),
)

species(
    label = '[CH]F(223)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'FC(F)=[C]O[CH]C(F)F(8440)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-510.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.939936,'amu*angstrom^2'), symmetry=1, barrier=(21.611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938853,'amu*angstrom^2'), symmetry=1, barrier=(21.5861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937141,'amu*angstrom^2'), symmetry=1, barrier=(21.5467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152683,0.0947474,-0.000137156,9.49642e-08,-2.54191e-11,-61280.6,27.9484], Tmin=(100,'K'), Tmax=(922.452,'K')), NASAPolynomial(coeffs=[18.17,0.0152957,-7.96037e-06,1.59336e-09,-1.14167e-13,-64661,-58.969], Tmin=(922.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOC(O)) + radical(Cdj(Cd-F1sF1s)(O2s-Cs))"""),
)

species(
    label = 'FC(F)=C1OC(C(F)F)C1F(7559)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {14,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-1094.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.919022,0.0939355,-0.000101318,5.25691e-08,-1.04021e-11,-131503,27.0517], Tmin=(100,'K'), Tmax=(1304.34,'K')), NASAPolynomial(coeffs=[24.5115,0.0124377,-3.55819e-06,5.39137e-10,-3.41338e-14,-137839,-101.248], Tmin=(1304.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1094.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCFF) + ring(2methyleneoxetane)"""),
)

species(
    label = 'FCC(OC=C(F)F)=C(F)F(7565)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {8,D}
11 C u0 p0 c0 {4,S} {5,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-942.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183926,0.104124,-0.000173348,1.62463e-07,-6.02349e-11,-113221,34.2699], Tmin=(100,'K'), Tmax=(795.585,'K')), NASAPolynomial(coeffs=[7.16105,0.0462707,-2.48193e-05,4.94347e-09,-3.48811e-13,-113728,4.67607], Tmin=(795.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-942.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'FC(F)[CH]O[C]1C(F)C1(F)F(8441)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {11,S} {13,S}
10 C u1 p0 c0 {6,S} {7,S} {8,S}
11 C u1 p0 c0 {6,S} {9,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-730.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.878867,0.121783,-0.000216587,1.98959e-07,-7.01087e-11,-87688.3,34.4294], Tmin=(100,'K'), Tmax=(841.49,'K')), NASAPolynomial(coeffs=[10.1486,0.0408659,-2.15494e-05,4.19628e-09,-2.89722e-13,-88535.2,-10.8734], Tmin=(841.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(C2CsJOCs) + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH][C]1OC(C(F)F)C1(F)F(8442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {13,S}
10 C u1 p0 c0 {6,S} {8,S} {11,S}
11 C u1 p0 c0 {5,S} {10,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-817.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.414631,0.106754,-0.000170055,1.49181e-07,-5.26125e-11,-98176.5,32.1203], Tmin=(100,'K'), Tmax=(766.982,'K')), NASAPolynomial(coeffs=[10.6805,0.0397527,-2.11503e-05,4.21858e-09,-2.98712e-13,-99609.7,-16.7117], Tmin=(766.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-817.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCsFFH) + group(CsCsFHH) + ring(Cs-O2s-Cs-Cs(F)(F)) + radical(C2CsJOCs) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)OC1C(F)F(8443)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {13,S}
10 C u1 p0 c0 {3,S} {7,S} {14,S}
11 C u1 p0 c0 {4,S} {5,S} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-790.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.554529,0.108202,-0.000158462,1.2174e-07,-3.7422e-11,-94949.2,32.9877], Tmin=(100,'K'), Tmax=(795.52,'K')), NASAPolynomial(coeffs=[14.3324,0.0333444,-1.73081e-05,3.44371e-09,-2.44816e-13,-97317.7,-35.4271], Tmin=(795.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-790.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFFH) + ring(Cs(C-FF)-Cs-O2s) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F[CH][C]=C(F)F(7343)',
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
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.805],'cm^-1')),
        HinderedRotor(inertia=(0.355657,'amu*angstrom^2'), symmetry=1, barrier=(8.17725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377451,-4.40201e-05,2.68133e-08,-6.58282e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.946,'K')), NASAPolynomial(coeffs=[8.46194,0.0130101,-6.31221e-06,1.26455e-09,-9.14158e-14,-25275.2,-12.1013], Tmin=(983.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'F[CH]C(OC=CF)=C(F)F(8444)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,D}
7  C u0 p0 c0 {5,S} {10,D} {11,S}
8  C u1 p0 c0 {1,S} {6,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 C u0 p0 c0 {4,S} {7,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-599.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,194,682,905,1196,1383,3221,180,180,1367.98],'cm^-1')),
        HinderedRotor(inertia=(0.215794,'amu*angstrom^2'), symmetry=1, barrier=(4.96153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11384,'amu*angstrom^2'), symmetry=1, barrier=(48.6013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0365845,'amu*angstrom^2'), symmetry=1, barrier=(48.6102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09996,0.0777031,-6.27283e-05,-4.97706e-08,8.90802e-11,-72061.8,29.3694], Tmin=(100,'K'), Tmax=(467.152,'K')), NASAPolynomial(coeffs=[8.35905,0.0401986,-2.14586e-05,4.29478e-09,-3.05099e-13,-73009.1,-3.00559], Tmin=(467.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
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
    label = 'F[CH]C(OC=C(F)F)=C(F)F(7525)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,D}
8  C u0 p0 c0 {6,S} {11,D} {12,S}
9  C u1 p0 c0 {1,S} {7,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {7,D}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-792.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,141,223,164,316,539,615,578,694,1133,1287,1372,1454,180,180,1380.48],'cm^-1')),
        HinderedRotor(inertia=(0.126819,'amu*angstrom^2'), symmetry=1, barrier=(2.91582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20995,'amu*angstrom^2'), symmetry=1, barrier=(50.811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0377984,'amu*angstrom^2'), symmetry=1, barrier=(50.8327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204164,0.103511,-0.000172138,1.58014e-07,-5.76318e-11,-95151.5,34.9053], Tmin=(100,'K'), Tmax=(785.979,'K')), NASAPolynomial(coeffs=[8.70032,0.0416961,-2.26832e-05,4.54103e-09,-3.21396e-13,-96041.6,-2.66755], Tmin=(785.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-792.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
)

species(
    label = '[O][CH]C(F)F(1636)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-255.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.637092,'amu*angstrom^2'), symmetry=1, barrier=(14.648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0921,0.0467254,-7.94794e-05,7.13576e-08,-2.49602e-11,-30696.7,15.5487], Tmin=(100,'K'), Tmax=(816.59,'K')), NASAPolynomial(coeffs=[6.94864,0.0153515,-7.91621e-06,1.55898e-09,-1.09043e-13,-31237,-5.34886], Tmin=(816.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = 'F[CH]C(OC[C](F)F)=C(F)F(8445)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {11,D}
9  C u1 p0 c0 {1,S} {2,S} {7,S}
10 C u1 p0 c0 {3,S} {8,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-798.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,190,488,555,1236,1407,234,589,736,816,1240,3237,182,240,577,636,1210,1413,316.522,316.531,316.534,2564.92],'cm^-1')),
        HinderedRotor(inertia=(0.165426,'amu*angstrom^2'), symmetry=1, barrier=(11.765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457021,'amu*angstrom^2'), symmetry=1, barrier=(32.4879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4569,'amu*angstrom^2'), symmetry=1, barrier=(32.487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457004,'amu*angstrom^2'), symmetry=1, barrier=(32.4868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.766134,0.114062,-0.000172703,1.3703e-07,-4.3564e-11,-95826.4,34.6122], Tmin=(100,'K'), Tmax=(769.289,'K')), NASAPolynomial(coeffs=[14.3903,0.0352477,-1.90137e-05,3.83094e-09,-2.73845e-13,-98158.1,-34.5323], Tmin=(769.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-798.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
)

species(
    label = 'FC[C](OC=C(F)F)[C](F)F(8446)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {11,D} {14,S}
10 C u1 p0 c0 {2,S} {3,S} {8,S}
11 C u0 p0 c0 {4,S} {5,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-776.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,360,370,350,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,182,240,577,636,1210,1413,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04201,0.119652,-0.000185669,1.45588e-07,-4.50569e-11,-93208.1,36.4429], Tmin=(100,'K'), Tmax=(793.843,'K')), NASAPolynomial(coeffs=[16.7598,0.0299529,-1.61783e-05,3.25065e-09,-2.31579e-13,-96034.5,-45.3307], Tmin=(793.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-776.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCsFHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C2CsJOC(O)) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(OC(F)[CH]F)=C(F)F(8447)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {10,S} {11,D}
9  C u1 p0 c0 {2,S} {7,S} {13,S}
10 C u1 p0 c0 {3,S} {8,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-791.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,350,440,435,1725,334,575,1197,1424,3202,234,589,736,816,1240,3237,182,240,577,636,1210,1413,292.79,292.799,292.82,769.238],'cm^-1')),
        HinderedRotor(inertia=(0.597754,'amu*angstrom^2'), symmetry=1, barrier=(36.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252256,'amu*angstrom^2'), symmetry=1, barrier=(15.347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597781,'amu*angstrom^2'), symmetry=1, barrier=(36.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597722,'amu*angstrom^2'), symmetry=1, barrier=(36.3647,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.605856,0.109154,-0.000146719,1.00358e-07,-2.74541e-11,-95031.4,33.3147], Tmin=(100,'K'), Tmax=(889.592,'K')), NASAPolynomial(coeffs=[16.426,0.0325727,-1.75922e-05,3.59121e-09,-2.60438e-13,-98061.8,-46.862], Tmin=(889.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-791.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
)

species(
    label = 'FC=CO[C]([C](F)F)C(F)F(8448)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 C u0 p0 c0 {5,S} {9,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-805.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0638,0.119804,-0.000182519,1.38977e-07,-4.1656e-11,-96707.5,36.2342], Tmin=(100,'K'), Tmax=(819.701,'K')), NASAPolynomial(coeffs=[17.6423,0.0285221,-1.5482e-05,3.12643e-09,-2.23767e-13,-99774.3,-50.2933], Tmin=(819.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-805.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C2CsJOC(O)) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]C(=CF)OC(F)C(F)F(8449)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {4,S} {9,D} {14,S}
11 C u2 p0 c0 {5,S} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-793.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,194,682,905,1196,1383,3221,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.352225,0.103977,-0.000136398,9.18791e-08,-2.4948e-11,-95303.1,30.9164], Tmin=(100,'K'), Tmax=(892.545,'K')), NASAPolynomial(coeffs=[15.3067,0.033802,-1.84665e-05,3.79489e-09,-2.76403e-13,-98098.5,-42.8492], Tmin=(892.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-793.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=C(O[CH]C(F)F)C(F)F(8450)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
9  C u0 p0 c0 {6,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {8,S} {14,S}
11 C u1 p0 c0 {5,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-720.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,522,611,926,1093,1137,1374,1416,3112,350,440,435,1725,3025,407.5,1350,352.5,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35726,0.129695,-0.000224233,1.96585e-07,-6.71397e-11,-86483,35.3301], Tmin=(100,'K'), Tmax=(810.389,'K')), NASAPolynomial(coeffs=[14.7279,0.0345579,-1.89979e-05,3.77698e-09,-2.64669e-13,-88573.1,-35.7002], Tmin=(810.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CCsJOC(O)) + radical(Cdj(Cd-CsO2s)(F1s))"""),
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
    label = 'FC=[C]O[CH]C(F)F(4905)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u1 p0 c0 {4,S} {5,S} {10,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-319.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,615,860,1140,1343,3152,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.864198,'amu*angstrom^2'), symmetry=1, barrier=(19.8696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864395,'amu*angstrom^2'), symmetry=1, barrier=(19.8742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864926,'amu*angstrom^2'), symmetry=1, barrier=(19.8863,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.130229,0.0871593,-0.00012213,8.22199e-08,-2.13835e-11,-38322.4,25.776], Tmin=(100,'K'), Tmax=(949.607,'K')), NASAPolynomial(coeffs=[17.4377,0.0142538,-6.96538e-06,1.36764e-09,-9.72952e-14,-41609.4,-56.8274], Tmin=(949.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CCsJOC(O)) + radical(Cdj(Cd-F1sH)(O2s-Cs))"""),
)

species(
    label = 'FC=C1OC(C(F)F)C1(F)F(7560)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {13,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {5,S} {10,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1121.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06031,0.0972683,-0.000108068,5.74016e-08,-1.15909e-11,-134676,28.0238], Tmin=(100,'K'), Tmax=(1282.41,'K')), NASAPolynomial(coeffs=[25.262,0.0113015,-2.99584e-06,4.29625e-10,-2.64998e-14,-141109,-104.275], Tmin=(1282.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1121.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCFH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'FC=C(OC=C(F)F)C(F)F(7566)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u0 p0 c0 {3,S} {8,D} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-967.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.256533,0.108118,-0.000189301,1.82511e-07,-6.83015e-11,-116216,34.0729], Tmin=(100,'K'), Tmax=(812.943,'K')), NASAPolynomial(coeffs=[5.97305,0.0486319,-2.63356e-05,5.2378e-09,-3.68135e-13,-116276,11.169], Tmin=(812.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-967.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'F[C](F)[C]1OC(C(F)F)C1F(8451)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {14,S}
10 C u1 p0 c0 {6,S} {8,S} {11,S}
11 C u1 p0 c0 {4,S} {5,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-802.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0956093,0.0985962,-0.000144212,1.18127e-07,-3.97041e-11,-96319.6,32.8226], Tmin=(100,'K'), Tmax=(723.152,'K')), NASAPolynomial(coeffs=[10.6866,0.0389623,-2.05293e-05,4.11758e-09,-2.94096e-13,-97879.2,-15.7016], Tmin=(723.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + group(CsCsFFH) + ring(Cs-O2s-Cs-Cs(F)) + radical(C2CsJOCs) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = 'F[C]C(OCC(F)F)=C(F)F(8452)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {14,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 C u2 p0 c0 {5,S} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-780.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,182,240,577,636,1210,1413,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.772798,0.114578,-0.000182773,1.55522e-07,-5.31323e-11,-93676.2,32.8391], Tmin=(100,'K'), Tmax=(740.116,'K')), NASAPolynomial(coeffs=[13.0704,0.0370413,-2.01163e-05,4.04281e-09,-2.87579e-13,-95650.8,-29.2771], Tmin=(740.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-780.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[CH][C](OC=C(F)F)C(F)F(8453)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u1 p0 c0 {3,S} {8,S} {14,S}
11 C u0 p0 c0 {4,S} {5,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-783.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,360,370,350,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34341,0.127624,-0.000213896,1.81613e-07,-6.05209e-11,-94078.4,36.0951], Tmin=(100,'K'), Tmax=(788.182,'K')), NASAPolynomial(coeffs=[16.1075,0.0316868,-1.72821e-05,3.44057e-09,-2.41946e-13,-96600.2,-42.4887], Tmin=(788.182,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCsFFH) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(C2CsJOC(O)) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH][C](OC=CF)C(F)(F)F(8454)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {11,D} {12,S}
10 C u1 p0 c0 {4,S} {8,S} {13,S}
11 C u0 p0 c0 {5,S} {9,D} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-826.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19068,0.123395,-0.000197375,1.59463e-07,-5.07897e-11,-99209.7,34.2827], Tmin=(100,'K'), Tmax=(772.249,'K')), NASAPolynomial(coeffs=[16.7829,0.0302898,-1.65131e-05,3.31491e-09,-2.35404e-13,-101985,-47.7824], Tmin=(772.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-826.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + group(CsCsFHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C2CsJOC(O)) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = '[CH]C(OC(F)C(F)F)=C(F)F(8455)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {13,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 C u2 p0 c0 {9,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-830.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425921,0.101552,-0.000114029,6.57457e-08,-1.5324e-11,-99690.4,32.6003], Tmin=(100,'K'), Tmax=(1032.16,'K')), NASAPolynomial(coeffs=[16.6792,0.0352627,-1.76938e-05,3.52303e-09,-2.52997e-13,-103221,-50.4636], Tmin=(1032.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-830.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O[CH]C(F)F)C(F)(F)F(8456)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
9  C u0 p0 c0 {6,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {8,S} {13,S}
11 C u1 p0 c0 {9,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-816.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3324,0.120154,-0.000169493,1.14164e-07,-2.96657e-11,-98065.1,33.5148], Tmin=(100,'K'), Tmax=(950.787,'K')), NASAPolynomial(coeffs=[22.8019,0.0186186,-9.30499e-06,1.84365e-09,-1.31808e-13,-102654,-81.7009], Tmin=(950.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-816.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(Cds_P)"""),
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
    E0 = (-275.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (3.31705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (235.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-267.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-250.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-44.3745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-150.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-210.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-185.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (43.0887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-49.4483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-89.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (74.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (159.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-254.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-202.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-79.1194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-76.8805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-69.9465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-2.86959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (245.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-267.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-250.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-150.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (7.56305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-204.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-219.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-87.0298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-89.6651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-53.0458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['O=CC(F)F(990)', 'FC=C=C(F)F(5948)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(F)F(506)', 'O=C([CH]F)[C](F)F(7925)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]F(223)', 'FC(F)=[C]O[CH]C(F)F(8440)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FC(F)=C1OC(C(F)F)C1F(7559)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FCC(OC=C(F)F)=C(F)F(7565)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FC(F)[CH]O[C]1C(F)C1(F)F(8441)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['F[CH][C]1OC(C(F)F)C1(F)F(8442)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(125.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['F[CH]C1([C](F)F)OC1C(F)F(8443)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(F)F(990)', 'F[CH][C]=C(F)F(7343)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(42.2569,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'F[CH]C(OC=CF)=C(F)F(8444)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.0507,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'F[CH]C(OC=C(F)F)=C(F)F(7525)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][CH]C(F)F(1636)', 'FC=C=C(F)F(5948)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][CH]C(F)F(1636)', 'F[CH][C]=C(F)F(7343)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHF(40)', 'FC(F)=[C]O[CH]C(F)F(8440)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH]C(OC[C](F)F)=C(F)F(8445)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.245e+06,'s^-1'), n=1.223, Ea=(12.1085,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 442 used for R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['FC[C](OC=C(F)F)[C](F)F(8446)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(112.053,'s^-1'), n=2.49718, Ea=(42.3838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH]C(OC(F)[CH]F)=C(F)F(8447)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(181.289,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FC=CO[C]([C](F)F)C(F)F(8448)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(198.71,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C]C(=CF)OC(F)C(F)F(8449)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(192.643,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]=C(O[CH]C(F)F)C(F)F(8450)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(186.649,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]F(744)', 'FC=[C]O[CH]C(F)F(4905)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FC=C1OC(C(F)F)C1(F)F(7560)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['FC=C(OC=C(F)F)C(F)F(7566)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['F[C](F)[C]1OC(C(F)F)C1F(8451)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(125.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CF2(43)', 'FC=[C]O[CH]C(F)F(4905)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]C(OCC(F)F)=C(F)F(8452)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH][C](OC=C(F)F)C(F)F(8453)'],
    products = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(10500,'s^-1'), n=2.14, Ea=(33.3465,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_noH] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['F[CH][C](OC=CF)C(F)(F)F(8454)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(188.561,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['[CH]C(OC(F)C(F)F)=C(F)F(8455)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(185.926,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(O[CH]C(F)F)=C(F)F(7556)'],
    products = ['[CH]=C(O[CH]C(F)F)C(F)(F)F(8456)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(222.545,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #2813',
    isomers = [
        'F[CH]C(O[CH]C(F)F)=C(F)F(7556)',
    ],
    reactants = [
        ('O=CC(F)F(990)', 'FC=C=C(F)F(5948)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2813',
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

