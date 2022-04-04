species(
    label = '[O]C(F)(C=O)C([CH]F)=C(F)F(5250)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u1 p0 c0 {2,S} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-721.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,234,589,736,816,1240,3237,182,240,577,636,1210,1413,217.854,218.41],'cm^-1')),
        HinderedRotor(inertia=(1.11572,'amu*angstrom^2'), symmetry=1, barrier=(37.6135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193346,'amu*angstrom^2'), symmetry=1, barrier=(6.5176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12022,'amu*angstrom^2'), symmetry=1, barrier=(37.6194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.499054,0.108072,-0.000176202,1.51682e-07,-5.17508e-11,-86608.6,34.6731], Tmin=(100,'K'), Tmax=(787.667,'K')), NASAPolynomial(coeffs=[12.4292,0.0336627,-1.7827e-05,3.52506e-09,-2.47606e-13,-88373.6,-22.8886], Tmin=(787.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-721.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = '[O]C(F)C([CH]F)=C(F)F(5521)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-613.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00285213,'amu*angstrom^2'), symmetry=1, barrier=(3.99289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.34005,'amu*angstrom^2'), symmetry=1, barrier=(53.8023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986169,0.0711419,-8.73976e-05,5.64251e-08,-1.47906e-11,-73719.1,27.1352], Tmin=(100,'K'), Tmax=(921.19,'K')), NASAPolynomial(coeffs=[11.501,0.0254841,-1.30516e-05,2.62071e-09,-1.88773e-13,-75656.3,-22.7298], Tmin=(921.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFHO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(O2sj(Cs-F1sCdH)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[O]C(F)(C=O)C(F)[C]=C(F)F(5254)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-597.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2782.5,750,1395,475,1775,1000,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.74735,'amu*angstrom^2'), symmetry=1, barrier=(17.1831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0764122,'amu*angstrom^2'), symmetry=1, barrier=(1.75687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80665,'amu*angstrom^2'), symmetry=1, barrier=(41.5384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4166.31,'J/mol'), sigma=(6.35837,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=650.77 K, Pc=36.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.191156,0.10176,-0.000165046,1.46684e-07,-5.24037e-11,-71735.4,35.4356], Tmin=(100,'K'), Tmax=(759.398,'K')), NASAPolynomial(coeffs=[10.2373,0.0377429,-2.06465e-05,4.15979e-09,-2.9613e-13,-73057.3,-10.2804], Tmin=(759.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C=C([O])C(=CF)C(F)(F)F(5522)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {4,S} {8,D} {12,S}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-921.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09645,0.104838,-0.000131316,7.54645e-08,-1.59558e-11,-110645,32.1235], Tmin=(100,'K'), Tmax=(977.339,'K')), NASAPolynomial(coeffs=[25.3868,0.00867492,-2.49136e-06,3.89862e-10,-2.60582e-14,-116405,-98.023], Tmin=(977.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-921.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCdFFF) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(CdCFH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=CC(F)(F)[C]([CH]F)C(=O)F(5523)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {3,S} {8,S} {13,S}
11 C u0 p0 c0 {4,S} {6,D} {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-831.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.955961,0.125556,-0.000236228,2.1799e-07,-7.54872e-11,-99842.6,34.4758], Tmin=(100,'K'), Tmax=(872.185,'K')), NASAPolynomial(coeffs=[10.6987,0.0355421,-1.85384e-05,3.53074e-09,-2.38159e-13,-100485,-12.1847], Tmin=(872.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-831.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + group(COCsFO) + group(Cds-OdCsH) + radical(C2CJCHO) + radical(CsCsF1sH)"""),
)

species(
    label = '[O]C(=COF)C([CH]F)=C(F)F(5524)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {8,S} {10,S} {11,D}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 C u1 p0 c0 {1,S} {7,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {7,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-498.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,295.805,296.925,297.849],'cm^-1')),
        HinderedRotor(inertia=(0.522045,'amu*angstrom^2'), symmetry=1, barrier=(33.0178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282216,'amu*angstrom^2'), symmetry=1, barrier=(17.6662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279464,'amu*angstrom^2'), symmetry=1, barrier=(17.6432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09353,0.113575,-0.000161869,1.10594e-07,-2.90446e-11,-59785.9,33.268], Tmin=(100,'K'), Tmax=(943.207,'K')), NASAPolynomial(coeffs=[21.8189,0.0164068,-7.33912e-06,1.37093e-09,-9.45319e-14,-64108.1,-75.9314], Tmin=(943.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-498.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(CdCFF) + radical(C=C(C)OJ) + radical(CsCdF1sH)"""),
)

species(
    label = 'O=C(F)[CH]OC([CH]F)=C(F)F(5251)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {9,S} {10,D}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  C u1 p0 c0 {1,S} {7,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {7,D}
11 C u0 p0 c0 {4,S} {6,D} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-716.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,234,589,736,816,1240,3237,182,240,577,636,1210,1413,611,648,830,1210,1753,391.741,392.288,392.358,393.349],'cm^-1')),
        HinderedRotor(inertia=(0.212193,'amu*angstrom^2'), symmetry=1, barrier=(23.1356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255337,'amu*angstrom^2'), symmetry=1, barrier=(27.8529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0785174,'amu*angstrom^2'), symmetry=1, barrier=(8.59226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255245,'amu*angstrom^2'), symmetry=1, barrier=(27.8468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3740.69,'J/mol'), sigma=(5.63003,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.29 K, Pc=47.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.59798,0.107743,-0.000149408,1.02142e-07,-2.75242e-11,-86061.1,34.0598], Tmin=(100,'K'), Tmax=(908.245,'K')), NASAPolynomial(coeffs=[17.9825,0.025908,-1.42473e-05,2.92613e-09,-2.12855e-13,-89436,-53.7912], Tmin=(908.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(COCsFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOC(O)) + radical(Csj(Cd-O2sCd)(F1s)(H))"""),
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
    label = 'O=CC(F)=C([CH]F)[C](F)F(5525)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,D} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {6,D} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {5,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-605.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,288,410,724,839,1320,234,589,736,816,1240,3237,161,297,490,584,780,1358,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00746976,'amu*angstrom^2'), symmetry=1, barrier=(7.07068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74546,'amu*angstrom^2'), symmetry=1, barrier=(40.1316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0424643,'amu*angstrom^2'), symmetry=1, barrier=(40.1374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0519,0.0726256,-7.88308e-05,4.54742e-08,-1.10685e-11,-72760.9,29.8096], Tmin=(100,'K'), Tmax=(964.819,'K')), NASAPolynomial(coeffs=[10.3999,0.0338708,-1.858e-05,3.84311e-09,-2.81445e-13,-74564.8,-14.9544], Tmin=(964.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-605.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)H) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
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
    label = '[O]C(F)([C]=C(F)F)C=O(5526)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {5,D} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-414.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.703965,'amu*angstrom^2'), symmetry=1, barrier=(16.1855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08131,'amu*angstrom^2'), symmetry=1, barrier=(24.8614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504084,0.086968,-0.000156682,1.43109e-07,-5.00589e-11,-49766.9,27.9875], Tmin=(100,'K'), Tmax=(840.854,'K')), NASAPolynomial(coeffs=[9.13991,0.0265112,-1.42691e-05,2.79355e-09,-1.93196e-13,-50534.2,-8.10544], Tmin=(840.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = 'O=CC1(F)OC(F)C1=C(F)F(5262)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-983.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452607,0.075567,-6.69761e-05,2.81387e-08,-4.70781e-12,-118175,25.4586], Tmin=(100,'K'), Tmax=(1411.68,'K')), NASAPolynomial(coeffs=[18.4003,0.024713,-1.29414e-05,2.62119e-09,-1.88888e-13,-123242,-67.3168], Tmin=(1411.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-983.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCFHO) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + ring(Cyclobutane)"""),
)

species(
    label = '[O]C(F)(C=O)[C]1C(F)C1(F)F(5527)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 C u0 p0 c0 {6,D} {9,S} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-661.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483308,0.0817396,-8.97182e-05,5.1197e-08,-1.19175e-11,-79446.1,30.5536], Tmin=(100,'K'), Tmax=(1027.88,'K')), NASAPolynomial(coeffs=[13.534,0.0309517,-1.56013e-05,3.12504e-09,-2.2524e-13,-82129,-32.7672], Tmin=(1027.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-661.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(C=OCOJ) + radical(CCJ(C)CO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = '[O][CH]C1(F)OC(F)(F)C1=CF(5294)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-690.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.131797,0.100169,-0.000144973,1.13701e-07,-3.64435e-11,-82862.6,28.3647], Tmin=(100,'K'), Tmax=(757.453,'K')), NASAPolynomial(coeffs=[11.7643,0.0373477,-2.05661e-05,4.2047e-09,-3.03959e-13,-84664.8,-25.7225], Tmin=(757.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C1(F)[CH]OO1(5528)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u1 p0 c0 {6,S} {7,S} {12,S}
10 C u1 p0 c0 {2,S} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-464.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184817,0.0916803,-0.000123824,9.11468e-08,-2.74833e-11,-55719.3,31.6218], Tmin=(100,'K'), Tmax=(803.085,'K')), NASAPolynomial(coeffs=[11.5011,0.035311,-1.85285e-05,3.72966e-09,-2.67923e-13,-57536.7,-20.4903], Tmin=(803.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCsJOO) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)C1=C(F)F(5400)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {5,S} {7,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-733.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.284064,0.0869293,-8.97924e-05,4.43062e-08,-8.50023e-12,-88047.6,27.3436], Tmin=(100,'K'), Tmax=(1272.05,'K')), NASAPolynomial(coeffs=[21.6581,0.0179315,-8.43011e-06,1.66513e-09,-1.19829e-13,-93629.9,-83.7945], Tmin=(1272.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
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
    label = '[O]C1OC1(F)C([CH]F)=C(F)F(5529)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {12,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 C u1 p0 c0 {2,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-659.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558249,0.0863715,-0.000107439,6.31223e-08,-7.2672e-12,-79195.4,30.2628], Tmin=(100,'K'), Tmax=(558.91,'K')), NASAPolynomial(coeffs=[8.76197,0.041125,-2.21457e-05,4.49151e-09,-3.23407e-13,-80322.8,-6.4243], Tmin=(558.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-659.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[O]C1C(F)C(=C(F)F)C1([O])F(5453)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-606.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296682,0.0898079,-0.000110194,7.15057e-08,-1.91125e-11,-72775.6,25.1136], Tmin=(100,'K'), Tmax=(896.76,'K')), NASAPolynomial(coeffs=[12.469,0.0355142,-1.93788e-05,3.99332e-09,-2.9162e-13,-74958.8,-32.2845], Tmin=(896.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = 'O=CC(=O)C([CH]F)=C(F)F(5530)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,D}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u1 p0 c0 {1,S} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 C u0 p0 c0 {5,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-653.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,375,552.5,462.5,1710,234,589,736,816,1240,3237,182,240,577,636,1210,1413,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0467442,'amu*angstrom^2'), symmetry=1, barrier=(2.23225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.129,'amu*angstrom^2'), symmetry=1, barrier=(48.9498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218804,'amu*angstrom^2'), symmetry=1, barrier=(48.9903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798578,0.0806441,-0.000101655,5.46071e-08,-1.06079e-12,-78524,28.609], Tmin=(100,'K'), Tmax=(553.437,'K')), NASAPolynomial(coeffs=[9.1951,0.0349951,-1.86869e-05,3.75907e-09,-2.68753e-13,-79683.7,-9.01247], Tmin=(553.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-653.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(CdCFF) + group(Cds-O2d(Cds-O2d)H) + radical(CsCdF1sH)"""),
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
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C(F)C([CH]F)=C(F)F(3254)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,D} {9,S}
7  C u1 p0 c0 {1,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-791.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,255,533,799,832,1228,620.783],'cm^-1')),
        HinderedRotor(inertia=(0.154492,'amu*angstrom^2'), symmetry=1, barrier=(3.55208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0129993,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979135,0.0725884,-0.000107444,8.59818e-08,-2.79109e-11,-95063.1,25.3328], Tmin=(100,'K'), Tmax=(750.717,'K')), NASAPolynomial(coeffs=[9.78045,0.0256917,-1.3738e-05,2.76507e-09,-1.97772e-13,-96384.5,-14.6047], Tmin=(750.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-791.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cd-CdCs(CO)) + group(COCFO) + group(CdCFF) + radical(CsCdF1sH)"""),
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
    label = 'O=[C]C(=O)C([CH]F)=C(F)F(5531)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {7,S} {8,S} {9,D}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u1 p0 c0 {1,S} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-493.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,375,552.5,462.5,1710,234,589,736,816,1240,3237,182,240,577,636,1210,1413,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.0011895,'amu*angstrom^2'), symmetry=1, barrier=(4.25494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84606,'amu*angstrom^2'), symmetry=1, barrier=(42.4446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85045,'amu*angstrom^2'), symmetry=1, barrier=(42.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251504,0.0906143,-0.000150066,1.31393e-07,-4.55144e-11,-59263.5,31.3128], Tmin=(100,'K'), Tmax=(790.27,'K')), NASAPolynomial(coeffs=[10.5631,0.0291145,-1.56688e-05,3.11345e-09,-2.19137e-13,-60602.6,-14.1688], Tmin=(790.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(CdCFF) + group(Cds-O2d(Cds-O2d)H) + radical(CsCdF1sH) + radical(CCCJ=O)"""),
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
    label = 'O=[C]C(O)(F)C([CH]F)=C(F)F(5532)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-805.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.779906,0.113531,-0.000185677,1.54989e-07,-5.08943e-11,-96673.2,35.7217], Tmin=(100,'K'), Tmax=(788.619,'K')), NASAPolynomial(coeffs=[14.9553,0.0288138,-1.52083e-05,2.99319e-09,-2.09376e-13,-99002.5,-35.4876], Tmin=(788.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-805.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)C(CF)=C(F)F(5533)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-707.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,182,240,577,636,1210,1413,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13796,'amu*angstrom^2'), symmetry=1, barrier=(26.164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542745,'amu*angstrom^2'), symmetry=1, barrier=(12.4788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5452,'amu*angstrom^2'), symmetry=1, barrier=(12.5352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729701,0.115724,-0.00020266,1.80179e-07,-6.17903e-11,-84896,35.1666], Tmin=(100,'K'), Tmax=(835.773,'K')), NASAPolynomial(coeffs=[12.3728,0.0327564,-1.73947e-05,3.39665e-09,-2.34921e-13,-86378.6,-21.4617], Tmin=(835.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-707.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=C(OF)C([CH]F)=C(F)F(5534)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,D}
8  C u0 p0 c0 {5,S} {7,S} {11,D}
9  C u1 p0 c0 {1,S} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {7,D}
11 C u0 p0 c0 {6,S} {8,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-494.945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,325,375,415,465,420,450,1700,1750,234,589,736,816,1240,3237,182,240,577,636,1210,1413,3010,987.5,1337.5,450,1655,280.994,281.02,281.095],'cm^-1')),
        HinderedRotor(inertia=(0.470582,'amu*angstrom^2'), symmetry=1, barrier=(26.3473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.469682,'amu*angstrom^2'), symmetry=1, barrier=(26.3514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.469797,'amu*angstrom^2'), symmetry=1, barrier=(26.3458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.956558,0.109476,-0.000146972,9.43472e-08,-2.3352e-11,-59349.9,32.6868], Tmin=(100,'K'), Tmax=(997.948,'K')), NASAPolynomial(coeffs=[22.0912,0.017094,-8.1116e-06,1.58207e-09,-1.12667e-13,-63949.9,-78.4577], Tmin=(997.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(CdCFF) + radical(C=COJ) + radical(CsCdF1sH)"""),
)

species(
    label = '[O]C=C([O])C(=C(F)F)C(F)F(5535)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u0 p0 c0 {6,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-873.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.981952,0.105473,-0.000137096,8.43022e-08,-1.97898e-11,-104920,33.1164], Tmin=(100,'K'), Tmax=(1057.31,'K')), NASAPolynomial(coeffs=[23.6721,0.0121992,-4.76477e-06,8.6076e-10,-5.94994e-14,-110133,-87.1983], Tmin=(1057.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-873.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(CdCFF) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=CC(F)(OF)C([C]F)=CF(5536)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {8,D} {13,S}
11 C u2 p0 c0 {4,S} {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.92624,0.120032,-0.00020884,1.85912e-07,-6.4567e-11,-48516.3,33.0524], Tmin=(100,'K'), Tmax=(806.012,'K')), NASAPolynomial(coeffs=[13.0022,0.0345224,-1.92089e-05,3.83981e-09,-2.70068e-13,-50229.3,-27.8383], Tmin=(806.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C(F)(C=O)C(=[C]F)C(F)F(5537)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {13,S}
11 C u1 p0 c0 {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-634.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2782.5,750,1395,475,1775,1000,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.522021,'amu*angstrom^2'), symmetry=1, barrier=(12.0023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521864,'amu*angstrom^2'), symmetry=1, barrier=(11.9987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46829,'amu*angstrom^2'), symmetry=1, barrier=(33.7589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.743534,0.117111,-0.000208025,1.88286e-07,-6.5655e-11,-76110.3,34.9933], Tmin=(100,'K'), Tmax=(832.197,'K')), NASAPolynomial(coeffs=[11.5184,0.0352178,-1.90409e-05,3.74592e-09,-2.60235e-13,-77356.3,-17.1355], Tmin=(832.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-634.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = '[O]C(F)(C=O)C(F)(F)[C]=CF(5258)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-622.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39522,'amu*angstrom^2'), symmetry=1, barrier=(32.0789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674921,'amu*angstrom^2'), symmetry=1, barrier=(15.5178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674734,'amu*angstrom^2'), symmetry=1, barrier=(15.5135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4155.65,'J/mol'), sigma=(6.65753,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=649.10 K, Pc=31.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.62062,0.111493,-0.000186417,1.62487e-07,-5.58116e-11,-74732.3,35.2404], Tmin=(100,'K'), Tmax=(793.325,'K')), NASAPolynomial(coeffs=[12.6953,0.0333978,-1.80428e-05,3.58707e-09,-2.5232e-13,-76500.4,-23.7455], Tmin=(793.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'O=CC(F)[C]([C](F)F)C(=O)F(5538)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u1 p0 c0 {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {13,S}
10 C u1 p0 c0 {2,S} {3,S} {8,S}
11 C u0 p0 c0 {4,S} {6,D} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-830.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.827148,0.121959,-0.000226694,2.09225e-07,-7.28011e-11,-99706.5,35.2606], Tmin=(100,'K'), Tmax=(865.471,'K')), NASAPolynomial(coeffs=[10.3281,0.0362565,-1.89779e-05,3.63655e-09,-2.46877e-13,-100359,-9.55717], Tmin=(865.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-830.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(COCsFO) + group(Cds-OdCsH) + radical(C2CJCHO) + radical(CsCsF1sF1s)"""),
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
    label = '[O]C(F)([C]=CF)C=O(5077)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {4,D} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-213.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0062,'amu*angstrom^2'), symmetry=1, barrier=(23.1346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788607,'amu*angstrom^2'), symmetry=1, barrier=(18.1316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.859501,0.0775521,-0.000134509,1.21601e-07,-4.24638e-11,-25599.7,25.7054], Tmin=(100,'K'), Tmax=(836.396,'K')), NASAPolynomial(coeffs=[8.40454,0.025749,-1.34133e-05,2.60803e-09,-1.80242e-13,-26312,-6.06014], Tmin=(836.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = 'O=CC1(F)OC(F)(F)C1=CF(5264)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1017.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.408671,0.0785493,-7.41011e-05,3.36164e-08,-6.09762e-12,-122225,26.2022], Tmin=(100,'K'), Tmax=(1307.65,'K')), NASAPolynomial(coeffs=[17.6482,0.0258147,-1.36092e-05,2.77619e-09,-2.01473e-13,-126733,-61.5924], Tmin=(1307.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1017.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + ring(Cyclobutane)"""),
)

species(
    label = '[O][CH]C1(F)OC(F)C1=C(F)F(5389)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-656.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.111103,0.0947693,-0.000129077,9.63818e-08,-2.97578e-11,-78820.9,26.9126], Tmin=(100,'K'), Tmax=(781.453,'K')), NASAPolynomial(coeffs=[11.0917,0.0385649,-2.1195e-05,4.3493e-09,-3.15816e-13,-80537.2,-23.3548], Tmin=(781.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)(F)C1=CF(5304)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {5,S} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-767.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.347882,0.090145,-9.76937e-05,5.0688e-08,-1.02214e-11,-92097,28.158], Tmin=(100,'K'), Tmax=(1212.65,'K')), NASAPolynomial(coeffs=[21.4863,0.0181221,-8.60209e-06,1.7079e-09,-1.23405e-13,-97392.3,-81.3885], Tmin=(1212.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-767.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1C([O])(F)C(=CF)C1(F)F(5353)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {4,S} {10,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-632.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117869,0.0935693,-0.000118376,7.80947e-08,-2.1006e-11,-75946.5,26.2213], Tmin=(100,'K'), Tmax=(895.78,'K')), NASAPolynomial(coeffs=[13.4612,0.0339843,-1.85961e-05,3.83292e-09,-2.79836e-13,-78337,-36.6833], Tmin=(895.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C(F)[CH]O(4283)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-466.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9486,0.00327696,8.62766e-05,-1.78856e-07,1.15095e-10,-56066.9,9.87278], Tmin=(10,'K'), Tmax=(486.336,'K')), NASAPolynomial(coeffs=[2.12805,0.0278262,-1.89743e-05,5.9068e-09,-6.93439e-13,-56003.1,16.1793], Tmin=(486.336,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-466.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)[CH]O""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)([C]=O)C(=CF)C(F)F(5539)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-724.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750394,0.115942,-0.00020239,1.79176e-07,-6.12416e-11,-86991.8,35.3484], Tmin=(100,'K'), Tmax=(834.597,'K')), NASAPolynomial(coeffs=[12.6631,0.0322583,-1.7128e-05,3.34524e-09,-2.31431e-13,-88555.3,-22.892], Tmin=(834.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-724.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC(O)(F)C([C]F)=C(F)F(5540)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {8,D}
11 C u2 p0 c0 {4,S} {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-740.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.788961,0.116882,-0.000203531,1.81388e-07,-6.28924e-11,-88877.1,33.3956], Tmin=(100,'K'), Tmax=(814.083,'K')), NASAPolynomial(coeffs=[12.5478,0.0339427,-1.86323e-05,3.7023e-09,-2.59346e-13,-90471.7,-24.6602], Tmin=(814.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-740.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]C(=C(F)F)C(F)(C=O)OF(5541)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {8,D}
11 C u2 p0 c0 {8,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-430.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,182,240,577,636,1210,1413,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.535213,0.110831,-0.000155556,1.07478e-07,-2.49117e-11,-51586.9,32.5574], Tmin=(100,'K'), Tmax=(617.34,'K')), NASAPolynomial(coeffs=[13.0808,0.0388688,-2.02139e-05,3.98933e-09,-2.81437e-13,-53577.9,-29.0741], Tmin=(617.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-430.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(CdCFF) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C(F)(F)F)C([O])(F)C=O(5542)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 C u1 p0 c0 {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-694.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,219,296,586,564,718,793,1177,1228,350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.648247,'amu*angstrom^2'), symmetry=1, barrier=(14.9045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64896,'amu*angstrom^2'), symmetry=1, barrier=(14.9209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28258,'amu*angstrom^2'), symmetry=1, barrier=(29.489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.682768,0.112444,-0.000187617,1.61351e-07,-5.44289e-11,-83355.3,33.4421], Tmin=(100,'K'), Tmax=(805.619,'K')), NASAPolynomial(coeffs=[13.4786,0.0313701,-1.66269e-05,3.27249e-09,-2.28582e-13,-85287.8,-29.6502], Tmin=(805.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-694.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
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
    E0 = (-233.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (93.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-15.0491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-15.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-121.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (197.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-70.5621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (125.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (288.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-224.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-2.06603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-107.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (23.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-179.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-186.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-130.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-118.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-70.3499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-166.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-154.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-76.6597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (72.9539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (4.90044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (212.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-158.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-113.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (97.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-74.6134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (159.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (39.9292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-8.67124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-120.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (308.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-224.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-107.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-179.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-144.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (8.08611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (70.6643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-127.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (186.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (140.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (7.83045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C([CH]F)=C(F)F(5521)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(88.9,'m^3/(mol*s)'), n=1.51, Ea=(338.022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(C=O)C(F)[C]=C(F)F(5254)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C=C([O])C(=CF)C(F)(F)F(5522)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(218.007,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC(F)(F)[C]([CH]F)C(=O)F(5523)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(112.132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=COF)C([CH]F)=C(F)F(5524)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(208.297,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=C(F)[CH]OC([CH]F)=C(F)F(5251)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(158.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', 'O=CC(F)=C([CH]F)[C](F)F(5525)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]F(804)', '[O]C(F)([C]=C(F)F)C=O(5526)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC1(F)OC(F)C1=C(F)F(5262)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C(F)(C=O)[C]1C(F)C1(F)F(5527)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O][CH]C1(F)OC(F)(F)C1=CF(5294)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]OO1(5528)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C1(F)[CH]OC(F)C1=C(F)F(5400)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.47116e+08,'s^-1'), n=0.669085, Ea=(53.4924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC1(F)OC1([CH]F)[C](F)F(5474)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C1OC1(F)C([CH]F)=C(F)F(5529)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C1C(F)C(=C(F)F)C1([O])F(5453)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(115.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 113.9 to 115.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'O=CC(=O)C([CH]F)=C(F)F(5530)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(22.4238,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CC(=O)F(4234)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(29.6835,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HCO(15)', 'O=C(F)C([CH]F)=C(F)F(3254)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(116.128,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][CH]C(=O)F(398)', 'FC=C=C(F)F(1325)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(14.4224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][CH]C(=O)F(398)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=[C]C(=O)C([CH]F)=C(F)F(5531)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(291.719,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CHF(40)', '[O]C(F)([C]=C(F)F)C=O(5526)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=[C]C(O)(F)C([CH]F)=C(F)F(5532)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(F)([C]=O)C(CF)=C(F)F(5533)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_1H] for rate rule [R4H_SS(Cd)S;CO_rad_out;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=C(OF)C([CH]F)=C(F)F(5534)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(103.954,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C=C([O])C(=C(F)F)C(F)F(5535)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(158.669,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=CC(F)(OF)C([C]F)=CF(5536)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(75.7412,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(F)(C=O)C(=[C]F)C(F)F(5537)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(C=O)C(F)(F)[C]=CF(5258)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC(F)[C]([C](F)F)C(=O)F(5538)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(112.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[C]F(138)', '[O]C(F)([C]=CF)C=O(5077)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC1(F)OC(F)(F)C1=CF(5264)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O][CH]C1(F)OC(F)C1=C(F)F(5389)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C1(F)[CH]OC(F)(F)C1=CF(5304)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.47116e+08,'s^-1'), n=0.669085, Ea=(53.4924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C1C([O])(F)C(=CF)C1(F)F(5353)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(88.819,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 88.0 to 88.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=C(F)[CH]O(4283)', 'F[C]=C=C(F)F(2804)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.01084e+12,'s^-1'), n=0.0883205, Ea=(241.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.021492990850914453, var=5.303057773824753, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4R!H-R_7BrCClFIOPSSi->O',), comment="""Estimated from node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4R!H-R_7BrCClFIOPSSi->O"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CF2(43)', '[O]C(F)([C]=CF)C=O(5077)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['[O]C(F)([C]=O)C(=CF)C(F)F(5539)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=2.04, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_S(Cd)SS;C_rad_out_single;XH_out] for rate rule [R4H_S(Cd)SS;C_rad_out_noH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    products = ['O=CC(O)(F)C([C]F)=C(F)F(5540)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_single]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C(=C(F)F)C(F)(C=O)OF(5541)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(82.8722,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C(C(F)(F)F)C([O])(F)C=O(5542)'],
    products = ['[O]C(F)(C=O)C([CH]F)=C(F)F(5250)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #1529',
    isomers = [
        '[O]C(F)(C=O)C([CH]F)=C(F)F(5250)',
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
    label = 'PDepNetwork #1529',
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

