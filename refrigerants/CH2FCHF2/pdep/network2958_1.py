species(
    label = '[O]C(C(=O)F)C(F)[C]=C(F)F(11539)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-633.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,486,617,768,1157,1926,562,600,623,1070,1265,1685,370,247.317,247.472,1254.51,1728.56],'cm^-1')),
        HinderedRotor(inertia=(0.00275407,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14596,'amu*angstrom^2'), symmetry=1, barrier=(6.33926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14595,'amu*angstrom^2'), symmetry=1, barrier=(6.33914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.531841,0.0844646,-0.000110676,8.01505e-08,-2.41253e-11,-76029.1,37.3543], Tmin=(100,'K'), Tmax=(799.075,'K')), NASAPolynomial(coeffs=[10.2512,0.035812,-1.93473e-05,3.95563e-09,-2.87019e-13,-77582.4,-7.35589], Tmin=(799.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-633.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'O=CC(=O)F(335)',
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
    label = 'FC=C=C(F)F(5101)',
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
    label = '[O]C(F)C(F)[C]=C(F)F(9746)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-490.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([164,312,561,654,898,1207,1299,3167,391,562,707,872,1109,1210,1289,3137,562,600,623,1070,1265,1685,370,180,1565.64],'cm^-1')),
        HinderedRotor(inertia=(0.283333,'amu*angstrom^2'), symmetry=1, barrier=(6.51438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283186,'amu*angstrom^2'), symmetry=1, barrier=(6.51101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3621.61,'J/mol'), sigma=(5.63348,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.69 K, Pc=45.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14028,0.0694917,-9.55762e-05,7.55495e-08,-2.49933e-11,-58919.3,28.5092], Tmin=(100,'K'), Tmax=(729.007,'K')), NASAPolynomial(coeffs=[8.13024,0.0311379,-1.66588e-05,3.37977e-09,-2.43675e-13,-59938.5,-3.0039], Tmin=(729.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-490.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(O2sj(Cs-CsF1sH)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C(C(=O)F)C([CH]F)=C(F)F(9286)',
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
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,486,617,768,1157,1926,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,626.816,627.586],'cm^-1')),
        HinderedRotor(inertia=(0.0114973,'amu*angstrom^2'), symmetry=1, barrier=(3.25819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14236,'amu*angstrom^2'), symmetry=1, barrier=(3.27313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011821,'amu*angstrom^2'), symmetry=1, barrier=(3.30363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269151,0.0881154,-0.000112605,7.54746e-08,-2.04419e-11,-89137.5,36.7913], Tmin=(100,'K'), Tmax=(894.936,'K')), NASAPolynomial(coeffs=[13.2379,0.0301492,-1.54461e-05,3.09615e-09,-2.22575e-13,-91458.7,-24.3354], Tmin=(894.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'O=C[C](F)OC(F)[C]=C(F)F(11545)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {5,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-679.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,280,501,1494,1531,2782.5,750,1395,475,1775,1000,562,600,623,1070,1265,1685,370,180,180,180,2098.84],'cm^-1')),
        HinderedRotor(inertia=(0.568524,'amu*angstrom^2'), symmetry=1, barrier=(13.0715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95824,'amu*angstrom^2'), symmetry=1, barrier=(45.0237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95675,'amu*angstrom^2'), symmetry=1, barrier=(44.9896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95823,'amu*angstrom^2'), symmetry=1, barrier=(45.0235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3700.2,'J/mol'), sigma=(5.66624,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=577.96 K, Pc=46.15 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.114716,0.0986773,-0.000146825,1.19235e-07,-3.93299e-11,-81574.8,33.3191], Tmin=(100,'K'), Tmax=(739.175,'K')), NASAPolynomial(coeffs=[11.5374,0.0356104,-1.88188e-05,3.76295e-09,-2.67942e-13,-83297,-19.3717], Tmin=(739.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(CsCOF1sO2s) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
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
    label = 'O=C(F)[CH]C(F)[C]=C(F)F(11811)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {6,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-544.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,3025,407.5,1350,352.5,611,648,830,1210,1753,562,600,623,1070,1265,1685,370,180,180,814.467],'cm^-1')),
        HinderedRotor(inertia=(0.184967,'amu*angstrom^2'), symmetry=1, barrier=(4.25276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185149,'amu*angstrom^2'), symmetry=1, barrier=(4.25694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104569,'amu*angstrom^2'), symmetry=1, barrier=(49.2599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847838,0.0772178,-0.000101363,7.57968e-08,-2.39008e-11,-65423,31.0092], Tmin=(100,'K'), Tmax=(760.276,'K')), NASAPolynomial(coeffs=[8.73728,0.0357107,-1.94736e-05,3.99169e-09,-2.89937e-13,-66622.7,-4.89073], Tmin=(760.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-544.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(CCJC=O) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[C]=C(F)F(3536)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u2 p0 c0 {3,D}
"""),
    E0 = (195.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([182,240,577,636,1210,1413],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41427e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.13,'K')), NASAPolynomial(coeffs=[7.28631,0.0032378,-1.64877e-06,3.54597e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[O]C([CH]F)C(=O)F(1655)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {5,S} {9,S}
7 C u0 p0 c0 {2,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-377.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,334,575,1197,1424,3202,486,617,768,1157,1926,255.4,997.392,4000],'cm^-1')),
        HinderedRotor(inertia=(0.238295,'amu*angstrom^2'), symmetry=1, barrier=(11.0764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598441,'amu*angstrom^2'), symmetry=1, barrier=(27.8119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78051,0.052316,-6.64848e-05,4.47469e-08,-1.21873e-11,-45272.2,25.505], Tmin=(100,'K'), Tmax=(889.955,'K')), NASAPolynomial(coeffs=[9.34319,0.0183245,-9.1925e-06,1.82886e-09,-1.30954e-13,-46618.2,-10.0988], Tmin=(889.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFHH) + group(COCsFO) + radical(C=OCOJ) + radical(CsCsF1sH)"""),
)

species(
    label = 'O=C(F)C1OC(=C(F)F)C1F(11547)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-1018.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0115179,0.0724509,-4.27473e-05,-1.2389e-08,1.41852e-11,-122378,27.3941], Tmin=(100,'K'), Tmax=(949.287,'K')), NASAPolynomial(coeffs=[23.7091,0.00972333,-2.29406e-06,4.00602e-10,-3.29915e-14,-128550,-94.5105], Tmin=(949.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1018.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(Cds-CdsCsOs) + group(COCsFO) + group(CdCFF) + ring(2methyleneoxetane)"""),
)

species(
    label = 'O=C(F)C(=O)C(F)C=C(F)F(11812)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1050.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160912,0.0952203,-0.000133133,8.91165e-08,-1.77587e-11,-126217,31.0646], Tmin=(100,'K'), Tmax=(588.158,'K')), NASAPolynomial(coeffs=[11.2415,0.0351436,-1.88901e-05,3.79957e-09,-2.71174e-13,-127785,-18.7588], Tmin=(588.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1050.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(COCFO) + group(CdCFF)"""),
)

species(
    label = 'O=C(F)C(O)C(F)=C=C(F)F(11813)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u0 p0 c0 {8,D} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-943.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.124949,0.0911416,-0.000120763,8.25057e-08,-2.25548e-11,-113328,34.9766], Tmin=(100,'K'), Tmax=(890.717,'K')), NASAPolynomial(coeffs=[14.1913,0.0279727,-1.43837e-05,2.88476e-09,-2.07301e-13,-115834,-31.2574], Tmin=(890.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-943.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CdCddCF) + group(COCsFO) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = 'F[C]1OOC1C(F)[C]=C(F)F(11814)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {13,S}
9  C u1 p0 c0 {2,S} {6,S} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-330.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925337,0.0798614,-8.44697e-05,1.2294e-08,3.3354e-11,-39594.9,32.9317], Tmin=(100,'K'), Tmax=(502.831,'K')), NASAPolynomial(coeffs=[8.39412,0.039581,-2.13861e-05,4.33048e-09,-3.10815e-13,-40587.9,-0.371268], Tmin=(502.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-330.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C1[C](F)OC(=C(F)F)C1F(11624)',
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
    label = '[O]C1(F)OC1C(F)[C]=C(F)F(11815)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
9  C u0 p0 c0 {1,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-510.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749238,0.0798032,-9.55005e-05,6.46987e-08,-1.86507e-11,-61336,31.7726], Tmin=(100,'K'), Tmax=(823.876,'K')), NASAPolynomial(coeffs=[9.19903,0.0387779,-2.08061e-05,4.25626e-09,-3.09501e-13,-62728.3,-7.35559], Tmin=(823.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFOO) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C1C(F)C(=C(F)F)C1([O])F(11775)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-606.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296682,0.0898079,-0.000110194,7.15057e-08,-1.91125e-11,-72775.6,25.1136], Tmin=(100,'K'), Tmax=(896.76,'K')), NASAPolynomial(coeffs=[12.469,0.0355142,-1.93788e-05,3.99332e-09,-2.9162e-13,-74958.8,-32.2845], Tmin=(896.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C([C](O)F)C(F)[C]=C(F)F(11816)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {6,D} {7,S} {9,S}
9  C u1 p0 c0 {2,S} {5,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-730.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,280,501,1494,1531,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06953,'amu*angstrom^2'), symmetry=1, barrier=(47.5827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336787,'amu*angstrom^2'), symmetry=1, barrier=(7.74339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336821,'amu*angstrom^2'), symmetry=1, barrier=(7.74417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0679,'amu*angstrom^2'), symmetry=1, barrier=(47.5451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.488811,0.112833,-0.00020339,1.8898e-07,-6.71747e-11,-87736,34.2654], Tmin=(100,'K'), Tmax=(841.088,'K')), NASAPolynomial(coeffs=[9.11687,0.0387683,-2.06845e-05,4.04251e-09,-2.79561e-13,-88347.9,-4.44601], Tmin=(841.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCOF1sO2s) + radical(Cds_S)"""),
)

species(
    label = 'F[CH][C]=C(F)F(5078)',
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
    label = 'CFO(51)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC(F)[C]=C(F)F(10277)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {4,D} {5,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-441.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,562,600,623,1070,1265,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.803413,'amu*angstrom^2'), symmetry=1, barrier=(18.472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807709,'amu*angstrom^2'), symmetry=1, barrier=(18.5708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19794,0.0689547,-0.000113583,1.03539e-07,-3.73202e-11,-52965.7,24.3135], Tmin=(100,'K'), Tmax=(801.967,'K')), NASAPolynomial(coeffs=[7.04714,0.0279297,-1.46838e-05,2.89944e-09,-2.03491e-13,-53522.8,-0.238646], Tmin=(801.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(Cds_S)"""),
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
    label = 'O=C(F)C(=O)C(F)[C]=C(F)F(11817)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {5,D} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {6,D} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-812.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,286,619,818,1246,1924,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.635043,'amu*angstrom^2'), symmetry=1, barrier=(14.6009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416814,'amu*angstrom^2'), symmetry=1, barrier=(9.58338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63871,'amu*angstrom^2'), symmetry=1, barrier=(37.6772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.304506,0.106036,-0.000186424,1.69304e-07,-5.974e-11,-97595.2,33.1455], Tmin=(100,'K'), Tmax=(814.528,'K')), NASAPolynomial(coeffs=[10.6505,0.0335437,-1.84973e-05,3.68372e-09,-2.58412e-13,-98759.6,-13.652], Tmin=(814.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-812.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(COCFO) + group(CdCFF) + radical(Cds_S)"""),
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
    label = '[O]C(C=C=C(F)F)C(=O)F(11818)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {6,S} {10,D} {12,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {7,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-508.692,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,94,120,354,641,825,1294,540,610,2055,326.011,1191.39,4000],'cm^-1')),
        HinderedRotor(inertia=(0.145586,'amu*angstrom^2'), symmetry=1, barrier=(10.4851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341066,'amu*angstrom^2'), symmetry=1, barrier=(25.0119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803836,0.0747115,-9.032e-05,5.72588e-08,-1.4686e-11,-61070.1,32.9618], Tmin=(100,'K'), Tmax=(942.467,'K')), NASAPolynomial(coeffs=[12.2436,0.0261588,-1.30449e-05,2.59711e-09,-1.86374e-13,-63226.4,-21.5508], Tmin=(942.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-508.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCddFF) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(509)',
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
    label = '[O]C(C(=O)F)C(F)=C=C(F)F(11819)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u0 p0 c0 {8,D} {10,D}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-699.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,145,326,398,834,1303,486,617,768,1157,1926,94,120,354,641,825,1294,540,610,2055,435.097,441.225,4000],'cm^-1')),
        HinderedRotor(inertia=(0.212203,'amu*angstrom^2'), symmetry=1, barrier=(28.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00222937,'amu*angstrom^2'), symmetry=1, barrier=(7.88285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205305,0.090189,-0.000131865,1.00796e-07,-3.08084e-11,-84017.9,35.4632], Tmin=(100,'K'), Tmax=(800.03,'K')), NASAPolynomial(coeffs=[12.7401,0.0275141,-1.43482e-05,2.86305e-09,-2.03943e-13,-86023.4,-22.2126], Tmin=(800.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CdCddCF) + group(COCsFO) + group(CdCddFF) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C(=O)F)C(F)C#CF(11820)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  C u0 p0 c0 {7,S} {10,T}
10 C u0 p0 c0 {3,S} {9,T}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-437.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,233,378,609,1068,1270,1314,3037,486,617,768,1157,1926,2175,525,239,401,1367,232.304,872.272,872.451,874.262],'cm^-1')),
        HinderedRotor(inertia=(0.143525,'amu*angstrom^2'), symmetry=1, barrier=(5.4867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143723,'amu*angstrom^2'), symmetry=1, barrier=(5.48694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144322,'amu*angstrom^2'), symmetry=1, barrier=(5.48567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597462,0.0823618,-0.000116504,8.49732e-08,-2.30169e-11,-52466.1,32.7368], Tmin=(100,'K'), Tmax=(643.167,'K')), NASAPolynomial(coeffs=[10.7201,0.0290767,-1.47835e-05,2.91179e-09,-2.05593e-13,-53968.2,-13.186], Tmin=(643.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-437.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(COCsFO) + group(Ct-CtCs) + group(CtCF) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)C(=O)C=[C][C](F)F(11821)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,D} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {10,D} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u1 p0 c0 {2,S} {3,S} {10,S}
10 C u1 p0 c0 {7,D} {9,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-479.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,286,619,818,1246,1924,161,297,490,584,780,1358,1685,370,180,777.042],'cm^-1')),
        HinderedRotor(inertia=(0.216452,'amu*angstrom^2'), symmetry=1, barrier=(4.97665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214552,'amu*angstrom^2'), symmetry=1, barrier=(4.93296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93436,'amu*angstrom^2'), symmetry=1, barrier=(44.4747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340052,0.0900014,-0.000154451,1.40405e-07,-5.00402e-11,-57545.1,31.7034], Tmin=(100,'K'), Tmax=(803.065,'K')), NASAPolynomial(coeffs=[9.15247,0.0312892,-1.71072e-05,3.41177e-09,-2.40177e-13,-58482.7,-5.90364], Tmin=(803.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]C(=O)C(F)[C]=C(F)F(11756)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-389.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,562,600,623,1070,1265,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.619774,'amu*angstrom^2'), symmetry=1, barrier=(14.2498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.620692,'amu*angstrom^2'), symmetry=1, barrier=(14.2709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44328,'amu*angstrom^2'), symmetry=1, barrier=(33.184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0172708,0.0993881,-0.00018098,1.66395e-07,-5.85439e-11,-46677.1,31.803], Tmin=(100,'K'), Tmax=(837.658,'K')), NASAPolynomial(coeffs=[9.76529,0.0301996,-1.65431e-05,3.25972e-09,-2.26234e-13,-47515.9,-8.75767], Tmin=(837.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)[C](O)C(F)[C]=C(F)F(11822)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {6,D} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-700.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.132861,0.099713,-0.000158643,1.36999e-07,-4.78418e-11,-84077.1,38.3583], Tmin=(100,'K'), Tmax=(728.684,'K')), NASAPolynomial(coeffs=[11.1409,0.0348955,-1.91805e-05,3.88361e-09,-2.77715e-13,-85642.3,-11.9285], Tmin=(728.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-700.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C2CsJOH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)=C[C](F)F(11823)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-760.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.129531,0.0923307,-0.000128287,9.44991e-08,-2.81053e-11,-91317,36.6859], Tmin=(100,'K'), Tmax=(818.294,'K')), NASAPolynomial(coeffs=[12.6253,0.0312479,-1.63158e-05,3.27469e-09,-2.34696e-13,-93362,-21.0931], Tmin=(818.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-760.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(F)C(O)C(F)=[C][C](F)F(11824)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u1 p0 c0 {3,S} {4,S} {11,S}
11 C u1 p0 c0 {8,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-744.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0941877,0.0931692,-0.000125919,8.83713e-08,-2.49754e-11,-89462.2,37.3605], Tmin=(100,'K'), Tmax=(859.854,'K')), NASAPolynomial(coeffs=[13.6154,0.0302676,-1.61857e-05,3.29089e-09,-2.37998e-13,-91787.4,-25.8298], Tmin=(859.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-744.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(COCsFO) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsF1s))"""),
)

species(
    label = '[O]C(F)=C([O])C(F)C=C(F)F(11825)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {10,D}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {8,D}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-789.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.357708,0.106988,-0.000181308,1.64378e-07,-5.86061e-11,-94792.1,33.6438], Tmin=(100,'K'), Tmax=(801.042,'K')), NASAPolynomial(coeffs=[9.98212,0.0382891,-2.07076e-05,4.12104e-09,-2.9003e-13,-95901.1,-10.5283], Tmin=(801.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-789.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C(F)C([CH][C]=C(F)F)OF(11826)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u1 p0 c0 {7,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 C u0 p0 c0 {2,S} {3,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-434.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,486,617,768,1157,1926,562,600,623,1070,1265,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.035276,0.089961,-0.000101457,5.59819e-08,-1.22283e-11,-52099.9,36.0816], Tmin=(100,'K'), Tmax=(1109.08,'K')), NASAPolynomial(coeffs=[18.1841,0.0242509,-1.25853e-05,2.56111e-09,-1.86512e-13,-56141.2,-53.7024], Tmin=(1109.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-434.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C=CCJCO) + radical(Cdj(Cs-CsHH)(Cd-F1sF1s))"""),
)

species(
    label = '[O]C([CH]C(F)=C(F)F)C(=O)F(11827)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {12,S}
8  C u1 p0 c0 {7,S} {9,S} {13,S}
9  C u0 p0 c0 {1,S} {8,S} {11,D}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-748.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133959,0.0962534,-0.000125365,8.26711e-08,-2.16853e-11,-89851.9,34.5071], Tmin=(100,'K'), Tmax=(929.156,'K')), NASAPolynomial(coeffs=[15.8781,0.0273231,-1.40879e-05,2.83184e-09,-2.04006e-13,-92827.5,-41.5654], Tmin=(929.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-748.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(CdCsCdF) + group(COCsFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(C=OCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O=[C]C(OF)C(F)[C]=C(F)F(11828)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u1 p0 c0 {8,S} {9,D}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-328.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,562,600,623,1070,1265,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.16413,'amu*angstrom^2'), symmetry=1, barrier=(3.77367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116204,'amu*angstrom^2'), symmetry=1, barrier=(3.39537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160822,'amu*angstrom^2'), symmetry=1, barrier=(3.69761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0689979,'amu*angstrom^2'), symmetry=1, barrier=(20.0823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211448,0.0913029,-0.000121002,8.38519e-08,-2.3601e-11,-39366.8,37.7899], Tmin=(100,'K'), Tmax=(859.608,'K')), NASAPolynomial(coeffs=[12.9545,0.0320056,-1.75291e-05,3.60373e-09,-2.62375e-13,-41557.6,-21.7604], Tmin=(859.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)C(F)C(F)=C(F)F(11829)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 C u1 p0 c0 {6,D} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-641.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.235668,0.101051,-0.000154006,1.20318e-07,-3.66572e-11,-76954.4,37.5918], Tmin=(100,'K'), Tmax=(696.098,'K')), NASAPolynomial(coeffs=[13.6861,0.0287574,-1.48251e-05,2.92203e-09,-2.058e-13,-79079.2,-25.87], Tmin=(696.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-641.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)C(OF)C(F)[C]=[C]F(11830)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-300.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,486,617,768,1157,1926,1685,370,167,640,1190,336.828,339.463,1120.94],'cm^-1')),
        HinderedRotor(inertia=(0.084354,'amu*angstrom^2'), symmetry=1, barrier=(6.74674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0838714,'amu*angstrom^2'), symmetry=1, barrier=(6.73995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0870535,'amu*angstrom^2'), symmetry=1, barrier=(6.73993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330927,'amu*angstrom^2'), symmetry=1, barrier=(27.1134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.21254,0.101459,-0.0001545,1.24617e-07,-4.05081e-11,-36038.2,37.8044], Tmin=(100,'K'), Tmax=(751.09,'K')), NASAPolynomial(coeffs=[12.5941,0.0332557,-1.8292e-05,3.71888e-09,-2.67252e-13,-37962,-20.3148], Tmin=(751.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)C(F)=[C]F(11831)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {7,S} {11,D}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 C u1 p0 c0 {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-607.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,1380,1390,370,380,2900,435,246,474,533,1155,486,617,768,1157,1926,167,640,1190,385.75,388.891,390.199,390.28,1897.75],'cm^-1')),
        HinderedRotor(inertia=(0.0805526,'amu*angstrom^2'), symmetry=1, barrier=(8.69626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0809744,'amu*angstrom^2'), symmetry=1, barrier=(8.68952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00339405,'amu*angstrom^2'), symmetry=1, barrier=(8.6795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0729838,0.0982728,-0.000146712,1.13228e-07,-3.37792e-11,-72963,37.2583], Tmin=(100,'K'), Tmax=(661.871,'K')), NASAPolynomial(coeffs=[12.4842,0.0316301,-1.66343e-05,3.31467e-09,-2.35252e-13,-74827.8,-19.6706], Tmin=(661.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-607.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COCsFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=OCOJ) + radical(Cdj(Cd-CsF1s)(F1s))"""),
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
    E0 = (-207.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (189.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-113.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-68.1476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (123.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (244.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-199.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-118.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-129.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (95.5373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-156.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-85.4732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-127.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-33.6427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-169.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-111.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-144.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (49.1071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-130.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-62.0383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (101.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (12.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-5.95176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (36.3964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-53.7953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-25.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-132.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-65.7264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (131.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-76.9122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (200.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (31.6432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (195.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-13.3281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=CC(=O)F(335)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C(F)[C]=C(F)F(9746)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(373.072,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C(C(=O)F)C([CH]F)=C(F)F(9286)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C[C](F)OC(F)[C]=C(F)F(11545)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(139.556,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[CH]C(F)[C]=C(F)F(11811)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=C(F)F(3536)', '[O]C([CH]F)C(=O)F(1655)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C(F)C1OC(=C(F)F)C1F(11547)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C(F)C(=O)C(F)C=C(F)F(11812)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C(F)C(O)C(F)=C=C(F)F(11813)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['F[C]1OOC1C(F)[C]=C(F)F(11814)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C1[C](F)OC(=C(F)F)C1F(11624)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.67521e+10,'s^-1'), n=0.355, Ea=(50.9402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C1(F)OC1C(F)[C]=C(F)F(11815)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(122.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 122.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C1C(F)C(=C(F)F)C1([O])F(11775)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_cddouble]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C([C](O)F)C(F)[C]=C(F)F(11816)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(271.651,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CC(=O)F(335)', 'F[CH][C]=C(F)F(5078)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(88.7532,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CFO(51)', 'O=CC(F)[C]=C(F)F(10277)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(94.759,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'O=C(F)C(=O)C(F)[C]=C(F)F(11817)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(30.7997,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', '[O]C(C=C=C(F)F)C(=O)F(11818)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.4909,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][CH]C(=O)F(509)', 'FC=C=C(F)F(5101)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(22.9381,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(C(=O)F)C(F)=C=C(F)F(11819)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0.390446,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F(37)', '[O]C(C(=O)F)C(F)C#CF(11820)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][CH]C(=O)F(509)', 'F[CH][C]=C(F)F(5078)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(2.49665,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=C(F)C(=O)C=[C][C](F)F(11821)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(329.222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=[C]C(=O)C(F)[C]=C(F)F(11756)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(281.287,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C(F)[C](O)C(F)[C]=C(F)F(11822)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C(C(=O)F)C(F)=C[C](F)F(11823)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['O=C(F)C(O)C(F)=[C][C](F)F(11824)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C(F)=C([O])C(F)C=C(F)F(11825)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=C(F)C([CH][C]=C(F)F)OF(11826)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(140.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C([CH]C(F)=C(F)F)C(=O)F(11827)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(130.791,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=[C]C(OF)C(F)[C]=C(F)F(11828)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    products = ['[O]C([C]=O)C(F)C(F)=C(F)F(11829)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(239.347,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=C(F)C(OF)C(F)[C]=[C]F(11830)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(71.3721,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(C(=O)F)C(F)C(F)=[C]F(11831)'],
    products = ['[O]C(C(=O)F)C(F)[C]=C(F)F(11539)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2958',
    isomers = [
        '[O]C(C(=O)F)C(F)[C]=C(F)F(11539)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'FC=C=C(F)F(5101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2958',
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

