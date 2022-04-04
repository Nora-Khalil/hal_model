species(
    label = 'O=C(F)[CH]OC(F)[C]=C(F)F(11543)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,D} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-663.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,3025,407.5,1350,352.5,611,648,830,1210,1753,562,600,623,1070,1265,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.586653,0.107085,-0.000153775,1.0966e-07,-3.07235e-11,-79644.8,35.0507], Tmin=(100,'K'), Tmax=(876.343,'K')), NASAPolynomial(coeffs=[17.4213,0.0248872,-1.30783e-05,2.62448e-09,-1.88008e-13,-82801,-49.4499], Tmin=(876.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-663.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(CCsJOCs) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
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
    label = '[CH]C(=O)F(8932)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.05,1655.52],'cm^-1')),
        HinderedRotor(inertia=(0.0191536,'amu*angstrom^2'), symmetry=1, barrier=(5.31745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10758e-05,3.91582e-09,-5.6142e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.41,'K')), NASAPolynomial(coeffs=[6.53552,0.00717479,-3.36979e-06,6.65147e-10,-4.72043e-14,-1600.48,-4.84321], Tmin=(1580.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C(F)[C]=C(F)F(9673)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-270.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,562,600,623,1070,1265,1685,370,180,180,758.456],'cm^-1')),
        HinderedRotor(inertia=(0.277283,'amu*angstrom^2'), symmetry=1, barrier=(6.37527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06534,0.0458699,-5.49839e-05,3.47006e-08,-8.89665e-12,-32523.4,21.7549], Tmin=(100,'K'), Tmax=(939.968,'K')), NASAPolynomial(coeffs=[8.9134,0.0167286,-8.48099e-06,1.71921e-09,-1.24821e-13,-33810.8,-10.8591], Tmin=(939.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCdH)) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
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
    label = 'O=C(F)[CH]O[CH]F(1613)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {3,S} {6,S} {8,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 C u1 p0 c0 {2,S} {3,S} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-391.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,580,1155,1237,1373,3147,329.164,329.274,329.29,329.722],'cm^-1')),
        HinderedRotor(inertia=(0.0949268,'amu*angstrom^2'), symmetry=1, barrier=(7.31554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951287,'amu*angstrom^2'), symmetry=1, barrier=(7.31657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270803,'amu*angstrom^2'), symmetry=1, barrier=(20.8614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22301,0.0655239,-0.000101563,8.06592e-08,-2.52935e-11,-47003,23.5677], Tmin=(100,'K'), Tmax=(784.081,'K')), NASAPolynomial(coeffs=[10.7363,0.01699,-8.71139e-06,1.70927e-09,-1.19954e-13,-48494.8,-20.0142], Tmin=(784.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-391.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFHHO) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(F1s)(O2s-Cs)(H))"""),
)

species(
    label = 'O=C(F)C1OC(F)C1=C(F)F(11549)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-992.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273353,0.0695839,-5.01269e-05,1.13198e-08,8.73053e-13,-119283,30.2175], Tmin=(100,'K'), Tmax=(1152.17,'K')), NASAPolynomial(coeffs=[19.9518,0.0199022,-9.70868e-06,1.97139e-09,-1.44554e-13,-125055,-72.875], Tmin=(1152.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-992.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(COCsFO) + group(CdCFF) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(F)COC(F)=C=C(F)F(11672)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-900.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.137125,0.0933279,-0.00011347,6.78479e-08,-1.59872e-11,-108112,31.8971], Tmin=(100,'K'), Tmax=(1034.61,'K')), NASAPolynomial(coeffs=[17.9805,0.0232816,-1.19158e-05,2.40966e-09,-1.74887e-13,-111861,-56.1265], Tmin=(1034.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-900.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCddFO) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = 'F[C]1OC1OC(F)[C]=C(F)F(11673)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
9  C u1 p0 c0 {2,S} {6,S} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-513.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163835,0.0890961,-0.000104565,6.2156e-08,-1.48562e-11,-61619.1,31.8499], Tmin=(100,'K'), Tmax=(1009.48,'K')), NASAPolynomial(coeffs=[15.3618,0.0288766,-1.50868e-05,3.06541e-09,-2.22628e-13,-64687.6,-41.6154], Tmin=(1009.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-513.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C]1[CH]OC(F)C(=C(F)F)O1(11622)',
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
    label = '[O]C1(F)[CH]OC(F)C1=C(F)F(11674)',
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
    label = 'O=C(F)[CH]OC=C=C(F)F(11675)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {4,S} {8,S} {11,S}
7  C u0 p0 c0 {4,S} {10,D} {12,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {7,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-540.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,611,648,830,1210,1753,94,120,354,641,825,1294,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98286,'amu*angstrom^2'), symmetry=1, barrier=(22.5979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.982667,'amu*angstrom^2'), symmetry=1, barrier=(22.5935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.981952,'amu*angstrom^2'), symmetry=1, barrier=(22.577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.765728,0.0954104,-0.000113157,5.82725e-08,-1.02804e-11,-64797.8,31.9094], Tmin=(100,'K'), Tmax=(987.094,'K')), NASAPolynomial(coeffs=[25.8253,0.00472553,-1.29083e-06,2.39256e-10,-1.94992e-14,-70878.9,-100.244], Tmin=(987.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-540.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-CdsOsH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CCsJOC(O))"""),
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
    label = 'O=C(F)[CH]OC(F)=C=C(F)F(11676)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {5,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {11,D}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u0 p0 c0 {8,D} {10,D}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-706.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,197,221,431,657,611,648,830,1210,1753,94,120,354,641,825,1294,540,610,2055,449.018,449.043,449.057,449.107,449.121],'cm^-1')),
        HinderedRotor(inertia=(0.000835874,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125803,'amu*angstrom^2'), symmetry=1, barrier=(18.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125785,'amu*angstrom^2'), symmetry=1, barrier=(18.0035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.572327,0.103503,-0.000142986,9.38284e-08,-2.38238e-11,-84772.9,32.9342], Tmin=(100,'K'), Tmax=(970.554,'K')), NASAPolynomial(coeffs=[20.4472,0.0168753,-9.10343e-06,1.86676e-09,-1.36167e-13,-88853,-67.8447], Tmin=(970.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(CdCddFO) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=C=COC(F)[C]=C(F)F(11677)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
7  C u0 p0 c0 {4,S} {10,D} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 C u0 p0 c0 {5,D} {7,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-462.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1232,'amu*angstrom^2'), symmetry=1, barrier=(25.8247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12359,'amu*angstrom^2'), symmetry=1, barrier=(25.8334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12365,'amu*angstrom^2'), symmetry=1, barrier=(25.8348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.700565,0.0943111,-0.000115687,6.5826e-08,-1.41853e-11,-55491.9,31.2235], Tmin=(100,'K'), Tmax=(1155.33,'K')), NASAPolynomial(coeffs=[24.5244,0.0069784,-2.30251e-06,3.99968e-10,-2.81293e-14,-61320.6,-94.1148], Tmin=(1155.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-462.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)OsH) + group(CdCFF) + missing(Cdd-CdO2d) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = 'O=C(F)[CH]OC(F)C#CF(11678)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
7  C u1 p0 c0 {4,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  C u0 p0 c0 {6,S} {10,T}
10 C u0 p0 c0 {3,S} {9,T}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-470.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([242,464,653,889,1056,1342,3790,3025,407.5,1350,352.5,611,648,830,1210,1753,2175,525,239,401,1367,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.122914,0.0981328,-0.00015072,1.16204e-07,-3.44703e-11,-56405.2,30.5408], Tmin=(100,'K'), Tmax=(696.605,'K')), NASAPolynomial(coeffs=[14.0688,0.0254235,-1.30628e-05,2.55857e-09,-1.79244e-13,-58595.4,-34.3242], Tmin=(696.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-470.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + group(COCsFO) + group(Ct-CtCs) + group(CtCF) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(F)[CH]OC(F)=C[C](F)F(11679)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {8,D} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u1 p0 c0 {5,S} {11,S} {13,S}
10 C u1 p0 c0 {3,S} {4,S} {7,S}
11 C u0 p0 c0 {2,S} {6,D} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-768.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612493,0.106102,-0.000144146,9.51645e-08,-2.4575e-11,-92282.5,33.9899], Tmin=(100,'K'), Tmax=(949.988,'K')), NASAPolynomial(coeffs=[19.2568,0.0224403,-1.2046e-05,2.46146e-09,-1.78971e-13,-96057.6,-60.8483], Tmin=(949.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-768.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CsCFFH) + group(Cds-CdsCsH) + group(CdCFO) + group(COCsFO) + radical(CCsJOC(O)) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(F)COC(F)=[C][C](F)F(11680)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u1 p0 c0 {3,S} {4,S} {11,S}
11 C u1 p0 c0 {9,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-695.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459658,0.106333,-0.000161022,1.25447e-07,-3.89004e-11,-83447.2,35.1114], Tmin=(100,'K'), Tmax=(789.933,'K')), NASAPolynomial(coeffs=[14.5914,0.0301137,-1.62789e-05,3.28258e-09,-2.34778e-13,-85824.9,-33.951], Tmin=(789.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CsCFFH) + group(Cds-CdsCsH) + group(CdCFO) + group(COCsFO) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-F1sO2s))"""),
)

species(
    label = 'O=C(F)C(F)O[CH][C]=C(F)F(11681)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {6,D} {7,S}
9  C u1 p0 c0 {5,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-736.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.718666,0.0844928,-9.50237e-05,2.62449e-08,2.45444e-11,-88439.3,34.326], Tmin=(100,'K'), Tmax=(516.411,'K')), NASAPolynomial(coeffs=[9.04633,0.0395369,-2.12233e-05,4.27375e-09,-3.05506e-13,-89560.1,-2.87056], Tmin=(516.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C=CCJ(O)C) + radical(Cdj(Cs-O2sHH)(Cd-F1sF1s))"""),
)

species(
    label = 'O=C(F)[CH]O[CH]C(F)=C(F)F(11682)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {10,D}
8  C u1 p0 c0 {5,S} {7,S} {13,S}
9  C u1 p0 c0 {5,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {7,D}
11 C u0 p0 c0 {2,S} {6,D} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-754.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.647879,0.11237,-0.00018945,1.65357e-07,-5.64669e-11,-90641.5,35.8619], Tmin=(100,'K'), Tmax=(813.104,'K')), NASAPolynomial(coeffs=[12.6589,0.0330486,-1.75486e-05,3.45028e-09,-2.40645e-13,-92347.2,-22.7648], Tmin=(813.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-754.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(CdCsCdF) + group(COCsFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(C=CCJ(O)C) + radical(CCsJOCs)"""),
)

species(
    label = 'O=[C]C(F)OC(F)[C]=C(F)F(11600)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {11,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {6,D} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-660.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,355,410,600,1181,1341,1420,3056,562,600,623,1070,1265,1685,370,1855,455,950,214.377,214.377,214.378,214.378],'cm^-1')),
        HinderedRotor(inertia=(0.430411,'amu*angstrom^2'), symmetry=1, barrier=(14.0368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287492,'amu*angstrom^2'), symmetry=1, barrier=(9.37599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287495,'amu*angstrom^2'), symmetry=1, barrier=(9.37597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35241,'amu*angstrom^2'), symmetry=1, barrier=(44.1057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.839524,0.115465,-0.00019207,1.62586e-07,-5.38787e-11,-79226.6,35.3696], Tmin=(100,'K'), Tmax=(807.014,'K')), NASAPolynomial(coeffs=[14.7246,0.0294007,-1.55213e-05,3.04399e-09,-2.12044e-13,-81448.2,-34.5812], Tmin=(807.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-660.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O][C]=COC(F)C(F)=C(F)F(11683)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u1 p0 c0 {6,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-687.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16845,0.105874,-0.000131819,7.64086e-08,-1.68309e-11,-82501.1,35.7764], Tmin=(100,'K'), Tmax=(1128.08,'K')), NASAPolynomial(coeffs=[26.0094,0.00950577,-3.67874e-06,6.81267e-10,-4.85647e-14,-88632.9,-98.6163], Tmin=(1128.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)C(F)OC(F)[C]=[C]F(11684)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-624.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,162,485,641,721,926,1292,1380,3063,486,617,768,1157,1926,1685,370,167,640,1190,257.708,257.717,257.723,1416.25],'cm^-1')),
        HinderedRotor(inertia=(0.204549,'amu*angstrom^2'), symmetry=1, barrier=(9.64192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204565,'amu*angstrom^2'), symmetry=1, barrier=(9.64192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796486,'amu*angstrom^2'), symmetry=1, barrier=(37.545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796601,'amu*angstrom^2'), symmetry=1, barrier=(37.545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.109703,0.101064,-0.000148061,1.08934e-07,-2.88085e-11,-74956.9,34.5667], Tmin=(100,'K'), Tmax=(611.668,'K')), NASAPolynomial(coeffs=[11.9402,0.035435,-1.94201e-05,3.9327e-09,-2.81691e-13,-76677.4,-19.6582], Tmin=(611.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-624.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'O=C(F)[CH]OC(F)C(F)=[C]F(11685)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {7,S} {11,D}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {9,S}
11 C u1 p0 c0 {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-640.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([236,527,855,1015,1182,1348,3236,246,474,533,1155,3025,407.5,1350,352.5,611,648,830,1210,1753,167,640,1190,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.847021,0.114769,-0.00018395,1.49161e-07,-4.76139e-11,-76899.7,35.2491], Tmin=(100,'K'), Tmax=(771.098,'K')), NASAPolynomial(coeffs=[15.9139,0.0278246,-1.48196e-05,2.9381e-09,-2.06894e-13,-79484.6,-41.256], Tmin=(771.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-640.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cs-(Cds-O2d)OsHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COCsFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCsJOCs) + radical(Cdj(Cd-CsF1s)(F1s))"""),
)

species(
    label = '[O]C(F)(C=O)C(F)[C]=C(F)F(11541)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.191156,0.10176,-0.000165046,1.46684e-07,-5.24037e-11,-71735.4,35.4356], Tmin=(100,'K'), Tmax=(759.398,'K')), NASAPolynomial(coeffs=[10.2373,0.0377429,-2.06465e-05,4.15979e-09,-2.9613e-13,-73057.3,-10.2804], Tmin=(759.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
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
    label = 'F[C]=COC(F)[C]=C(F)F(11686)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {10,D} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 C u1 p0 c0 {4,S} {7,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-341.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,167,640,1190,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03047,'amu*angstrom^2'), symmetry=1, barrier=(23.6924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02988,'amu*angstrom^2'), symmetry=1, barrier=(23.6789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03029,'amu*angstrom^2'), symmetry=1, barrier=(23.6885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.307709,0.0958437,-0.000126076,7.97011e-08,-1.95264e-11,-40927.8,30.9833], Tmin=(100,'K'), Tmax=(1005.12,'K')), NASAPolynomial(coeffs=[19.4972,0.0170271,-8.45199e-06,1.6841e-09,-1.21321e-13,-44909,-64.6647], Tmin=(1005.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FC1=COC(F)C(=C(F)F)O1(11639)',
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
    label = 'OC(F)=COC(F)=C=C(F)F(11687)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-791.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331314,0.102841,-0.000150965,1.14387e-07,-3.45307e-11,-95036.8,33.6968], Tmin=(100,'K'), Tmax=(810.615,'K')), NASAPolynomial(coeffs=[14.4843,0.0297369,-1.56964e-05,3.14542e-09,-2.24591e-13,-97438.8,-34.6698], Tmin=(810.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-791.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCddFO) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O][C](F)C1OC(F)C1=C(F)F(11688)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {5,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {2,S} {6,S} {7,S}
11 C u0 p0 c0 {3,S} {4,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-622.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387115,0.080089,-8.19153e-05,4.13829e-08,-8.36597e-12,-74726.4,30.5451], Tmin=(100,'K'), Tmax=(1185.02,'K')), NASAPolynomial(coeffs=[16.5079,0.025673,-1.30346e-05,2.63166e-09,-1.90636e-13,-78547.1,-49.9649], Tmin=(1185.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cyclobutane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'O=C=[C]OC(F)[C]=C(F)F(11611)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u1 p0 c0 {6,S} {7,D}
9  C u1 p0 c0 {4,S} {10,D}
10 C u0 p0 c0 {5,D} {9,D}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-223.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,562,600,623,1070,1265,1670,1700,300,440,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.891591,'amu*angstrom^2'), symmetry=1, barrier=(20.4994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892044,'amu*angstrom^2'), symmetry=1, barrier=(20.5098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892392,'amu*angstrom^2'), symmetry=1, barrier=(20.5178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.112013,0.0915086,-0.000126229,8.29517e-08,-2.09937e-11,-26687.1,32.2473], Tmin=(100,'K'), Tmax=(976.704,'K')), NASAPolynomial(coeffs=[18.9098,0.0136074,-6.59151e-06,1.29209e-09,-9.20726e-14,-30402.9,-59.0739], Tmin=(976.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)=[C]OC(F)[C]=C(F)F(11689)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {6,S} {11,D}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {5,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-563.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,162,485,641,721,926,1292,1380,3063,293,496,537,1218,562,600,623,1070,1265,1670,1700,300,440,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.763615,0.10782,-0.00014718,9.68131e-08,-2.47271e-11,-67544.4,35.8038], Tmin=(100,'K'), Tmax=(963.958,'K')), NASAPolynomial(coeffs=[20.4599,0.0197517,-1.01376e-05,2.03572e-09,-1.46776e-13,-71636.1,-65.8085], Tmin=(963.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)=[C]OC(F)C=C(F)F(11690)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-687.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0647417,0.0965062,-0.000126367,8.46892e-08,-2.28028e-11,-82507,33.795], Tmin=(100,'K'), Tmax=(901.374,'K')), NASAPolynomial(coeffs=[14.8275,0.0304198,-1.63915e-05,3.35067e-09,-2.43394e-13,-85191.7,-36.5052], Tmin=(901.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)=CO[C](F)[C]=C(F)F(11691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  C u1 p0 c0 {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-597.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16553,0.109636,-0.000141275,8.46829e-08,-1.93792e-11,-71673,34.7957], Tmin=(100,'K'), Tmax=(1082.67,'K')), NASAPolynomial(coeffs=[25.3478,0.0116789,-5.55612e-06,1.11085e-09,-8.11744e-14,-77414,-95.2214], Tmin=(1082.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + radical(Cs_P) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = 'FOC(F)=CO[CH][C]=C(F)F(11692)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  C u1 p0 c0 {5,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-308.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00175824,0.0959561,-0.000131528,9.5754e-08,-2.83019e-11,-36968.1,35.5799], Tmin=(100,'K'), Tmax=(821.042,'K')), NASAPolynomial(coeffs=[12.7444,0.0338782,-1.81194e-05,3.67183e-09,-2.64704e-13,-39060.6,-23.3836], Tmin=(821.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF) + radical(C=CCJ(O)C) + radical(Cdj(Cs-O2sHH)(Cd-F1sF1s))"""),
)

species(
    label = 'FO[C]=COC(F)[C]=C(F)F(11693)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {11,D} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 C u1 p0 c0 {6,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-246.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,162,485,641,721,926,1292,1380,3063,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1670,1700,300,440,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07541,0.108208,-0.000141172,8.64528e-08,-2.02276e-11,-29495.4,37.3628], Tmin=(100,'K'), Tmax=(1059.52,'K')), NASAPolynomial(coeffs=[24.3212,0.0123274,-5.43075e-06,1.04186e-09,-7.43409e-14,-34877,-86.6296], Tmin=(1059.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(C=CJO)"""),
)

species(
    label = 'F[C]=[C]C(F)OC=C(F)OF(11694)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {7,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-196.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,162,485,641,721,926,1292,1380,3063,3010,987.5,1337.5,450,1655,326,540,652,719,1357,1685,370,167,640,1190,344.419,344.438,344.444,344.47],'cm^-1')),
        HinderedRotor(inertia=(0.415389,'amu*angstrom^2'), symmetry=1, barrier=(34.9733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0949659,'amu*angstrom^2'), symmetry=1, barrier=(7.98558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177377,'amu*angstrom^2'), symmetry=1, barrier=(14.9359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415308,'amu*angstrom^2'), symmetry=1, barrier=(34.9738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506924,0.107198,-0.000155115,1.13775e-07,-3.31807e-11,-23499.4,34.7547], Tmin=(100,'K'), Tmax=(838.387,'K')), NASAPolynomial(coeffs=[15.7715,0.0295359,-1.61717e-05,3.29568e-09,-2.37915e-13,-26229,-40.9104], Tmin=(838.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    E0 = (-238.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (149.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (229.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-230.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-145.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-54.9717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-176.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-144.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-185.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (17.7326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-154.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-69.2802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (75.0779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (69.2854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (9.96324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-55.8456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-43.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-46.5539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-92.5869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-48.6559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-52.2423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-19.6783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-44.8938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-41.9471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (326.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-230.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-220.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-163.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (211.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (12.4496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-137.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-77.7809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (183.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (230.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (270.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=CC(=O)F(335)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=O)F(8932)', '[O]C(F)[C]=C(F)F(9673)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C(F)F(3536)', 'O=C(F)[CH]O[CH]F(1613)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)C1OC(F)C1=C(F)F(11549)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)COC(F)=C=C(F)F(11672)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['F[C]1OC1OC(F)[C]=C(F)F(11673)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['F[C]1[CH]OC(F)C(=C(F)F)O1(11622)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.71916e+10,'s^-1'), n=0.520838, Ea=(62.1316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['[O]C1(F)[CH]OC(F)C1=C(F)F(11674)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.83275e+09,'s^-1'), n=0.690186, Ea=(94.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_cddouble] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra;radadd_intra_cddouble]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(=O)F(335)', 'F[CH][C]=C(F)F(5078)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(73.2924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O=C(F)[CH]OC=C=C(F)F(11675)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.9906,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]C(=O)F(509)', 'FC=C=C(F)F(5101)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C(F)[CH]OC(F)=C=C(F)F(11676)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=C=COC(F)[C]=C(F)F(11677)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.9385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C(F)[CH]OC(F)C#CF(11678)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(41.4438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(509)', 'F[CH][C]=C(F)F(5078)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)[CH]OC(F)=C[C](F)F(11679)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)COC(F)=[C][C](F)F(11680)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.08189e+06,'s^-1'), n=1.81713, Ea=(195.293,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_O;C_rad_out_H/OneDe;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)C(F)O[CH][C]=C(F)F(11681)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(191.873,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['O=C(F)[CH]O[CH]C(F)=C(F)F(11682)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(145.84,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C]C(F)OC(F)[C]=C(F)F(11600)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['[O][C]=COC(F)C(F)=C(F)F(11683)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(186.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C(F)C(F)OC(F)[C]=[C]F(11684)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(179.589,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(F)[CH]OC(F)C(F)=[C]F(11685)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(170.768,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(C=O)C(F)[C]=C(F)F(11541)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(130.567,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(6)', 'F[C]=COC(F)[C]=C(F)F(11686)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['FC1=COC(F)C(=C(F)F)O1(11639)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['OC(F)=COC(F)=C=C(F)F(11687)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['[O][C](F)C1OC(F)C1=C(F)F(11688)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(74.6464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', 'O=C=[C]OC(F)[C]=C(F)F(11611)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(290.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction30',
    reactants = ['OC(F)=[C]OC(F)[C]=C(F)F(11689)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    products = ['[O]C(F)=[C]OC(F)C=C(F)F(11690)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.8284271247461903
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OC(F)=CO[C](F)[C]=C(F)F(11691)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FOC(F)=CO[CH][C]=C(F)F(11692)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(66.4027,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction34',
    reactants = ['FO[C]=COC(F)[C]=C(F)F(11693)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(52.0458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F[C]=[C]C(F)OC=C(F)OF(11694)'],
    products = ['O=C(F)[CH]OC(F)[C]=C(F)F(11543)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(41.5946,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #2964',
    isomers = [
        'O=C(F)[CH]OC(F)[C]=C(F)F(11543)',
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

