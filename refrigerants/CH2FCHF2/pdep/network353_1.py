species(
    label = '[O]C(O[C](F)C=O)C(=O)F(1908)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {10,S}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-720.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,280,501,1494,1531,486,617,768,1157,1926,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521126,0.0852815,-0.000137449,1.24736e-07,-4.52168e-11,-86585.7,33.7939], Tmin=(100,'K'), Tmax=(788.6,'K')), NASAPolynomial(coeffs=[7.69414,0.0359485,-1.89813e-05,3.76303e-09,-2.65161e-13,-87314.3,3.44472], Tmin=(788.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C(F)O[C](F)C=O(1896)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
7  C u1 p0 c0 {2,S} {3,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-564.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180,180,1019.25],'cm^-1')),
        HinderedRotor(inertia=(2.46014,'amu*angstrom^2'), symmetry=1, barrier=(56.5635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124129,'amu*angstrom^2'), symmetry=1, barrier=(2.85396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4576,'amu*angstrom^2'), symmetry=1, barrier=(56.505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3887.65,'J/mol'), sigma=(6.08201,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=607.24 K, Pc=39.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67532,0.0561309,-6.6877e-05,4.60525e-08,-1.34585e-11,-67766.1,24.0309], Tmin=(100,'K'), Tmax=(815.891,'K')), NASAPolynomial(coeffs=[7.55062,0.0273264,-1.39201e-05,2.78105e-09,-1.99493e-13,-68724.8,-3.1185], Tmin=(815.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sO2sH)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OO[C](F)C=O(1914)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u1 p0 c0 {2,S} {4,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-534.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,253,307,474,528,1484,1504,1502,1560,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,277.548],'cm^-1')),
        HinderedRotor(inertia=(0.941161,'amu*angstrom^2'), symmetry=1, barrier=(51.3828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940658,'amu*angstrom^2'), symmetry=1, barrier=(51.3846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94094,'amu*angstrom^2'), symmetry=1, barrier=(51.3849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939978,'amu*angstrom^2'), symmetry=1, barrier=(51.3831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941187,'amu*angstrom^2'), symmetry=1, barrier=(51.3827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3968.58,'J/mol'), sigma=(6.11625,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=619.88 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493724,0.0857462,-0.000132079,1.20146e-07,-4.47372e-11,-64193.4,30.3466], Tmin=(100,'K'), Tmax=(755.742,'K')), NASAPolynomial(coeffs=[6.98817,0.0407898,-2.18447e-05,4.37587e-09,-3.11206e-13,-64872.8,2.83287], Tmin=(755.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-534.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C]F(284)',
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
    label = '[O]C([O])C(=O)F(2004)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-304.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,180,1798.23,3800.51],'cm^-1')),
        HinderedRotor(inertia=(0.291432,'amu*angstrom^2'), symmetry=1, barrier=(6.70059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38172,0.0399157,-6.5594e-05,5.93607e-08,-2.10029e-11,-36549.3,21.4878], Tmin=(100,'K'), Tmax=(830.335,'K')), NASAPolynomial(coeffs=[5.67999,0.0160365,-8.02185e-06,1.54742e-09,-1.06924e-13,-36821.6,7.84743], Tmin=(830.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ) + radical(C=OCOJ)"""),
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
    label = 'O=C[C](F)O[CH]C(=O)F(2345)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {3,S} {8,S}
7  C u1 p0 c0 {3,S} {9,S} {10,S}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-569.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,611,648,830,1210,1753,243.643,243.68,243.691,695.249],'cm^-1')),
        HinderedRotor(inertia=(1.23607,'amu*angstrom^2'), symmetry=1, barrier=(52.0716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405083,'amu*angstrom^2'), symmetry=1, barrier=(17.0575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284079,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23643,'amu*angstrom^2'), symmetry=1, barrier=(52.072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721521,0.0779985,-0.00011018,8.28115e-08,-2.50977e-11,-68383.6,27.8968], Tmin=(100,'K'), Tmax=(804.27,'K')), NASAPolynomial(coeffs=[11.1393,0.0261821,-1.35331e-05,2.69356e-09,-1.9185e-13,-70059.2,-20.0931], Tmin=(804.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-569.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=CC1(F)OC(C(=O)F)O1(2346)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-990.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43563,0.044361,2.06039e-05,-1.5831e-07,1.61322e-10,-119068,25.8771], Tmin=(100,'K'), Tmax=(389.565,'K')), NASAPolynomial(coeffs=[3.97122,0.0420932,-2.26425e-05,4.64966e-09,-3.39186e-13,-119290,18.6017], Tmin=(389.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-990.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(COCsFO) + ring(O2s-Cs-O2s-Cs(F))"""),
)

species(
    label = 'O=CC(F)OC(=O)C(=O)F(2380)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1060.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970076,0.0740418,-8.47003e-05,5.21553e-08,-1.3471e-11,-127438,28.2163], Tmin=(100,'K'), Tmax=(917.634,'K')), NASAPolynomial(coeffs=[10.3298,0.0332418,-1.80064e-05,3.70129e-09,-2.70073e-13,-129156,-16.1345], Tmin=(917.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1060.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-OdCsH) + group(COCFO)"""),
)

species(
    label = 'O=C=C(F)OC(O)C(=O)F(1920)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-936.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287873,0.085065,-0.000108002,6.78984e-08,-1.68481e-11,-112476,31.9244], Tmin=(100,'K'), Tmax=(983.545,'K')), NASAPolynomial(coeffs=[15.9007,0.0215689,-1.11643e-05,2.25954e-09,-1.63773e-13,-115547,-43.1392], Tmin=(983.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-936.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-CsOsOsH) + group(COCsFO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = '[O]C(OC1(F)[CH]O1)C(=O)F(2381)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-574.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759245,0.0775665,-9.75052e-05,6.55105e-08,-1.80196e-11,-69022.9,32.1014], Tmin=(100,'K'), Tmax=(876.412,'K')), NASAPolynomial(coeffs=[11.2556,0.0296592,-1.55086e-05,3.13579e-09,-2.26496e-13,-70862.7,-17.1524], Tmin=(876.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-574.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-CsOsOsH) + group(COCsFO) + ring(Cs(F)(O2)-O2s-Cs) + radical(C=OCOJ) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring)"""),
)

species(
    label = 'O=C[C](F)OC1OO[C]1F(2382)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {1,S} {5,S} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-416.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731369,0.0793555,-0.000118962,1.05542e-07,-3.86963e-11,-50042.5,31.4979], Tmin=(100,'K'), Tmax=(738.947,'K')), NASAPolynomial(coeffs=[7.37457,0.0368839,-1.95309e-05,3.91279e-09,-2.78794e-13,-50846.6,2.66095], Tmin=(738.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C1OC(F)(C=O)O[C]1F(2325)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {4,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-692.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.339788,0.0793367,-0.000109904,7.59883e-08,-1.98574e-11,-83148.9,26.7234], Tmin=(100,'K'), Tmax=(1059.66,'K')), NASAPolynomial(coeffs=[15.3507,0.0132202,-2.9315e-06,2.6917e-10,-7.13024e-15,-85799.4,-44.061], Tmin=(1059.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-692.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFOO) + group(CsCFHO) + group(Cds-OdCsH) + ring(1,3-Dioxolane) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C(F)C1OO[CH][C](F)O1(2304)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u1 p0 c0 {1,S} {3,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-507.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.460808,0.079806,-8.61734e-05,4.59025e-08,-9.04763e-12,-60823,31.6232], Tmin=(100,'K'), Tmax=(1473.6,'K')), NASAPolynomial(coeffs=[19.1448,0.0108192,1.00415e-07,-3.89965e-10,3.79476e-14,-64889.1,-64.7545], Tmin=(1473.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(COCsFO) + ring(124trioxane) + radical(CsCsF1sO2s) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C1(F)OC1O[C](F)C=O(2383)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
9  C u1 p0 c0 {2,S} {4,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-597.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15609,0.0696947,-7.94867e-05,5.23929e-08,-1.48948e-11,-71809.3,28.3165], Tmin=(100,'K'), Tmax=(830.026,'K')), NASAPolynomial(coeffs=[8.0931,0.0362648,-1.9074e-05,3.87075e-09,-2.80353e-13,-72960.9,-3.85829], Tmin=(830.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFOO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C1OC(F)(C=O)C1([O])F(2372)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {5,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-598.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10137,0.0709799,-9.15025e-05,7.18319e-08,-2.43945e-11,-71875.6,30.2124], Tmin=(100,'K'), Tmax=(702.761,'K')), NASAPolynomial(coeffs=[6.97677,0.0375349,-2.01093e-05,4.09898e-09,-2.96887e-13,-72701.3,3.94002], Tmin=(702.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(O2s-Cs-Cs-Cs(F)) + radical(O2sj(Cs-F1sCsCs)) + radical(CCOJ) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C1OC(C(=O)F)O[C]1F(2384)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {5,S} {9,S} {12,S}
9  C u1 p0 c0 {1,S} {4,S} {8,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-712.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579544,0.0721599,-9.34436e-05,6.12539e-08,-1.51968e-11,-85568.6,29.0307], Tmin=(100,'K'), Tmax=(1115.85,'K')), NASAPolynomial(coeffs=[14.7469,0.0129598,-2.55174e-06,1.92599e-10,-2.32578e-15,-88206.5,-38.5241], Tmin=(1115.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-712.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(COCsFO) + ring(1,3-Dioxolane) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OC(=O)[C](O)F(2385)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {9,S} {12,S}
5  O u0 p2 c0 {7,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {5,D} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {10,S}
9  C u1 p0 c0 {1,S} {4,S} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-854.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,307,474,528,1484,1504,1502,1560,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.842938,0.0764853,-0.000109619,9.64312e-08,-3.54166e-11,-102627,30.025], Tmin=(100,'K'), Tmax=(746.354,'K')), NASAPolynomial(coeffs=[6.57878,0.038822,-2.00112e-05,3.96369e-09,-2.80795e-13,-103290,5.32262], Tmin=(746.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-854.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(CsCOF1sO2s)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CO[C](F)C=O(2108)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-563.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,421.486,421.686],'cm^-1')),
        HinderedRotor(inertia=(0.00094774,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401861,'amu*angstrom^2'), symmetry=1, barrier=(50.5001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407332,'amu*angstrom^2'), symmetry=1, barrier=(50.5155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.446,0.0344798,-2.01998e-05,5.03353e-09,-4.76771e-13,-67733.6,19.6715], Tmin=(100,'K'), Tmax=(2466.45,'K')), NASAPolynomial(coeffs=[17.94,0.00935137,-4.91706e-06,9.02554e-10,-5.80385e-14,-75376.3,-69.0649], Tmin=(2466.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CsCOF1sO2s)"""),
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
    label = 'O=C[C](F)OC(=O)C(=O)F(2386)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {1,S} {3,S} {10,S}
8  C u0 p0 c0 {3,S} {4,D} {9,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-921.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,286,619,818,1246,1924,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995617,0.0732643,-0.000101105,8.1961e-08,-2.81247e-11,-110695,28.3289], Tmin=(100,'K'), Tmax=(701.182,'K')), NASAPolynomial(coeffs=[7.72528,0.034872,-1.89705e-05,3.8659e-09,-2.79306e-13,-111638,-1.7486], Tmin=(701.182,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-921.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-OdCsH) + group(COCFO) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(OC(F)=C=O)C(=O)F(2387)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-692.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,197,221,431,657,2120,512.5,787.5,403.789,403.798,403.799,403.801,403.804,403.806],'cm^-1')),
        HinderedRotor(inertia=(0.00103385,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00103388,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00103388,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463717,0.0828674,-0.000114177,7.89361e-08,-2.15983e-11,-83169.2,32.0764], Tmin=(100,'K'), Tmax=(893.793,'K')), NASAPolynomial(coeffs=[14.1398,0.0216614,-1.14569e-05,2.31719e-09,-1.67113e-13,-85613.9,-32.3672], Tmin=(893.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-692.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-CsOsOsH) + group(COCsFO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(C=OCOJ)"""),
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
    label = 'O=[C]C(=O)O[C](F)C=O(2328)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {2,S} {8,S}
7  C u0 p0 c0 {2,S} {3,D} {9,S}
8  C u0 p0 c0 {4,D} {6,S} {10,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-499.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41989,0.0629373,-8.18573e-05,6.17164e-08,-1.9727e-11,-59951.3,26.5795], Tmin=(100,'K'), Tmax=(749.56,'K')), NASAPolynomial(coeffs=[7.57061,0.0301123,-1.61648e-05,3.28532e-09,-2.37411e-13,-60873.3,-1.32071], Tmin=(749.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-499.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(CsCOF1sO2s) + radical(OC=OCJ=O)"""),
)

species(
    label = 'O=C[C]F-2(1596)',
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
    label = 'O=C[C](F)O[C](O)C(=O)F(2388)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {3,S} {4,S} {9,S}
8  C u1 p0 c0 {1,S} {3,S} {10,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-759.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359542,0.0884287,-0.000142249,1.26378e-07,-4.50177e-11,-91208.3,34.2823], Tmin=(100,'K'), Tmax=(774.89,'K')), NASAPolynomial(coeffs=[9.06317,0.0339054,-1.81308e-05,3.61474e-09,-2.55666e-13,-92269.1,-3.62917], Tmin=(774.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-759.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(Cs_P) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(OC(F)[C]=O)C(=O)F(2389)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-701.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,355,410,600,1181,1341,1420,3056,486,617,768,1157,1926,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.086368,0.100501,-0.000176242,1.58124e-07,-5.46903e-11,-84242.4,35.4348], Tmin=(100,'K'), Tmax=(836.8,'K')), NASAPolynomial(coeffs=[10.6081,0.0302405,-1.59892e-05,3.11907e-09,-2.15667e-13,-85362.1,-10.2505], Tmin=(836.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-701.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O][C](F)C(=O)OC(F)C=O(2390)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {2,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-766.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54433,0.0875992,-0.000150055,1.42453e-07,-5.22908e-11,-92114.4,31.4368], Tmin=(100,'K'), Tmax=(832.519,'K')), NASAPolynomial(coeffs=[5.47643,0.0398139,-2.05567e-05,3.99759e-09,-2.77018e-13,-92100.9,13.5598], Tmin=(832.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-766.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sCOH)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C][C](F)OC(O)C(=O)F(2391)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {10,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-806.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0272646,0.0985677,-0.000173478,1.57827e-07,-5.52365e-11,-96816.4,34.4638], Tmin=(100,'K'), Tmax=(840.281,'K')), NASAPolynomial(coeffs=[9.43399,0.032347,-1.69908e-05,3.30428e-09,-2.28062e-13,-97640.2,-4.77698], Tmin=(840.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-806.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C]C(OF)O[C](F)C=O(2392)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u1 p0 c0 {1,S} {3,S} {9,S}
9  C u0 p0 c0 {5,D} {8,S} {12,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-421.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,280,501,1494,1531,2782.5,750,1395,475,1775,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997569,0.0786379,-8.43074e-05,2.8411e-09,4.48661e-11,-50593,31.3861], Tmin=(100,'K'), Tmax=(496.094,'K')), NASAPolynomial(coeffs=[9.07686,0.0352917,-1.91512e-05,3.84896e-09,-2.73868e-13,-51662.8,-4.63173], Tmin=(496.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-421.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C([C]=O)OC(F)(F)C=O(2393)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-674.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,267.341,267.356,267.386,2157.84],'cm^-1')),
        HinderedRotor(inertia=(0.136985,'amu*angstrom^2'), symmetry=1, barrier=(6.94984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13705,'amu*angstrom^2'), symmetry=1, barrier=(6.95011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235859,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564666,'amu*angstrom^2'), symmetry=1, barrier=(28.6471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19108,0.0923184,-0.000151942,1.34049e-07,-4.69535e-11,-80967.2,35.1248], Tmin=(100,'K'), Tmax=(785.318,'K')), NASAPolynomial(coeffs=[10.1758,0.0315502,-1.69403e-05,3.37328e-09,-2.37975e-13,-82229.8,-8.68683], Tmin=(785.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-674.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]C=O(1542)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,S} {5,S}
4 C u1 p0 c0 {2,D} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619047,'amu*angstrom^2'), symmetry=1, barrier=(0.260393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58699e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=C[C](F)OF(1556)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-259.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.060644,'amu*angstrom^2'), symmetry=1, barrier=(56.7934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.97821,'amu*angstrom^2'), symmetry=1, barrier=(68.4748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53489,0.035878,-5.03765e-05,4.63806e-08,-1.82263e-11,-31159.2,17.0234], Tmin=(100,'K'), Tmax=(714.201,'K')), NASAPolynomial(coeffs=[4.51401,0.0209835,-1.10922e-05,2.24131e-09,-1.61033e-13,-31344.7,8.8217], Tmin=(714.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(C(=O)F)C([O])C(=O)F(1905)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u1 p2 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-670.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,430,542,558,676,687,849,1103,1211,1906,1946,180,863.045,865.44,865.826,4000],'cm^-1')),
        HinderedRotor(inertia=(0.168896,'amu*angstrom^2'), symmetry=1, barrier=(3.88326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0268743,'amu*angstrom^2'), symmetry=1, barrier=(14.3119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6245,'amu*angstrom^2'), symmetry=1, barrier=(14.3585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.800564,0.0700854,-7.53059e-05,4.03237e-08,-8.57302e-12,-80489.2,37.5001], Tmin=(100,'K'), Tmax=(1137.9,'K')), NASAPolynomial(coeffs=[14.9392,0.0203835,-9.78669e-06,1.9369e-09,-1.39143e-13,-83706.8,-32.537], Tmin=(1137.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(COCsFO) + group(COCsFO) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=C(F)OC([O])C(=O)F(2394)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-399.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,293,496,537,1218,3120,650,792.5,1650,180,180,180,810.01,810.823],'cm^-1')),
        HinderedRotor(inertia=(0.00496985,'amu*angstrom^2'), symmetry=1, barrier=(2.31355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694072,'amu*angstrom^2'), symmetry=1, barrier=(15.9581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6936,'amu*angstrom^2'), symmetry=1, barrier=(15.9472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.581685,0.0790279,-0.000104448,6.87062e-08,-1.78348e-11,-47960.3,30.2638], Tmin=(100,'K'), Tmax=(941.826,'K')), NASAPolynomial(coeffs=[14.4973,0.0199284,-1.03244e-05,2.08254e-09,-1.5039e-13,-50581.5,-36.0373], Tmin=(941.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CdCFO) + group(COCsFO) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'O=C(F)C1OOC=C(F)O1(2309)',
    structure = adjacencyList("""1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 C u0 p0 c0 {1,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-704.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728013,0.0622415,-4.56912e-05,9.00091e-09,2.32384e-12,-84571.1,27.1471], Tmin=(100,'K'), Tmax=(1023.89,'K')), NASAPolynomial(coeffs=[16.9933,0.0171761,-6.74001e-06,1.26454e-09,-9.06468e-14,-88870.4,-56.4371], Tmin=(1023.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(CdCFO) + group(Cds-CdsOsH) + group(COCsFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
)

species(
    label = 'O=C(F)C(=O)OC(F)=CO(2395)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {9,S} {12,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {3,S} {5,D} {10,S}
9  C u0 p0 c0 {4,S} {7,D} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-1018.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0684076,0.09189,-0.000121627,7.95184e-08,-2.05704e-11,-122419,28.5889], Tmin=(100,'K'), Tmax=(942.946,'K')), NASAPolynomial(coeffs=[16.1784,0.0235524,-1.29208e-05,2.66384e-09,-1.94579e-13,-125457,-48.1863], Tmin=(942.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1018.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(COCFO)"""),
)

species(
    label = '[O]C1OC(F)=COO[C]1F(2295)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-332.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433154,0.0387376,8.88167e-05,-1.80366e-07,8.22083e-11,-39842.1,30.1855], Tmin=(100,'K'), Tmax=(917.045,'K')), NASAPolynomial(coeffs=[37.446,-0.0190746,1.38696e-05,-2.65197e-09,1.66862e-13,-50988.1,-168.934], Tmin=(917.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CdCFO) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCOJ) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O][C](F)C1OC(F)(C=O)O1(2361)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {5,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-619.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881633,0.0795063,-0.000134253,1.32025e-07,-5.09423e-11,-74440,31.9321], Tmin=(100,'K'), Tmax=(799.543,'K')), NASAPolynomial(coeffs=[3.99961,0.0419318,-2.25321e-05,4.49511e-09,-3.17469e-13,-74236.1,21.9798], Tmin=(799.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFOO) + group(CsCFHO) + group(Cds-OdCsH) + ring(O2s-Cs-O2s-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1OC(F)=COC1([O])F(2362)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {4,S} {6,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u0 p0 c0 {2,S} {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-580.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515345,0.074069,-8.38182e-05,4.95193e-08,-1.15206e-11,-69693.9,21.3062], Tmin=(100,'K'), Tmax=(1054.44,'K')), NASAPolynomial(coeffs=[14.6525,0.0204399,-7.52776e-06,1.28488e-09,-8.45496e-14,-72675.3,-47.6469], Tmin=(1054.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(CdCFO) + ring(23dihydro14dioxin) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCOJ) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C(O[C]=C=O)C(=O)F(2396)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {6,S}
8  C u1 p0 c0 {2,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-264.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,1685,370,2120,512.5,787.5,388.102,388.14,388.146,388.257,388.286],'cm^-1')),
        HinderedRotor(inertia=(0.115302,'amu*angstrom^2'), symmetry=1, barrier=(12.3307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111801,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115396,'amu*angstrom^2'), symmetry=1, barrier=(12.3296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675399,0.0760904,-0.000108604,7.59287e-08,-2.06507e-11,-31704.4,32.1931], Tmin=(100,'K'), Tmax=(906.335,'K')), NASAPolynomial(coeffs=[14.5392,0.0149016,-7.3311e-06,1.43308e-09,-1.01273e-13,-34217.3,-33.3278], Tmin=(906.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-CsOsOsH) + group(COCsFO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(OC(F)=[C]O)C(=O)F(2397)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {10,S} {12,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-604.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,293,496,537,1218,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00664869,0.0926,-0.00013023,9.06505e-08,-2.47495e-11,-72560.9,35.8114], Tmin=(100,'K'), Tmax=(899.39,'K')), NASAPolynomial(coeffs=[16.097,0.021037,-1.0874e-05,2.1763e-09,-1.55964e-13,-75455.1,-40.1084], Tmin=(899.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COCsFO) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)C(=O)O[C](F)[CH]O(2398)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {9,S} {12,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u0 p0 c0 {3,S} {5,D} {10,S}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-775.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.480889,0.112385,-0.000209118,1.92248e-07,-6.69927e-11,-93104.3,34.5117], Tmin=(100,'K'), Tmax=(852.654,'K')), NASAPolynomial(coeffs=[10.756,0.0314629,-1.71357e-05,3.3441e-09,-2.29817e-13,-93995.2,-11.896], Tmin=(852.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-775.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(O2s-H)(H))"""),
)

species(
    label = '[O]C=[C]OC(OF)C(=O)F(2399)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-413.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.592489,0.0908827,-0.000103107,5.07029e-08,-8.45838e-12,-49596.5,36.7315], Tmin=(100,'K'), Tmax=(1012.87,'K')), NASAPolynomial(coeffs=[25.0319,0.00623351,-2.24983e-06,4.47677e-10,-3.50695e-14,-55636.1,-91.4088], Tmin=(1012.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(COCsFO) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(O[C]=COF)C(=O)F(2400)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {4,S} {10,D} {12,S}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-288.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.217279,0.0919517,-0.000120606,7.5601e-08,-1.82435e-11,-34515.7,37.0553], Tmin=(100,'K'), Tmax=(1024.09,'K')), NASAPolynomial(coeffs=[19.752,0.0139527,-6.35831e-06,1.22671e-09,-8.71424e-14,-38605.7,-59.7599], Tmin=(1024.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(COCsFO) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C([C]=O)OC(F)=COF(2401)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {12,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-277.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,326,540,652,719,1357,3010,987.5,1337.5,450,1655,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.115858,0.0923964,-0.00013367,9.90948e-08,-2.92792e-11,-33270.9,34.5148], Tmin=(100,'K'), Tmax=(827.3,'K')), NASAPolynomial(coeffs=[13.734,0.0265527,-1.42869e-05,2.8921e-09,-2.07959e-13,-35524.2,-28.603], Tmin=(827.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(Cs-CsOsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
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
    E0 = (-253.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (151.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (81.2566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (186.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (141.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-244.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-164.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-228.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-69.5943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (50.8391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-199.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-39.3423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-130.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-130.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-148.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-91.8692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-179.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-171.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-200.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-12.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (38.8815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-17.6789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (199.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-99.1414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-96.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-112.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-159.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (151.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-31.1666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (314.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-64.8786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (311.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-236.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-228.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (135.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-131.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-112.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (205.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (26.1475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-164.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (132.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (247.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (235.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=CC(=O)F(335)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)O[C](F)C=O(1896)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(366.588,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C[C](F)OO[C](F)C=O(1914)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.25418e+20,'s^-1'), n=-1.9758, Ea=(148.136,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C]F(284)', '[O]C([O])C(=O)F(2004)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C[C](F)O[CH]C(=O)F(2345)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=CC1(F)OC(C(=O)F)O1(2346)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=CC(F)OC(=O)C(=O)F(2380)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C=C(F)OC(O)C(=O)F(1920)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C(OC1(F)[CH]O1)C(=O)F(2381)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C[C](F)OC1OO[C]1F(2382)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 303.2 to 303.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1OC(F)(C=O)O[C]1F(2325)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.47116e+08,'s^-1'), n=0.669085, Ea=(53.4924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C(F)C1OO[CH][C](F)O1(2304)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.334e+11,'s^-1'), n=0.280863, Ea=(213.707,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic
Ea raised from 212.4 to 213.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1(F)OC1O[C](F)C=O(2383)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(123.029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 123.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1OC(F)(C=O)C1([O])F(2372)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(122.465,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.0 to 122.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1OC(C(=O)F)O[C]1F(2384)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.16887e+10,'s^-1'), n=0.492453, Ea=(104.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_O] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[C](F)OC(=O)[C](O)F(2385)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(294.472,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=CC(=O)F(335)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(50.7403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(51)', 'O=CO[C](F)C=O(2108)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(114.517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C[C](F)OC(=O)C(=O)F(2386)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(40.9014,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(OC(F)=C=O)C(=O)F(2387)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][CH]C(=O)F(509)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'O=[C]C(=O)O[C](F)C=O(2328)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(294.788,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[C]F-2(1596)', '[O]C([O])C(=O)F(2004)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(257654,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C[C](F)O[C](O)C(=O)F(2388)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(OC(F)[C]=O)C(=O)F(2389)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_noH] for rate rule [R2H_S;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O][C](F)C(=O)OC(F)C=O(2390)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=[C][C](F)OC(O)C(=O)F(2391)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(36450.5,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_3;O_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=[C]C(OF)O[C](F)C=O(2392)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(105.471,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C([C]=O)OC(F)(F)C=O(2393)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(175.269,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]C=O(1542)', 'O=C[C](F)OF(1556)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.52781e-06,'m^3/(mol*s)'), n=3.38614, Ea=(181.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OR] for rate rule [Cdd_Cd;H_OR]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(C(=O)F)C([O])C(=O)F(1905)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(137.471,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(6)', '[CH]=C(F)OC([O])C(=O)F(2394)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C(F)C1OOC=C(F)O1(2309)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(16.6709,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination
Ea raised from 13.8 to 16.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C(F)C(=O)OC(F)=CO(2395)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1OC(F)=COO[C]1F(2295)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(388.244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra;radadd_intra] for rate rule [R7_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 383.3 to 388.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O][C](F)C1OC(F)(C=O)O1(2361)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['[O]C1OC(F)=COC1([O])F(2362)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(140.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_O] for rate rule [R7_SMSS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 139.1 to 140.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['HF(38)', '[O]C(O[C]=C=O)C(=O)F(2396)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(282.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C(OC(F)=[C]O)C(=O)F(2397)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    products = ['O=C(F)C(=O)O[C](F)[CH]O(2398)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]C=[C]OC(OF)C(=O)F(2399)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(78.4925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]C(O[C]=COF)C(=O)F(2400)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(68.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C([C]=O)OC(F)=COF(2401)'],
    products = ['[O]C(O[C](F)C=O)C(=O)F(1908)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(45.8647,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #353',
    isomers = [
        '[O]C(O[C](F)C=O)C(=O)F(1908)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #353',
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

