species(
    label = '[O]C(C=[C]F)C(=O)F(5025)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-222.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,167,640,1190,269.43,270.135,3467.62],'cm^-1')),
        HinderedRotor(inertia=(0.140591,'amu*angstrom^2'), symmetry=1, barrier=(7.32104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0406701,'amu*angstrom^2'), symmetry=1, barrier=(19.1767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.551,0.0570591,-6.80761e-05,4.29776e-08,-1.1e-11,-26647.4,28.1078], Tmin=(100,'K'), Tmax=(943.885,'K')), NASAPolynomial(coeffs=[10.1701,0.0205332,-1.00306e-05,1.9804e-09,-1.41491e-13,-28274.5,-12.9768], Tmin=(943.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'C2HF(58)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (95.331,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0062,'amu')),
        LinearRotor(inertia=(51.6236,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([429.793,429.793,596.357,596.357,1107.96,2365.05,3506.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1870.76,'J/mol'), sigma=(4.25,'angstroms'), dipoleMoment=(1,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4498,0.0030263,3.99146e-05,-8.9615e-08,5.74336e-11,11468.6,5.90915], Tmin=(10,'K'), Tmax=(555.749,'K')), NASAPolynomial(coeffs=[4.23833,0.0086714,-5.87678e-06,1.96876e-09,-2.53031e-13,11206.2,0.995297], Tmin=(555.749,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.331,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C=[C]F(5065)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u1 p0 c0 {2,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-85.4923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.619121,'amu*angstrom^2'), symmetry=1, barrier=(14.2348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92895,0.050125,-7.24405e-05,5.46067e-08,-1.56392e-11,-10212,17.6888], Tmin=(100,'K'), Tmax=(649.489,'K')), NASAPolynomial(coeffs=[8.13552,0.0172054,-8.66369e-06,1.71849e-09,-1.22011e-13,-11130.1,-10.4371], Tmin=(649.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.4923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCdH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'O=C[C](F)OC=[C]F(5028)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u1 p0 c0 {1,S} {3,S} {7,S}
6  C u0 p0 c0 {3,S} {8,D} {9,S}
7  C u0 p0 c0 {4,D} {5,S} {10,S}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-247.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.51438,'amu*angstrom^2'), symmetry=1, barrier=(34.8186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5144,'amu*angstrom^2'), symmetry=1, barrier=(34.819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51508,'amu*angstrom^2'), symmetry=1, barrier=(34.8348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3525.94,'J/mol'), sigma=(5.52777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=550.74 K, Pc=47.37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09027,0.0656548,-7.83977e-05,4.71626e-08,-1.12597e-11,-29670.4,23.5101], Tmin=(100,'K'), Tmax=(1019.19,'K')), NASAPolynomial(coeffs=[13.1019,0.0185128,-9.01568e-06,1.77865e-09,-1.27257e-13,-32118.8,-34.6672], Tmin=(1019.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCOF1sO2s) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    label = 'O=C(F)C=C[C]F(5082)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,D} {6,S} {8,S}
5 C u0 p0 c0 {4,D} {7,S} {9,S}
6 C u0 p0 c0 {1,S} {3,D} {4,S}
7 C u2 p0 c0 {2,S} {5,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-161.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,255,533,799,832,1228,180,180,180,1353.49],'cm^-1')),
        HinderedRotor(inertia=(0.150135,'amu*angstrom^2'), symmetry=1, barrier=(3.4519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150543,'amu*angstrom^2'), symmetry=1, barrier=(3.46127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77635,0.0540438,-8.24987e-05,7.35472e-08,-2.68764e-11,-19289.4,20.8637], Tmin=(100,'K'), Tmax=(748.572,'K')), NASAPolynomial(coeffs=[6.40201,0.0243089,-1.28608e-05,2.57445e-09,-1.83206e-13,-19841.4,0.826162], Tmin=(748.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(COCFO) + radical(CsCCl_triplet)"""),
)

species(
    label = 'O=C(F)C1C=C(F)O1(5083)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,D}
8  C u0 p0 c0 {1,S} {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-548.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73752,0.0367757,2.77058e-06,-3.45511e-08,1.64698e-11,-65834.8,22.4879], Tmin=(100,'K'), Tmax=(997.979,'K')), NASAPolynomial(coeffs=[15.6393,0.0114349,-4.80194e-06,1.00964e-09,-7.93025e-14,-70122.4,-52.1316], Tmin=(997.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(CdCFO) + group(COCsFO) + ring(Cs-Cd-Cd(F)-O2s)"""),
)

species(
    label = 'O=C(F)C(=O)C=CF(5084)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,D} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  C u0 p0 c0 {2,S} {6,D} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-637.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60795,0.0575884,-7.13463e-05,4.95014e-08,-1.43188e-11,-76633.5,22.8382], Tmin=(100,'K'), Tmax=(829.633,'K')), NASAPolynomial(coeffs=[8.30369,0.0253066,-1.29818e-05,2.60311e-09,-1.87077e-13,-77744.5,-8.21451], Tmin=(829.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-637.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(COCFO) + group(CdCFH)"""),
)

species(
    label = 'F[C]=CC1OO[C]1F(5085)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {1,S} {4,S} {5,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (81.6185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7258,0.0516577,-5.17916e-05,2.68237e-08,-5.66862e-12,9896.97,25.9307], Tmin=(100,'K'), Tmax=(1124.95,'K')), NASAPolynomial(coeffs=[10.5804,0.0201727,-9.80899e-06,1.94366e-09,-1.39382e-13,7904.79,-17.8302], Tmin=(1124.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.6185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C1C=C(F)O[C]1F(5086)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u1 p0 c0 {1,S} {3,S} {5,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-287.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5657,0.0429859,-1.30583e-05,-2.20815e-08,1.39435e-11,-34460.5,19.8408], Tmin=(100,'K'), Tmax=(938.814,'K')), NASAPolynomial(coeffs=[15.5236,0.00994574,-2.49698e-06,4.06248e-10,-3.03141e-14,-38246,-52.8198], Tmin=(938.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFO) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)OC1C=[C]F(5087)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-99.2416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45658,0.0505933,-4.4527e-05,1.90975e-08,-3.23695e-12,-11840.1,25.2127], Tmin=(100,'K'), Tmax=(1415.2,'K')), NASAPolynomial(coeffs=[14.0789,0.0149163,-6.7118e-06,1.28348e-09,-8.99983e-14,-15412.7,-40.066], Tmin=(1415.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.2416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFOO) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C1C=C(F)C1([O])F(5088)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-181.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27921,0.0494074,-2.92612e-05,-2.06361e-09,4.88095e-12,-21780.8,22.0189], Tmin=(100,'K'), Tmax=(1063.92,'K')), NASAPolynomial(coeffs=[16.3951,0.0127774,-6.09846e-06,1.26911e-09,-9.58625e-14,-26140.5,-57.2162], Tmin=(1063.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd(F)-Cd-Cs-Cs) + radical(CC(C)OJ) + radical(CC(C)OJ) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'O=C(C=[C]F)[C](O)F(5089)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {5,D}
5  C u0 p0 c0 {4,D} {6,S} {7,S}
6  C u1 p0 c0 {1,S} {3,S} {5,S}
7  C u0 p0 c0 {5,S} {8,D} {9,S}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-302.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,280,501,1494,1531,3010,987.5,1337.5,450,1655,167,640,1190,180,2372.4],'cm^-1')),
        HinderedRotor(inertia=(0.441537,'amu*angstrom^2'), symmetry=1, barrier=(10.1518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85005,'amu*angstrom^2'), symmetry=1, barrier=(42.5364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8478,'amu*angstrom^2'), symmetry=1, barrier=(42.4847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07554,0.0737244,-0.000129138,1.21909e-07,-4.45142e-11,-36322.1,25.1986], Tmin=(100,'K'), Tmax=(824.284,'K')), NASAPolynomial(coeffs=[6.23475,0.0301408,-1.60745e-05,3.16673e-09,-2.20987e-13,-36542.5,5.12751], Tmin=(824.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(CdCFH) + radical(CsCOF1sO2s) + radical(Cdj(Cd-COH)(F1s))"""),
)

species(
    label = '[CH]=[C]F(489)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (350.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([167,640,1190,1142.58,1502.03,3807.5],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55592,0.00767983,-5.29841e-07,-4.54684e-09,2.16672e-12,42166.8,9.46361], Tmin=(100,'K'), Tmax=(1043.41,'K')), NASAPolynomial(coeffs=[5.86815,0.00351561,-1.29999e-06,2.62257e-10,-1.9892e-14,41428.4,-3.01582], Tmin=(1043.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
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
    label = 'O=CC=[C]F(564)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {7,S}
5 C u1 p0 c0 {1,S} {3,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-22.5172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(1.4261,'amu*angstrom^2'), symmetry=1, barrier=(32.7888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0457,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91323,0.00578592,9.15613e-05,-2.18108e-07,1.55212e-10,-2704.49,9.7028], Tmin=(10,'K'), Tmax=(471.623,'K')), NASAPolynomial(coeffs=[3.6152,0.0240295,-1.64471e-05,5.22424e-09,-6.26736e-13,-2851.16,9.06364], Tmin=(471.623,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-22.5172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCCD[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(=O)C=[C]F(5090)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {3,D} {6,S} {7,S}
6 C u0 p0 c0 {5,S} {8,D} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {5,S}
8 C u1 p0 c0 {2,S} {6,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-387.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,286,619,818,1246,1924,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.117461,'amu*angstrom^2'), symmetry=1, barrier=(2.70067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32383,'amu*angstrom^2'), symmetry=1, barrier=(30.4374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23249,0.0663245,-0.000105675,9.00907e-08,-3.07212e-11,-46533.2,24.3037], Tmin=(100,'K'), Tmax=(765.32,'K')), NASAPolynomial(coeffs=[9.14418,0.0216122,-1.14524e-05,2.2749e-09,-1.60596e-13,-47645.7,-11.1064], Tmin=(765.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(COCFO) + group(CdCFH) + radical(Cdj(Cd-COH)(F1s))"""),
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
    label = '[O]C(C#CF)C(=O)F(5091)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 C u0 p0 c0 {5,S} {8,T}
8 C u0 p0 c0 {2,S} {7,T}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-214.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,2175,525,239,401,1367,310.657,1062.51,1063.6],'cm^-1')),
        HinderedRotor(inertia=(0.0777803,'amu*angstrom^2'), symmetry=1, barrier=(5.58367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768593,'amu*angstrom^2'), symmetry=1, barrier=(5.5392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57856,0.0566946,-7.91799e-05,5.92591e-08,-1.7774e-11,-25703.8,26.5977], Tmin=(100,'K'), Tmax=(815.218,'K')), NASAPolynomial(coeffs=[9.37616,0.0184335,-8.77784e-06,1.68463e-09,-1.17426e-13,-26975.1,-9.42798], Tmin=(815.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(COCsFO) + group(Ct-CtCs) + group(CtCF) + radical(C=OCOJ)"""),
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
    label = 'O=[C]C(=O)C=[C]F(5076)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {5,S} {6,D} {8,S}
5 C u0 p0 c0 {2,D} {4,S} {7,S}
6 C u1 p0 c0 {1,S} {4,D}
7 C u1 p0 c0 {3,D} {5,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (35.7727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,167,640,1190,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.35735,'amu*angstrom^2'), symmetry=1, barrier=(8.21618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28529,'amu*angstrom^2'), symmetry=1, barrier=(29.5513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52784,0.0600164,-0.000101556,8.90842e-08,-3.04106e-11,4386.06,23.0542], Tmin=(100,'K'), Tmax=(828.521,'K')), NASAPolynomial(coeffs=[8.3427,0.0181167,-9.40705e-06,1.82872e-09,-1.26535e-13,3565.65,-6.67775], Tmin=(828.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.7727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(CdCFH) + group(Cds-O2d(Cds-O2d)H) + radical(Cdj(Cd-COH)(F1s)) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)C(O)=C[C]F(5092)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  C u2 p0 c0 {2,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-367.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742951,0.0732176,-9.94319e-05,6.53032e-08,-1.66276e-11,-44091.6,23.2493], Tmin=(100,'K'), Tmax=(967.983,'K')), NASAPolynomial(coeffs=[15.2605,0.0132262,-6.46738e-06,1.27629e-09,-9.12395e-14,-46902.2,-46.3168], Tmin=(967.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(COCFO) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C([C]=CF)C(=O)F(5093)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u0 p0 c0 {2,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-234.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49889,0.05923,-7.72596e-05,5.47755e-08,-1.58093e-11,-28123.1,28.2787], Tmin=(100,'K'), Tmax=(840.197,'K')), NASAPolynomial(coeffs=[9.30143,0.0220845,-1.09453e-05,2.15857e-09,-1.53512e-13,-29434.2,-8.00564], Tmin=(840.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = 'O=C(F)C(O)[C]=[C]F(5094)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u1 p0 c0 {5,S} {8,D}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-228.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24622,0.0649434,-9.03195e-05,6.57952e-08,-1.91782e-11,-27346.2,28.8952], Tmin=(100,'K'), Tmax=(837.725,'K')), NASAPolynomial(coeffs=[10.6202,0.0201831,-1.01715e-05,2.01153e-09,-1.4297e-13,-28916.8,-14.6689], Tmin=(837.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C(F)=C([O])C=CF(4926)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {1,S} {4,S} {5,D}
8  C u0 p0 c0 {2,S} {6,D} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-399.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180,632.754],'cm^-1')),
        HinderedRotor(inertia=(0.581043,'amu*angstrom^2'), symmetry=1, barrier=(13.3593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06386,0.0681805,-9.16553e-05,6.36643e-08,-1.75801e-11,-47970.2,23.4291], Tmin=(100,'K'), Tmax=(885.4,'K')), NASAPolynomial(coeffs=[11.7866,0.0197381,-9.58693e-06,1.87071e-09,-1.32236e-13,-49869,-26.9968], Tmin=(885.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]C(C=[C]F)OF(5095)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {1,S} {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u1 p0 c0 {4,D} {5,S}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (82.4579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1855,455,950,167,640,1190,439.039,439.047],'cm^-1')),
        HinderedRotor(inertia=(0.08588,'amu*angstrom^2'), symmetry=1, barrier=(11.7484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858916,'amu*angstrom^2'), symmetry=1, barrier=(11.7482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858798,'amu*angstrom^2'), symmetry=1, barrier=(11.7482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1414,0.064993,-8.2453e-05,5.23111e-08,-1.30847e-11,10018.6,28.8604], Tmin=(100,'K'), Tmax=(977.382,'K')), NASAPolynomial(coeffs=[12.9763,0.0165575,-8.11741e-06,1.6066e-09,-1.15067e-13,7705.21,-27.965], Tmin=(977.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.4579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(CCCJ=O) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C([C]=O)C=C(F)F(5096)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,D}
8  C u1 p0 c0 {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-259.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49054,0.0566636,-6.622e-05,3.99776e-08,-9.63807e-12,-31117.4,28.7302], Tmin=(100,'K'), Tmax=(1007.53,'K')), NASAPolynomial(coeffs=[11.3052,0.0176974,-8.20618e-06,1.58991e-09,-1.12674e-13,-33095,-18.6934], Tmin=(1007.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (-88.4512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (304.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (60.8356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (215.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-80.1669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-25.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (215.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-40.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (34.5774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-34.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (93.8631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (23.2748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (12.6679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-16.0175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (24.8217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (135.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (269.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (163.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (37.0089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (81.2958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-13.1706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (96.8386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (319.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (135.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['C2HF(58)', 'O=CC(=O)F(4234)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C=[C]F(5065)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(375.222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['O=C[C](F)OC=[C]F(5028)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(149.287,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'O=C(F)C=C[C]F(5082)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['O=C(F)C1C=C(F)O1(5083)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Cdsinglepri_rad_out]
Euclidian distance = 1.7320508075688772
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['O=C(F)C(=O)C=CF(5084)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['F[C]=CC1OO[C]1F(5085)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 303.2 to 303.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C1C=C(F)O[C]1F(5086)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.12592e+09,'s^-1'), n=0.421077, Ea=(48.3272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonyl_intra;radadd_intra_cdsingle]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C1(F)OC1C=[C]F(5087)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(123.029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 123.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C1C=C(F)C1([O])F(5088)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra;radadd_intra_cdsingle]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C(C=[C]F)[C](O)F(5089)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(262.843,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]F(489)', 'O=CC(=O)F(4234)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(22.1254,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CFO(51)', 'O=CC=[C]F(564)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(91.7248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=C(F)C(=O)C=[C]F(5090)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(26.0428,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2HF(58)', '[O][CH]C(=O)F(398)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(10.1487,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]C(C#CF)C(=O)F(5091)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2440,'m^3/(mol*s)'), n=1.64, Ea=(4.16086,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]F(489)', '[O][CH]C(=O)F(398)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'O=[C]C(=O)C=[C]F(5076)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(274.692,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['O=C(F)C(O)=C[C]F(5092)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C([C]=CF)C(=O)F(5093)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['O=C(F)C(O)[C]=[C]F(5094)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C(F)=C([O])C=CF(4926)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(323312,'s^-1'), n=2.11016, Ea=(185.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]C(C=[C]F)OF(5095)'],
    products = ['[O]C(C=[C]F)C(=O)F(5025)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(C=[C]F)C(=O)F(5025)'],
    products = ['[O]C([C]=O)C=C(F)F(5096)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(223.763,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #1504',
    isomers = [
        '[O]C(C=[C]F)C(=O)F(5025)',
    ],
    reactants = [
        ('C2HF(58)', 'O=CC(=O)F(4234)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1504',
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

