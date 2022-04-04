species(
    label = '[O]C(=O)C(F)[C]=C(F)F(9687)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-481.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,562,600,623,1070,1265,1685,370,180,180,180,180,316.803,1440.3,1440.61,1442.22],'cm^-1')),
        HinderedRotor(inertia=(2.23332,'amu*angstrom^2'), symmetry=1, barrier=(51.3485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.22919,'amu*angstrom^2'), symmetry=1, barrier=(51.2534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06299,0.0765455,-0.000139724,1.39641e-07,-5.34161e-11,-57774,26.2701], Tmin=(100,'K'), Tmax=(819.219,'K')), NASAPolynomial(coeffs=[3.79863,0.0373178,-2.05284e-05,4.09403e-09,-2.87714e-13,-57354.1,18.9161], Tmin=(819.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFF) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[O]C(=O)C([CH]F)=C(F)F(7600)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,D} {9,S}
7  C u1 p0 c0 {1,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-499.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1431.6,1431.88,1431.95,1432.03,1432.12],'cm^-1')),
        HinderedRotor(inertia=(0.28773,'amu*angstrom^2'), symmetry=1, barrier=(6.61548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287712,'amu*angstrom^2'), symmetry=1, barrier=(6.61507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00646,0.0806205,-0.000158354,1.62219e-07,-6.2119e-11,-59963.5,27.5456], Tmin=(100,'K'), Tmax=(835.736,'K')), NASAPolynomial(coeffs=[2.73442,0.0383312,-2.13938e-05,4.25984e-09,-2.97766e-13,-59064.3,26.627], Tmin=(835.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-499.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCFHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF) + radical(CCOJ) + radical(CsCdF1sH)"""),
)

species(
    label = 'O=C([CH][C]=C(F)F)OF(10513)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {7,D}
6  C u1 p0 c0 {7,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u0 p0 c0 {1,S} {2,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-262.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4263,0.063093,-6.08645e-05,2.80384e-08,-5.30952e-12,-31470.1,21.9207], Tmin=(100,'K'), Tmax=(1221.34,'K')), NASAPolynomial(coeffs=[12.7013,0.0261666,-1.55135e-05,3.28382e-09,-2.42486e-13,-34224.3,-34.7294], Tmin=(1221.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFF) + radical(C=CCJCO) + radical(Cds_S)"""),
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
    label = '[O]C([O])=CF(524)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {1,S} {4,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-216.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,194,682,905,1196,1383,3221,971.27,971.394],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65012,0.017516,-7.47389e-06,3.12837e-10,1.86785e-13,-26022.4,12.2475], Tmin=(100,'K'), Tmax=(2155.49,'K')), NASAPolynomial(coeffs=[14.5246,0.0047619,-3.76599e-06,7.6433e-10,-5.09564e-14,-32435.4,-52.569], Tmin=(2155.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'O=[C]C(F)[C]=C(F)F(10514)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 C u1 p0 c0 {4,D} {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-281.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,562,600,623,1070,1265,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.711412,'amu*angstrom^2'), symmetry=1, barrier=(16.3568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710501,'amu*angstrom^2'), symmetry=1, barrier=(16.3358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.101,0.0721779,-0.000129327,1.18221e-07,-4.13953e-11,-33724.4,25.4812], Tmin=(100,'K'), Tmax=(842.363,'K')), NASAPolynomial(coeffs=[8.07503,0.0226949,-1.2068e-05,2.35427e-09,-1.62555e-13,-34318.7,-3.52153], Tmin=(842.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1OC(=C(F)F)C1F(9691)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-825.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862286,0.0671481,-8.1126e-05,4.86324e-08,-1.09385e-11,-99116.4,24.7035], Tmin=(100,'K'), Tmax=(899.795,'K')), NASAPolynomial(coeffs=[14.312,0.0143726,-4.84027e-06,7.75399e-10,-4.90196e-14,-101821,-40.3412], Tmin=(899.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-825.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(CsCCFH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(CdCFF) + longDistanceInteraction_cyclic(Cs(F)-CO) + ring(4-Methylene-2-oxetanone)"""),
)

species(
    label = 'O=C(O)C(F)=C=C(F)F(7857)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {7,S} {9,D}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-667.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,180,180,1400.46,1400.84],'cm^-1')),
        HinderedRotor(inertia=(0.224529,'amu*angstrom^2'), symmetry=1, barrier=(5.16236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224591,'amu*angstrom^2'), symmetry=1, barrier=(5.16379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02111,0.0743013,-0.000129635,1.18695e-07,-4.17391e-11,-80199.8,27.196], Tmin=(100,'K'), Tmax=(846.527,'K')), NASAPolynomial(coeffs=[7.41982,0.0264133,-1.35e-05,2.60084e-09,-1.78607e-13,-80650.6,1.12803], Tmin=(846.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-667.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCddCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = 'FC(F)=[C]C(F)[C]1OO1(10515)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-90.5999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90806,0.0522831,-5.48621e-05,3.29031e-08,-8.69093e-12,-10826.6,28.012], Tmin=(100,'K'), Tmax=(876.576,'K')), NASAPolynomial(coeffs=[6.89854,0.0295105,-1.58934e-05,3.26604e-09,-2.3842e-13,-11701.5,4.5932], Tmin=(876.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.5999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(CdCFF) + ring(O2s-O2s-Cs(C-F)) + radical(Cs_P) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'O=C1O[C]([C](F)F)C1F(10497)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u1 p0 c0 {4,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-521.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34478,0.0598273,-7.05261e-05,4.49148e-08,-1.15127e-11,-62599.7,27.7402], Tmin=(100,'K'), Tmax=(948.825,'K')), NASAPolynomial(coeffs=[10.584,0.0208773,-8.95034e-06,1.65033e-09,-1.13229e-13,-64353,-16.3483], Tmin=(948.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-521.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C2CsJOC(O)) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = '[O]C1([O])C(=C(F)F)C1F(10516)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {3,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-317.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56848,0.0785257,-0.000115532,8.23073e-08,-2.20125e-11,-38042.3,26.0966], Tmin=(100,'K'), Tmax=(770.1,'K')), NASAPolynomial(coeffs=[14.5337,0.0141003,-5.84492e-06,1.03043e-09,-6.76772e-14,-40433.7,-39.1912], Tmin=(770.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cd(Cd-FF)-Cs) + radical(C=CC(C)(O)OJ) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'O=C(O)[C](F)[C]=C(F)F(10517)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {7,D}
6  C u1 p0 c0 {1,S} {7,S} {9,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-501.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768478,0.0796275,-0.00013628,1.2437e-07,-4.4514e-11,-60220.7,27.8502], Tmin=(100,'K'), Tmax=(804.157,'K')), NASAPolynomial(coeffs=[8.25733,0.0286549,-1.56044e-05,3.10721e-09,-2.18579e-13,-60981.4,-3.88826], Tmin=(804.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-501.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFF) + radical(Cs_P) + radical(Cds_S)"""),
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
    label = '[O]C(=O)C=C=C(F)F(10518)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {6,S} {8,D} {9,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u0 p0 c0 {1,S} {2,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-256.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,180,1820.34,1820.51,1820.58,1820.73],'cm^-1')),
        HinderedRotor(inertia=(0.360523,'amu*angstrom^2'), symmetry=1, barrier=(8.28914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60836,0.0654706,-0.000129348,1.34642e-07,-5.2046e-11,-30744.2,24.2325], Tmin=(100,'K'), Tmax=(839.924,'K')), NASAPolynomial(coeffs=[1.88578,0.0337919,-1.85597e-05,3.67567e-09,-2.56156e-13,-29719.9,29.3169], Tmin=(839.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(517)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.00186761,1.03198e-05,-1.52369e-08,5.8052e-12,3804.5,8.40409], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.36181,0.000422681,-4.06569e-07,1.52493e-10,-1.51968e-14,2816.74,-6.4394], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = '[O]C(=O)C(F)=C=C(F)F(10519)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {7,S} {9,D}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {9,D}
9 C u0 p0 c0 {6,D} {8,D}
"""),
    E0 = (-441.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,180,180,1609.37,1612.16,1613.88],'cm^-1')),
        HinderedRotor(inertia=(0.197712,'amu*angstrom^2'), symmetry=1, barrier=(4.54579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35803,0.0718748,-0.000143161,1.45564e-07,-5.47007e-11,-53070.6,26.1527], Tmin=(100,'K'), Tmax=(856.016,'K')), NASAPolynomial(coeffs=[2.77354,0.032751,-1.76381e-05,3.441e-09,-2.36801e-13,-52121.9,26.5009], Tmin=(856.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCddCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CCOJ)"""),
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
    label = '[O]C(=O)C(F)C#CF(10520)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {8,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u0 p0 c0 {5,S} {8,T}
8 C u0 p0 c0 {2,S} {7,T}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-259.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,239,401,1367,406.222,406.225,406.227,406.227,406.229,406.233,406.233],'cm^-1')),
        HinderedRotor(inertia=(0.00102155,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00102156,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59319,0.0633569,-0.000114538,1.152e-07,-4.40846e-11,-31191.6,21.973], Tmin=(100,'K'), Tmax=(831.121,'K')), NASAPolynomial(coeffs=[2.97538,0.0331976,-1.76811e-05,3.47736e-09,-2.42301e-13,-30609.5,20.4447], Tmin=(831.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(CtCF) + radical(CCOJ)"""),
)

species(
    label = '[O]C([O])=C(F)C=C(F)F(10521)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {9,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {4,S} {5,S} {7,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-550.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.934043,0.0736125,-0.000100673,7.1317e-08,-2.03774e-11,-66069.1,24.5249], Tmin=(100,'K'), Tmax=(849.861,'K')), NASAPolynomial(coeffs=[11.4637,0.0240558,-1.32104e-05,2.71123e-09,-1.96995e-13,-67858.9,-24.5621], Tmin=(849.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-550.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(CdCFF) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=O)[CH]C(F)=C(F)F(10522)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,D}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  C u0 p0 c0 {4,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-578.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53904,0.0631,-7.2345e-05,4.80399e-08,-1.43055e-11,-69454.5,21.4826], Tmin=(100,'K'), Tmax=(777.872,'K')), NASAPolynomial(coeffs=[6.6527,0.036805,-2.16405e-05,4.5852e-09,-3.39944e-13,-70250.1,-1.90345], Tmin=(777.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-578.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CdCsCdF) + group(Cds-OdCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O=C(OF)C(F)[C]=[C]F(10523)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u1 p0 c0 {6,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-121.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2850,1437.5,1250,1305,750,350,1685,370,167,640,1190,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561807,0.0850621,-0.000148733,1.36613e-07,-4.89442e-11,-14461.3,27.5588], Tmin=(100,'K'), Tmax=(807.999,'K')), NASAPolynomial(coeffs=[8.61683,0.0294571,-1.63065e-05,3.25783e-09,-2.29261e-13,-15249.5,-6.40757], Tmin=(807.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C(=O)C(F)C(F)=[C]F(10524)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-430.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,246,474,533,1155,167,640,1190,180,180,180,180,316.853,1376.54,1376.55,1376.71],'cm^-1')),
        HinderedRotor(inertia=(0.0391622,'amu*angstrom^2'), symmetry=1, barrier=(52.6629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.29129,'amu*angstrom^2'), symmetry=1, barrier=(52.6813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936959,0.079023,-0.000143383,1.40551e-07,-5.28108e-11,-51689.1,26.4475], Tmin=(100,'K'), Tmax=(823.479,'K')), NASAPolynomial(coeffs=[4.7821,0.0356759,-1.94875e-05,3.86956e-09,-2.71061e-13,-51485.9,13.7224], Tmin=(823.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-430.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCOJ) + radical(Cdj(Cd-CsF1s)(F1s))"""),
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
    E0 = (-217.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-109.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (211.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (242.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (225.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-209.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-139.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (172.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (24.2086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-53.4392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-6.95236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (114.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-52.0972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (33.0515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-201.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (111.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (94.0514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-35.3799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-80.5495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (206.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-9.59089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['CO2(14)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['[O]C(=O)C([CH]F)=C(F)F(7600)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.77375e+11,'s^-1'), n=0.335, Ea=(108.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HR!H)CJ;CJ;CO] + [cCs(-HR!H)CJ;CdsJ;C] for rate rule [cCs(-HR!H)CJ;CdsJ;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C([CH][C]=C(F)F)OF(10513)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(210.38,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=C(F)F(3536)', '[O]C([O])=CF(524)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=[C]C(F)[C]=C(F)F(10514)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['O=C1OC(=C(F)F)C1F(9691)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['O=C(O)C(F)=C=C(F)F(7857)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['FC(F)=[C]C(F)[C]1OO1(10515)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(390.543,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 388.1 to 390.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['O=C1O[C]([C](F)F)C1F(10497)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['[O]C1([O])C(=C(F)F)C1F(10516)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.241e+10,'s^-1'), n=0.754, Ea=(164.522,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cddouble] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cddouble]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C(O)[C](F)[C]=C(F)F(10517)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(231.468,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(=O)C=C=C(F)F(10518)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(34.6022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=O(517)', 'FC=C=C(F)F(5101)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(17.8875,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', '[O]C(=O)C(F)=C=C(F)F(10519)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(58.9142,'m^3/(mol*s)'), n=1.72001, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2675934934080712, var=0.0038470988293168424, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CO2(14)', 'F[CH][C]=C(F)F(5078)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(139.45,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', '[O]C(=O)C(F)C#CF(10520)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(35.3767,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]=O(517)', 'F[CH][C]=C(F)F(5078)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['[O]C([O])=C(F)C=C(F)F(10521)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    products = ['[O]C(=O)[CH]C(F)=C(F)F(10522)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(137.412,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C(OF)C(F)[C]=[C]F(10523)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(64.4673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(=O)C(F)C(F)=[C]F(10524)'],
    products = ['[O]C(=O)C(F)[C]=C(F)F(9687)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(157.818,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2575',
    isomers = [
        '[O]C(=O)C(F)[C]=C(F)F(9687)',
    ],
    reactants = [
        ('CO2(14)', 'FC=C=C(F)F(5101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2575',
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

