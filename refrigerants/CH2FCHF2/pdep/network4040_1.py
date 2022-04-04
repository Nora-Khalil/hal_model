species(
    label = '[O]C(=O)C1(F)[CH]C(F)=C1(14207)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-255.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,600,960.85,960.85,960.85,960.85,960.85,960.85,960.85,960.85,960.85,2336.45],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09683,0.0464608,-2.99293e-05,8.67302e-09,-1.00567e-12,-30689,24.0425], Tmin=(100,'K'), Tmax=(1876.27,'K')), NASAPolynomial(coeffs=[12.5672,0.0241392,-1.20841e-05,2.33236e-09,-1.60825e-13,-34618.1,-33.0597], Tmin=(1876.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(CCJCC=O)"""),
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
    label = 'FC1=CC(F)=C1(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (29.0416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,180,1170.89,2106.02,2106.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89011,0.00686793,0.000107406,-2.40335e-07,1.58082e-10,3500.24,9.34061], Tmin=(10,'K'), Tmax=(520.912,'K')), NASAPolynomial(coeffs=[4.12106,0.0279909,-1.93509e-05,6.26873e-09,-7.66374e-13,3165.53,5.39524], Tmin=(520.912,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(29.0416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCC(F)DC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(=O)[C](F)C1C=C1F(14277)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {6,D}
8  C u1 p0 c0 {2,S} {5,S} {9,S}
9  C u0 p0 c0 {3,S} {4,D} {8,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-113.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,234,347,1316,1464,180,1080.17,1080.17,1080.17,1080.17,1080.17,1080.17,1080.17,1080.17,1080.17,1080.17,2276.4],'cm^-1')),
        HinderedRotor(inertia=(0.0197189,'amu*angstrom^2'), symmetry=1, barrier=(0.453377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197189,'amu*angstrom^2'), symmetry=1, barrier=(0.453377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97232,0.0513121,-4.31888e-05,2.12579e-08,-4.93394e-12,-13638.1,25.2136], Tmin=(100,'K'), Tmax=(939.103,'K')), NASAPolynomial(coeffs=[5.59941,0.035863,-1.85127e-05,3.74054e-09,-2.70661e-13,-14319.4,7.94288], Tmin=(939.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cd-Cd-Cs(C-F)) + radical(CCOJ) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O]C(=O)C1[C](F)C=C1F(14208)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-287.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,228.383,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,1079.01,2256.44],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4369.4,'J/mol'), sigma=(6.16145,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.49 K, Pc=42.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01883,0.0509553,-4.0285e-05,1.75587e-08,-3.55112e-12,-34563.2,22.0219], Tmin=(100,'K'), Tmax=(1054.46,'K')), NASAPolynomial(coeffs=[6.07771,0.0355584,-1.83825e-05,3.71126e-09,-2.68088e-13,-35419.2,2.22484], Tmin=(1054.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C(=O)C1=CC(F)(F)[CH]1(14278)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {7,S} {8,D} {9,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {5,S} {6,D} {10,S}
9  C u0 p0 c0 {3,S} {4,D} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-169.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73075,0.0555036,-7.12168e-05,6.44954e-08,-2.61385e-11,-20279.1,26.1823], Tmin=(100,'K'), Tmax=(677.566,'K')), NASAPolynomial(coeffs=[4.2259,0.0364132,-1.93015e-05,3.91777e-09,-2.8295e-13,-20517.1,15.8545], Tmin=(677.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(CCOJ) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C(OF)[C]1[CH]C(F)=C1(14279)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u1 p0 c0 {6,S} {7,S} {9,S}
6  C u1 p0 c0 {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (40.2447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,3150,900,1100,271,519,563,612,1379,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35125,0.0501707,-3.26815e-05,8.25966e-09,-4.15081e-13,4942.31,25.7683], Tmin=(100,'K'), Tmax=(1371.21,'K')), NASAPolynomial(coeffs=[14.8463,0.0199103,-9.54089e-06,1.85237e-09,-1.29958e-13,385.279,-46.7197], Tmin=(1371.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.2447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cd(F)-Cd-Cs-Cs) + radical(CCJ(C)CO) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1O[CH]C(F)=C[C]1F(14280)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {8,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,D} {10,S}
6  C u1 p0 c0 {2,S} {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,D} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-407.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24157,0.0301831,1.88725e-05,-3.39193e-08,1.09616e-11,-48920.7,20.4105], Tmin=(100,'K'), Tmax=(1233.5,'K')), NASAPolynomial(coeffs=[9.88991,0.0326291,-1.72371e-05,3.50515e-09,-2.53039e-13,-52880.5,-26.496], Tmin=(1233.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cyclohexane) + radical(Cs_P) + radical(C=CCJ(O)C)"""),
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
    label = 'O=[C]C1(F)[CH]C(F)=C1(14281)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u1 p0 c0 {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u1 p0 c0 {3,D} {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-57.0185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,1855,455,950,180,396.704,1209.82,1209.83,1209.83],'cm^-1')),
        HinderedRotor(inertia=(0.0207939,'amu*angstrom^2'), symmetry=1, barrier=(21.597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60074,0.0453994,-3.25916e-05,8.58589e-09,-3.13429e-14,-6765.2,25.5321], Tmin=(100,'K'), Tmax=(1159.83,'K')), NASAPolynomial(coeffs=[13.2825,0.0152814,-6.79278e-06,1.31682e-09,-9.4041e-14,-10159,-35.5066], Tmin=(1159.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.0185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=C1OC2C(F)=CC12F(14212)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-476.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77077,0.0481377,-3.16836e-05,9.20888e-09,-1.04382e-12,-57273,19.102], Tmin=(100,'K'), Tmax=(2013.52,'K')), NASAPolynomial(coeffs=[16.9342,0.018015,-9.2437e-06,1.77926e-09,-1.21373e-13,-63379.5,-64.6657], Tmin=(2013.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-476.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + longDistanceInteraction_cyclic(Cs(F)-CO) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = '[O]C(=O)C1(F)C2[C](F)C21(14282)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,S} {4,D} {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-208.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60382,0.0589006,-7.63989e-05,7.12448e-08,-2.99663e-11,-24986.3,22.5966], Tmin=(100,'K'), Tmax=(663.681,'K')), NASAPolynomial(coeffs=[3.8717,0.039984,-2.17839e-05,4.46947e-09,-3.24785e-13,-25171.7,13.4559], Tmin=(663.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(Cds-OdCsOs) + polycyclic(s2_3_3_ane) + radical(CCOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FC1=CC(F)([C]2OO2)[CH]1(14283)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u1 p0 c0 {5,S} {9,S} {11,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u0 p0 c0 {5,S} {9,D} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (85.7416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38098,0.0490241,-2.93925e-05,3.11192e-09,1.91366e-12,10413.8,25.046], Tmin=(100,'K'), Tmax=(1139.03,'K')), NASAPolynomial(coeffs=[13.7611,0.0194229,-8.6824e-06,1.6849e-09,-1.2041e-13,6693.49,-40.2434], Tmin=(1139.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.7416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(O2s-O2s-Cs(C-F)) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Cs_P)"""),
)

species(
    label = 'O=C1OC2[C](F)[CH]C12F(14223)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
7  C u1 p0 c0 {5,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-204.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33919,0.0433915,-2.51018e-05,5.97511e-09,-5.32734e-13,-24566.6,21.1137], Tmin=(100,'K'), Tmax=(2640.97,'K')), NASAPolynomial(coeffs=[24.7845,0.00939529,-5.79249e-06,1.10072e-09,-7.13053e-14,-36421.9,-108.969], Tmin=(2640.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-OdCsOs) + polycyclic(s2_4_4_ane) + radical(CCJCC=O) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = 'O=C1OC2(F)[CH]C1(F)[CH]2(14238)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-290.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90106,0.0328393,1.86116e-05,-4.40726e-08,1.72465e-11,-34874.5,21.2819], Tmin=(100,'K'), Tmax=(1065.62,'K')), NASAPolynomial(coeffs=[13.228,0.0220227,-1.07861e-05,2.23627e-09,-1.67373e-13,-39088.4,-42.5293], Tmin=(1065.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsOs) + polycyclic(s3_4_5_ane) + radical(CCJCO) + radical(CCJCO) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = '[CH]=C(F)C=C(F)C([O])=O(14284)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,D} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {5,D} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {6,S}
9  C u1 p0 c0 {7,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-105.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,288,410,724,839,1320,250,446,589,854,899,3120,650,792.5,1650,180,180,180,1559.66,1559.8,1559.89],'cm^-1')),
        HinderedRotor(inertia=(0.198295,'amu*angstrom^2'), symmetry=1, barrier=(4.5592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198101,'amu*angstrom^2'), symmetry=1, barrier=(4.55473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882818,0.0773977,-0.000128816,1.19606e-07,-4.37807e-11,-12607.5,27.4587], Tmin=(100,'K'), Tmax=(804.022,'K')), NASAPolynomial(coeffs=[6.66125,0.0332083,-1.7567e-05,3.47583e-09,-2.44154e-13,-13037.6,3.94539], Tmin=(804.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(CdCCF) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = '[O]C1([O])C2C(F)=CC21F(14266)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {7,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-41.7487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74751,0.0587248,-5.18333e-05,-2.62932e-09,2.90273e-11,-4949.29,19.4412], Tmin=(100,'K'), Tmax=(488.046,'K')), NASAPolynomial(coeffs=[5.96326,0.0362191,-1.9687e-05,4.03422e-09,-2.92983e-13,-5504.25,0.656945], Tmin=(488.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.7487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCCF) + group(Cs-CsCsOsOs) + group(CdCsCdF) + group(Cds-CdsCsH) + polycyclic(s2_3_4_ene_1) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = '[O]C(=O)C1=CC(F)=C1(14285)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u1 p2 c0 {8,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,D} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {3,D} {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (154.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54387,0.0591795,-8.19647e-05,6.56844e-08,-2.19142e-11,18622,21.8794], Tmin=(100,'K'), Tmax=(724.974,'K')), NASAPolynomial(coeffs=[7.5461,0.0260647,-1.34537e-05,2.68786e-09,-1.91977e-13,17751.6,-5.14788], Tmin=(724.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cd-Cd-Cd(F)) + radical(CCOJ)"""),
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
    label = '[O]C(=O)C1(F)C=C=C1(14286)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,D} {9,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {2,S} {3,D} {4,S}
8  C u0 p0 c0 {5,D} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (121.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80731,0.0595849,-0.000105783,1.15636e-07,-4.8187e-11,14673.8,20.2473], Tmin=(100,'K'), Tmax=(801.357,'K')), NASAPolynomial(coeffs=[-0.318833,0.0429194,-2.35285e-05,4.72866e-09,-3.35245e-13,15890.5,35.4988], Tmin=(801.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(CCOJ)"""),
)

species(
    label = 'F[C]1[CH]C(F)=C1(615)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (86.8055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,271,519,563,612,1379,180,180,1104.81,1106.27,1801.34],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06134,0.0383117,-3.38918e-05,1.45559e-08,-2.45109e-12,10513.6,15.9062], Tmin=(100,'K'), Tmax=(1431.45,'K')), NASAPolynomial(coeffs=[12.0907,0.0102863,-4.52462e-06,8.78888e-10,-6.24417e-14,7642.3,-36.0769], Tmin=(1431.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdH)_ring) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring)"""),
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
    label = '[O]C(=O)C1=CC(F)=[C]1(14287)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {7,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u0 p0 c0 {4,D} {6,S} {9,S}
6 C u0 p0 c0 {1,S} {5,S} {8,D}
7 C u0 p0 c0 {2,S} {3,D} {4,S}
8 C u1 p0 c0 {4,S} {6,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (373.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,250,446,589,854,899,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36542,0.0656676,-0.000114205,1.05883e-07,-3.812e-11,45017.3,23.6755], Tmin=(100,'K'), Tmax=(819.861,'K')), NASAPolynomial(coeffs=[6.85773,0.0246632,-1.31901e-05,2.60518e-09,-1.8214e-13,44594.2,1.18129], Tmin=(819.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cd-Cd-Cd(F)) + radical(CCOJ) + radical(cyclobutadiene-C1)"""),
)

species(
    label = '[O]C(=O)C1(F)[C]=C=C1(14288)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {4,S} {8,D} {9,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 C u1 p0 c0 {4,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (359.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.2927,0.0263142,-1.05556e-05,1.99269e-10,2.30632e-13,43108.9,8.60383], Tmin=(100,'K'), Tmax=(2718.04,'K')), NASAPolynomial(coeffs=[53.005,-0.0225958,4.67833e-06,-6.5328e-10,4.3788e-14,9302.1,-283.766], Tmin=(2718.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C(=O)C1(F)[C]=C(F)C1(14289)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {5,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-203.105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,246,474,533,1155,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53422,0.0600394,-7.66487e-05,6.68891e-08,-2.60279e-11,-24344.5,24.0626], Tmin=(100,'K'), Tmax=(676.78,'K')), NASAPolynomial(coeffs=[4.83916,0.0374095,-1.96289e-05,3.9608e-09,-2.85051e-13,-24720.9,9.93243], Tmin=(676.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'O=C(O)C1(F)[C]C(F)=C1(14290)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {11,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  C u0 p0 c0 {2,S} {6,D} {9,S}
9  C u2 p0 c0 {5,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-233.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2950,1000,323,467,575,827,1418,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21895,0.0604183,-5.54925e-05,2.52534e-08,-4.63527e-12,-27932.5,24.4569], Tmin=(100,'K'), Tmax=(1289.42,'K')), NASAPolynomial(coeffs=[13.611,0.0219762,-1.07725e-05,2.13199e-09,-1.52368e-13,-31128.2,-38.4776], Tmin=(1289.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C(=O)[C]1C=C(F)C1F(14291)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-258.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29009,0.0432561,-2.42319e-05,5.59689e-09,-4.81741e-13,-30993.3,23.1163], Tmin=(100,'K'), Tmax=(2556.9,'K')), NASAPolynomial(coeffs=[26.066,0.00828215,-5.01744e-06,9.26792e-10,-5.83394e-14,-43877.9,-115.329], Tmin=(2556.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(CCJ(C)CO) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'O=C(OF)C1(F)[CH][C]=C1(14292)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u1 p0 c0 {5,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {10,S}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (114.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42792,0.0551987,-4.49174e-05,1.76585e-08,-2.80277e-12,13816.2,26.1511], Tmin=(100,'K'), Tmax=(1464.05,'K')), NASAPolynomial(coeffs=[13.6134,0.0219066,-1.0808e-05,2.12667e-09,-1.50578e-13,10248.1,-37.2819], Tmin=(1464.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCJCC=O) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]C(=O)C1(F)C=[C]C1F(14293)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-168.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2950,1000,180,800.674,801.124,801.285,801.336,801.356,801.388,801.462,801.498,802.008],'cm^-1')),
        HinderedRotor(inertia=(0.00390891,'amu*angstrom^2'), symmetry=1, barrier=(1.78057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90064,0.0525657,-4.9564e-05,3.06608e-08,-9.29154e-12,-20223.7,23.0245], Tmin=(100,'K'), Tmax=(737.665,'K')), NASAPolynomial(coeffs=[4.56214,0.0381331,-2.02149e-05,4.13529e-09,-3.01479e-13,-20616.4,10.9942], Tmin=(737.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-174.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (129.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-51.8741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (178.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (302.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-80.1709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (267.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-166.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (39.0381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (166.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-47.8651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-19.2001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (113.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (41.3756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (316.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (141.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (286.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-114.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (199.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (296.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (378.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (19.3784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-50.8389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (6.18556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (257.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (59.7935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['CO2(14)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(=O)[C](F)C1C=C1F(14277)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.6691e+10,'s^-1'), n=0.663, Ea=(162.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDe;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.7778e+12,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ-CdH;CO] for rate rule [cCs(-R!HR!H)CJ;CsJ-CdH;CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=O)C1=CC(F)(F)[CH]1(14278)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(266.863,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(OF)[C]1[CH]C(F)=C1(14279)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(181.237,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['O=C1O[CH]C(F)=C[C]1F(14280)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.75746e+12,'s^-1'), n=0.324012, Ea=(94.4534,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', 'O=[C]C1(F)[CH]C(F)=C1(14281)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['O=C1OC2C(F)=CC12F(14212)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.2e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 4.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['[O]C(=O)C1(F)C2[C](F)C21(14282)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['FC1=CC(F)([C]2OO2)[CH]1(14283)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(341.43,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 340.9 to 341.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['O=C1OC2[C](F)[CH]C12F(14223)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.75944e+08,'s^-1'), n=1.13751, Ea=(126.759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCd]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 4.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['O=C1OC2(F)[CH]C1(F)[CH]2(14238)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.34643e+11,'s^-1'), n=0.14784, Ea=(155.424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(F)C=C(F)C([O])=O(14284)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['[O]C1([O])C2C(F)=CC21F(14266)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHCd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(=O)C1=CC(F)=C1(14285)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(8.62981,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=O(517)', 'FC1=CC(F)=C1(307)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(306.062,'m^3/(mol*s)'), n=1.16366, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', '[O]C(=O)C1(F)C=C=C1(14286)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(11.5893,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CO2(14)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(120.712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=O(517)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[O]C(=O)C1=CC(F)=[C]1(14287)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(123.132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(=O)C1(F)[C]=C=C1(14288)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(219.202,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(=O)C1(F)[C]=C(F)C1(14289)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(O)C1(F)[C]C(F)=C1(14290)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    products = ['[O]C(=O)[C]1C=C(F)C1F(14291)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(180.81,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C(OF)C1(F)[CH][C]=C1(14292)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.1856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(=O)C1(F)C=[C]C1F(14293)'],
    products = ['[O]C(=O)C1(F)[CH]C(F)=C1(14207)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(147.46,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #4040',
    isomers = [
        '[O]C(=O)C1(F)[CH]C(F)=C1(14207)',
    ],
    reactants = [
        ('CO2(14)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4040',
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

